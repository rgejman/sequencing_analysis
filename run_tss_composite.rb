#!/usr/bin/env ruby -KU
require 'constants'
require 'mysql'

def count_scores(tss_coords_file, folder, output_file)
  tss_coords = {}
  File.readlines(tss_coords_file).each {|line|
    tokens = line.split("\t")
    tss_coords[tokens[0]] ||= []
    tss_coords[tokens[0]] << [tokens[1].to_i, tokens[2].to_i, tokens[5]] # i.e. TSS_COORDS[chr] << [start, end, strand]
  }
  scores = Array.new(2001, 0.0) #including the "0" position there are 2000 positions.
  n = 0
  Dir.chdir(folder)
  for file in Dir["chr*.wig.gz"]
    chr = file.gsub(".wig.gz","")
    lines = `gunzip -c #{file}`
    for line in lines
      n+=1
      t = line.split(" ") #0 = pos, 1 = score
      pos = t[0].to_i
      next if tss_coords[chr].empty? # skip to the next line if we are passed any genes in the TSS file.
      need_to_delete = false
      for l,r,s in tss_coords[chr]
        if r < pos # no need to keep looking through the TSS if we have passed the pos.
          need_to_delete = true
          break
        end
        if l <= pos and r >= pos
          if s == "+" #strand
            scores[2000 - (r - pos)] += t[1].to_f  # get the pos relative to the coord_pairs (0-based) and add the score to this pos.
          else
            scores[r - pos] += t[1].to_f
          end
          break
        end
      end
       # this is much faster than using delete_if because we don't have to iterate through every entry in tss_coords
       # and we only delete if we know there is stuff to delete.
      tss_coords[chr].delete_if{|a| a[1] < pos } if need_to_delete
      if n % 100 == 0
        File.open(output_file, "w") do |f|
          for score in scores
            f.puts score
          end
        end
      end
    end
  end
  File.open(output_file, "w") do |f|
    for score in scores
      f.puts score
    end
  end
end


conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM analysis_pairs WHERE active=1 ORDER BY created_at asc")
res.each_hash do |row|
  f_name                      = row["foreground"]
  b_name                      = row["background"]
  analysis_folder_name        = f_name + "_" + b_name
  quest_analysis_folder_path  = "#{QUEST_FOLDER}/#{analysis_folder_name}"
  composite_plot_path         = "#{COMPOSITE_PLOTS_FOLDER}/#{analysis_folder_name}/tss.png"
  running_file                = running_file(analysis_folder_name, "make_tss_composite")
  tmp_folder                  = "#{TMP_FOLDER}/#{analysis_folder_name}"
  final_folder_path           = "#{COMPOSITE_PLOTS_FOLDER}/#{analysis_folder_name}"
  f_wig_path                  = "#{quest_analysis_folder_path}/tracks/wig_profiles/by_chr/ChIP_normalized"
  b_wig_path                  = "#{quest_analysis_folder_path}/tracks/wig_profiles/by_chr/background_normalized"
  tss_coords_file             = "#{USEFUL_BED_FILES}/mm9.tss.2kb.bed"
  next if File.exists? composite_plot_path # the composite plot has already been generated, so this is done.
  next if File.exists? running_file #This is being processed
  next unless File.exists? f_wig_path
  next unless File.exists? b_wig_path

  Dir.chdir(TMP_FOLDER)
  `mkdir -p #{tmp_folder}`
  `touch #{running_file}`

  begin
    # Put the TSS coordinates into a data structure (array of start/end pairs in hashmap keyed on chromosome)
    child1 = fork
    count_scores(tss_coords_file, f_wig_path, "#{tmp_folder}/scores_f.txt") if child1.nil? # child1 is nil if the thread is the child.
    child2 = fork unless child1.nil? # fork if we are the parent.
    count_scores(tss_coords_file, b_wig_path, "#{tmp_folder}/scores_b.txt") if child2.nil? # child2 is nil if the thread is the 2nd fork.
    Process.waitall
    exit if child1.nil? or child2.nil? #if you are either one of the children, exit here.
    #control script continues here.

    #Dir.chdir(tmp_folder)


    #{}`mv #{tmp_folder} #{COMPOSITE_PLOTS_FOLDER}/`
  ensure
    Process.kill child1 unless child1.nil?
    Process.kill child2 unless child2.nil?
    FileUtils.rm(tmp_folder,      :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end