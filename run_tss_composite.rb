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
  #n = 0
  Dir.chdir(folder)
  child_id = nil
  files = Dir["chr*.wig.gz"]
  pipes = []
  for file in files
    rd, wr = IO.pipe # only the parent should do this.
    child_id = fork #fork the process
    if !child_id.nil? #if we are the parent
      pipes << [rd,wr]
    else
      chr = file.gsub(".wig.gz","")
      next if tss_coords[chr].nil? or tss_coords[chr].empty?
      lines = `gunzip -c #{file}`
      for line in lines
        #n+=1
        t = line.split(" ") #0 = pos, 1 = score
        pos = t[0].to_i
        break if tss_coords[chr].empty? # ignore the rest of the file if we are passed any genes in the TSS file.
        for l,r,s in tss_coords[chr]
          if r < pos # no need to keep looking through the TSS if we have passed the pos AND we have TSSes to delete.
            tss_coords[chr].delete_if{|a| a[1] < pos }
            break
          end
          if l <= pos and r >= pos
            score = t[1].to_f
            if s == "+" #strand
              scores[2000 - (r - pos)] += score  # get the pos relative to the coord_pairs (0-based) and add the score to this pos.
            else
              scores[r - pos] += score
            end
            #don't break here because there may be multiple TSS for which this pos matches.
          end
        end
      end
      if child_id.nil? # if we are the child
        rd.close
        data = Marshal.dump(scores)
        wr.print [data.length].pack('l')
        wr.print data
        wr.close
        exit(0)
      end
    end
  end
  #only parent should reach here
  Process.waitall()
  for rd,wr in pipes
    wr.close
    data = rd.read.unpack('l').shift
    data = rd.read( data )
    loc_scores = Marshal.load(data)
    rd.close
    for i in 0...scores
      scores[i] += loc_scores[i]
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
    exit(0) if child1.nil? or child2.nil? #if you are either one of the children, exit here.
    #control script continues here.
    a = Process.waitall()
    raise "Forked process failed." if a.any?{|uid, ps| !ps.success? }

    Dir.chdir(tmp_folder)
    `r --vanilla < #{SCRIPTS_FOLDER}/make_composite_tss_plot.r`

    `mv -f #{tmp_folder} #{COMPOSITE_PLOTS_FOLDER}/`
  ensure
    FileUtils.rm(tmp_folder,      :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end