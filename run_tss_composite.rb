#!/usr/bin/env ruby -KU
require 'constants'
require 'mysql'

def count_scores(scores, file)
  Open3.popen3(`gunzip -c #{file}`) { |stdin, stdout, stderr|
    t = Thread.new(stdout) do |terr|
      while (line = terr.gets)
        puts "stderr: #{line}"
      end
    end
    while (line = stderr.gets) #for some reason gunzip output to stderr via Open3.popen3
      if line[0,5] == "track" #the 1st header line starts with "track"
        h2_line = stdin.gets #variableStep header line
        chr = h2_line.split(" ")[1].split("=")[1] #the chromosome #.
        next
      end          
      tokens = line.split(" ") #0 = pos, 1 = score
      pos = tokens[0].to_i
      for coord_pair in TSS_COORDS[chr]
        if coord_pair[0] <= pos and coord_pair[1] >= pos
          scores[coord_pair[1] - pos - 2000] += tokens[1]  # get the pos relative to the coord_pairs (0-based) and add the score to this pos.
        end
      end
    end
    t.join()
  }
end


conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM analysis_pairs WHERE active=1 ORDER BY created_at asc")
res.each_hash do |row|
  f_name = row["foreground"]
  b_name = row["background"]
  analysis_folder_name = f_name + "_" + b_name
  analysis_folder_path = "#{QUEST_FOLDER}/#{analysis_folder_name}"
  next unless File.exists? analysis_folder_path #has not yet been analyzed.
  composite_plot_path = "#{COMPOSITE_PLOTS_FOLDER}/#{analysis_folder_name}_tss.png"
  next if File.exists? composite_plot_path # the composite plot has already been generated, so this is done.
  running_file        = running_file(analysis_folder_name, "make_tss_composite")
  tmp_folder          = "#{TMP_FOLDER}/#{analysis_folder_name}"
  next if File.exists? running_file #This is being processed
  # Use normalized profiles to compare the two!
  f_wig_path = "#{analysis_folder_path}/tracks/wig_profiles/ChIP_normalized.profile.wig.gz"
  b_wig_path = "#{analysis_folder_path}/tracks/wig_profiles/background_normalized.profile.wig.gz"
  next unless File.exists? f_wig_path
  next unless File.exists? b_wig_path
  begin
    Dir.chdir(TMP_FOLDER)
    `touch #{running_file}`
    TSS_COORDS = {}
    # Put the TSS coordinates into a data structure (array of start/end pairs in hashmap keyed on chromosome)
    File.readlines("#{USEFUL_BED_FILES}/mm9.tss.2kb.bed").each {|line|
      tokens = line.split("\t")
      TSS_COORDS[tokens[0]] ||= []
      TSS_COORDS[tokens[0]] << [tokens[1].to_i, tokens[2].to_i] # i.e. TSS_COORDS[chr] << [start, end]
    }
    pp TSS_COORDS
    TSS_SCORES_F = Array.new(2000, 0) #Initialized w/ 2000 "0" objects
    TSS_SCORES_B = Array.new(2000, 0) #Initialized w/ 2000 "0" objects
  
    t1 = Thread.new(TSS_SCORES_F, f_wig_path) do |tss_scores, file|
      count_scores(tss_scores, file)
    end
    t2 = Thread.new(TSS_SCORES_B, b_wig_path) do |tss_scores, file|
      count_scores(tss_scores, file)
    end
    t1.join()
    t2.join()
    pp TSS_SCORES_F
    pp TSS_SCORES_B
    
    #`mv #{tmp_folder} #{analysis_folder_path}`
  ensure
    #FileUtils.rm(tmp_folder,      :force=>true)
    FileUtils.rm(running_file,    :force=>true)
    t1.kill()
    t2.kill()
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end