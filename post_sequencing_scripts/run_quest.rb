#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
require 'open3'
conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM analysis_pairs WHERE active=1 ORDER BY created_at asc")
res.each_hash do |row|
  f_name = row["foreground"]
  b_name = row["background"]
  f_user = row["foreground_user"]
  b_user = row["background_user"]
  
  f_path = "#{ALIGNMENTS_FOLDER}/#{f_user}/#{f_name}.sorted.bam"
  b_path = "#{ALIGNMENTS_FOLDER}/#{b_user}/#{b_name}.sorted.bam"
  analysis_folder_name = f_name + "_" + b_name
  analysis_folder_path = "#{QUEST_FOLDER}/#{f_user}/#{analysis_folder_name}"
  running_file        = running_file(analysis_folder_name, "run_quest")
  
  #Since quest and macs output folders with the same name, we must differentiate between them in the tmp folder.
  tmp_folder          = "#{TMP_FOLDER}/#{analysis_folder_name}_quest"
  
  next if File.exists? running_file #This is being processed
  next if File.exists? analysis_folder_path #this has already been analyzed.
  next unless File.exists? "#{f_path}"
  next unless File.exists? "#{b_path}"
  Dir.chdir(TMP_FOLDER)
  `touch #{running_file}`
  begin
    GENOME = "mm9"
    x = "generate_QuEST_parameters.pl -silent -bam_align_ChIP #{f_path} -bam_align_RX_noIP #{b_path} -gt #{QUEST_GENOME_TABLES}/#{GENOME} -ap #{tmp_folder}"
    puts x
    Open3.popen3(x) {
      |stdin, stdout, stderr|
      t = Thread.new(stderr) do |terr|
        while (line = terr.gets)
          puts "stderr: #{line}"
        end
      end
      t = Thread.new(stdout) do |terr|
         while (line = terr.gets)
           puts "stdout: #{line}"
         end
       end
      #stdin.puts "" #ONLY CALL THIS IF WE DO NOT HAVE A CONTROL.
      chip_config = 3
      if f_name =~ /PolII/
        chip_config = 2
      elsif f_name =~ /(H3)|(H4)/
        chip_config = 3
      end
      stdin.puts chip_config # options: 1:TF, 2:polII-like, 3:histone ChIP, 4:individual config. # in the future this must be parsed from the filename
      stdin.puts 2 #recommended peak calling params.
      stdin.puts "y" #to run QuEST analysis now.
      t.join()
    }
    `mkdir -p #{QUEST_FOLDER}/#{f_user}`
    `mv #{tmp_folder} #{analysis_folder_path}`
  ensure
    FileUtils.rm(tmp_folder,      :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end