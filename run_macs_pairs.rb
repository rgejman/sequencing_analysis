#!/usr/bin/env ruby -KU
require 'constants'
require 'mysql'
conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM analysis_pairs WHERE active=1 ORDER BY created_at asc")
res.each_hash do |row|
  f_name = row["foreground"]
  b_name = row["background"]
  f_path = "#{ALIGNMENTS_FOLDER}/#{f_name}.sorted.bam"
  b_path = "#{ALIGNMENTS_FOLDER}/#{b_name}.sorted.bam"
  analysis_folder_name = f_name + "_" + b_name
  analysis_folder_path = "#{MACS_FOLDER}/#{analysis_folder_name}"
  next if File.exists? analysis_folder_path #this has already been analyzed.
  next unless File.exists? "#{f_path}"
  next unless File.exists? "#{b_path}"
  running_file        = running_file(analysis_folder_name, "run_macs_pair")
  output_folder       = "#{MACS_FOLDER}/#{analysis_folder_name}"
  complete_file       = "#{output_folder}/#{analysis_folder_name}_model.pdf"
  next if File.exists? complete_file # We use the model pdf as evidence that the run has completed.
  next if File.exists? running_file #This is being processed
  model_file          = "#{TMP_FOLDER}/#{analysis_folder_name}_model.r"
  
  extensions = ["summits.bed", "negative_peaks.xls", "peaks.xls", "peaks.bed", "model.r", "model.pdf"]
  
  puts analysis_folder_name
  Dir.chdir(TMP_FOLDER)
  `touch #{running_file}`
  `mkdir #{analysis_folder_path}`
  begin
    puts `macs14 -t #{f_path} -c #{b_path} --g mm -n #{analysis_folder_name}`
    `r --vanilla < #{model_file}`
    for ext in extensions
      `mv #{analysis_folder_name}_#{ext} #{output_folder}/`
    end
  rescue => e
    throw e
  ensure
    for ext in extensions
       FileUtils.rm("#{analysis_folder_name}_#{ext}", :force=>true)
     end
    FileUtils.rm(running_file,    :force=>true)
    conn.close
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end
conn.close