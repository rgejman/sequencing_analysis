#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
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
  analysis_folder_path = "#{MACS_FOLDER}/#{f_user}/#{analysis_folder_name}"
  running_file        = running_file(analysis_folder_name, "run_macs")
  puts "macs: Checking #{analysis_folder_name}"
  #Since quest and macs output folders with the same name, we must differentiate between them in the tmp folder.
  tmp_folder          = "#{TMP_FOLDER}/#{analysis_folder_name}_macs"
  
  running_file        = running_file(analysis_folder_name, "run_macs_pair")
  complete_file       = "#{analysis_folder_path}/#{analysis_folder_name}_model.pdf"
  model_file          = "#{analysis_folder_path}/#{analysis_folder_name}_model.r"
  
  puts "\t Running file exists?"
  next if File.exists? running_file #This is being processed
  puts "\t Analysis folder exists?"
  next if File.exists? analysis_folder_path #this has already been analyzed.
  puts "\t foreground exists?"
  next unless File.exists? "#{f_path}"
  puts "\t background exists?"
  next unless File.exists? "#{b_path}"
  puts "\t Model exists?"
  next if File.exists? complete_file # We use the model pdf as evidence that the run has completed.
    
  puts analysis_folder_name
  `touch #{running_file}`
  `mkdir -p #{analysis_folder_path}`
  Dir.chdir(analysis_folder_path)
  begin
    puts `macs14 -t #{f_path} -c #{b_path} --g mm -n #{analysis_folder_name}`
    `r --vanilla < #{model_file}`
    `intersectBed -u -wa -a #{USEFUL_BED_FILES}/mm9.ensembl_with_symbols.tss.2kb.prot_coding.bed -b #{analysis_folder_name}_peaks.bed > #{analysis_folder_name}_peaks.intersect_tss.bed`
  rescue => e
    throw e
  ensure
    FileUtils.rm(running_file,    :force=>true)
    conn.close
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end
conn.close