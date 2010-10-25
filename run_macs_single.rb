#!/usr/bin/env ruby -KU
require 'constants'
require 'mysql'
Dir.foreach("#{ALIGNMENTS_FOLDER}/") do |file|
  next unless file =~ /\.sorted.bam$/
  file_path             = "#{ALIGNMENTS_FOLDER}/#{file}"
  analysis_folder_name  = file.gsub(".sorted.bam","")
  analysis_folder_path  = "#{MACS_FOLDER}/#{analysis_folder_name}"
  next if File.exists? analysis_folder_path #this has already been analyzed.
  next unless File.exists? "#{file_path}"
  running_file        = running_file(analysis_folder_name, "run_macs_single")
  output_folder       = "#{MACS_FOLDER}/#{analysis_folder_name}"
  complete_file       = "#{output_folder}/#{analysis_folder_name}_model.pdf"
  next if File.exists? complete_file # We use the model pdf as evidence that the run has completed.
  next if File.exists? running_file #This is being processed
  model_file          = "#{TMP_FOLDER}/#{analysis_folder_name}_model.r"
  
  extensions = ["summits.bed", "negative_peaks.xls", "peaks.xls", "peaks.bed", "model.r", "model.pdf"]
  
  puts analysis_folder_name
  `touch #{running_file}`
  `mkdir #{analysis_folder_path}`
  Dir.chdir(TMP_FOLDER)
  begin
    puts `macs14 -t #{file_path} --g mm -n #{analysis_folder_name}`
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