#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
files = Dir.glob("#{ALIGNMENTS_FOLDER}/**/*.sorted.bam")

for file_path in files
  next unless file_path =~ /\.sorted.bam$/
  analysis_folder_name  = file.split("/").last.gsub(".sorted.bam","")
  analysis_folder_path  = "#{MACS_FOLDER}/#{analysis_folder_name}"
  next if File.exists? analysis_folder_path #this has already been analyzed.
  next unless File.exists? "#{file_path}"
  running_file        = running_file(analysis_folder_name, "run_macs_single")
  output_folder       = "#{MACS_FOLDER}/#{analysis_folder_name}"
  complete_file       = "#{output_folder}/#{analysis_folder_name}_model.pdf"
  next if File.exists? complete_file # We use the model pdf as evidence that the run has completed.
  next if File.exists? running_file #This is being processed
  model_file          = "#{MACS_FOLDER}/#{analysis_folder_name}_model.r"
    
  puts analysis_folder_name
  `touch #{running_file}`
  `mkdir #{analysis_folder_path}`
  Dir.chdir(analysis_folder_path)
  begin
    puts `macs14 -t #{file_path} --g mm -n #{analysis_folder_name}`
    `r --vanilla < #{model_file}`
    `intersectBed -u -wa -a #{USEFUL_BED_FILES}/mm9.ensembl_with_symbols.tss.2kb.prot_coding.bed -b #{analysis_folder_name}_peaks.bed > #{analysis_folder_name}_peaks.overlap_tss.bed`
  rescue => e
    throw e
  ensure
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end