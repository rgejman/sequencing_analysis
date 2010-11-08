#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
files = Dir.entries("#{FASTQ_CHIP_FOLDER}/") + Dir.entries("#{FASTQ_RNA_SEQ_FOLDER}/")
for file in files
  next unless file =~ /\.txt$/
  running_file = running_file(file, "fastqc")
  fastqc_output_folder_name = "#{file.gsub('.txt','')}_fastqc"
  fastqc_output_folder_path = "#{FASTQC_FOLDER}/#{fastqc_output_folder_name}"

  next if File.exists? fastqc_output_folder_path # The file has been processed in the past
  next if File.exists? running_file #This is being processed

  puts file
  Dir.chdir(TMP_FOLDER)
  `touch #{running_file}`
  begin
    `fastqc #{FASTQ_FOLDER}/#{file}`
    `mv #{TMP_FOLDER}/#{fastqc_output_folder_name} #{FASTQC_FOLDER}/`
  rescue => e
    FileUtils.rm(fastqc_output_folder_path,     :force=>true)
    throw e
  ensure
    FileUtils.rm("#{TMP_FOLDER}/#{fastqc_output_folder_name}.zip",  :force=>true)
    FileUtils.rm("#{TMP_FOLDER}/#{fastqc_output_folder_name}",      :force=>true)
    FileUtils.rm(running_file,                                      :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end