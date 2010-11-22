#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
EXTEND_LENGTH = 100

Dir.foreach("#{ALIGNMENTS_FOLDER}/") do |file|
  next unless file =~ /\.sorted.bam$/
  running_file          = running_file(file, "igvtools")
  tmp_file1             = "#{TMP_FOLDER}/#{file}.cov.tdf"
  tmp_file2             = "#{TMP_FOLDER}/#{file}.wig"
  output_file_1         = "#{IGVTOOLS_OUTPUT_FOLDER}/#{file}.cov.tdf"
  output_file_2         = "#{IGVTOOLS_OUTPUT_FOLDER}/#{file}.wig"
  
  input_path            = "#{ALIGNMENTS_FOLDER}/#{file}"
  next if File.exists? output_file_2 #Already processed
  next if File.exists? running_file #Being processed
  puts file
  `touch #{running_file}`
  begin
    `igvtools count -e #{EXTEND_LENGTH} "#{input_path}" "#{tmp_file1},#{tmp_file2}" "#{GENOMES_FOLDER}/#{GENOME}.genome"`
    FileUtils.mv(tmp_file1, output_file_1)
    FileUtils.mv(tmp_file2, output_file_2)
    
  rescue => e
    FileUtils.rm(output_file_1,     :force=>true)
    FileUtils.rm(output_file_2,     :force=>true)
    
    throw e
  ensure
    FileUtils.rm(tmp_file1,        :force=>true)
    FileUtils.rm(tmp_file2,        :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end