#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
Dir.foreach("#{ALIGNMENTS_FOLDER}/") do |file|
  running_file  = running_file(file, "make_bed")
  tmp_file      = "#{TMP_FOLDER}/#{file.gsub(".bam",".bed")}"
  output_file   = "#{ALIGNMENTS_FOLDER}/#{file.gsub(".bam",".bed")}"
  input_file    = "#{ALIGNMENTS_FOLDER}/#{file}"
  next unless file =~ /\.bam$/
  next if File.exists? output_file # The file has been processed in the past
  next if File.exists? running_file #This is being processed
  puts file
  `touch #{running_file}`
  Dir.chdir(TMP_FOLDER)
  begin
    `bamToBed -i #{input_file} > #{tmp_file}`
    FileUtils.mv(tmp_file, output_file)
  rescue => e
    FileUtils.rm(output_file,     :force=>true)
    raise e
  ensure
    FileUtils.rm(tmp_file,        :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end
