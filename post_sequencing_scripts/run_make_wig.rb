#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'

Dir.foreach("#{ALIGNMENTS_FOLDER}/") do |file|
  next unless file =~ /\.sorted\.bam$/
  running_file        = running_file(file, "make_wig")
  tmp_wig_file        = "#{TMP_FOLDER}/#{file.gsub(".bam",".wig")}"
  tmp_bigwig_file     = "#{TMP_FOLDER}/#{file.gsub(".bam",".bigwig")}"
  output_wig_file     = "#{WIG_FOLDER}/#{file.gsub(".bam",".wig")}" ### Must be changed to temporary outfile.
  
  input_file		= "#{ALIGNMENTS_FOLDER}/#{file}"
  next if File.exists? output_wig_file # The file has been processed in the past
  next if File.exists? running_file #This is being processed
  puts file
  Dir.chdir(TMP_FOLDER)
  `touch #{running_file}`
  begin
    puts `bam_to_wiggle.py #{input_file} -w #{tmp_wig_file}`
    FileUtils.mv(tmp_wig_file, output_wig_file)
  rescue => e
    FileUtils.rm(output_wig_file,         :force=>true)
    throw e
  ensure
    FileUtils.rm(tmp_wig_file,    :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end
