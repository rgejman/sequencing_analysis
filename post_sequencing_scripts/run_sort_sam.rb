#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'

Dir.foreach("#{ALIGNMENTS_FOLDER}/") do |file|  
  running_file = running_file(file, "sort_sam")
  tmp_file = "#{TMP_FOLDER}/#{file.gsub(".sam",".sorted.sam")}"
  output_file = "#{ALIGNMENTS_FOLDER}/#{file.gsub(".sam",".sorted.sam")}"
  input_path  = "#{ALIGNMENTS_FOLDER}/#{file}"
  next unless file =~ /\.sam$/ and file !~ /\.sorted\./
  next if File.exists? output_file #Already processed
  next if File.exists? running_file #Being processed
  puts file
  begin
    `touch #{running_file}`
    `igvtools sort  "#{input_path}" "#{tmp_file}"`
    ### Let's check the filesize of the sorted and unsorted to make sure that the sort worked correctly
    s1 = File.size(input_path)
    s2 = File.size(tmp_file)
    if Math.abs(s1 - s2) > 1048576 #1 megabyte in bytes
      throw "Original file was #{s1} and sorted file is #{s2}. Sort failed. Removing sorted file."
    end
    FileUtils.mv(tmp_file, output_file)
    FileUtils.rm(input_path) # Removes the unsorted SAM file.
  rescue => e
    FileUtils.rm(output_file,     :force=>true)
    throw e
  ensure
    FileUtils.rm(tmp_file, :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end