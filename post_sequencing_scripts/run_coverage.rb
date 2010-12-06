#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
EXTEND_LENGTH = 100
files = Dir.glob("#{ALIGNMENTS_FOLDER}/**/*.bam")

for file in files
  running_file          = running_file(file, "coverage")
  base                  = file.split("/").last.gsub('.bam','.cov.tdf') #.bam is added automatically
  user                    = base.split("_").first
  tmp_file1             = "#{TMP_FOLDER}/#{base}.cov.tdf"
  tmp_file2             = "#{TMP_FOLDER}/#{base}.wig"
  output_file_1         = "#{COVERAGE_FOLDER}/#{user}/#{base}.cov.tdf"
  output_file_2         = "#{WIG_FOLDER}/#{user}/#{base}.wig"
  
  input_path            = file
  
  next if File.exists? output_file_1 #Already processed
  next if File.exists? output_file_2 #Already processed
  
  next if File.exists? running_file #Being processed
  puts file
  `touch #{running_file}`
  begin
    `igvtools count -e #{EXTEND_LENGTH} "#{input_path}" "#{tmp_file1},#{tmp_file2}" "#{GENOMES_FOLDER}/#{GENOME}.genome"`
    `mkdir -p #{COVERAGE_FOLDER}/#{user}`
    `mkdir -p #{WIG_FOLDER}/#{user}`
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