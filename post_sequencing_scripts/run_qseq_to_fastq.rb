#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'

Dir.foreach("#{QSEQ_FOLDER}/") do |file|
  next unless file =~ /\_qseq.txt$/
  output_filename         = file.gsub("qseq","fastq")
  running_file            = running_file(file, "qseq_to_fastq")
  output_file             = "#{TMP_FOLDER}/#{output_filename}" #We manually move this.

  next if File.exists? running_file
  next if File.exists? output_file
  
  base = file.gsub("_qseq.txt","")
  
  puts file
  `touch #{running_file}`

  Dir.chdir(QSEQ_FOLDER)
  begin
    cat  = "cat #{file}"
    awk1 = "gawk '{gsub(/\./, \"N\", $9);print}'"
    awk2 = "awk '{print \"@\"$1\":#{base}:\" $11 \":\" $3 \":\" $4 \":\" $5 \":\" $6\"#0/1\"; print $9; print \"+\"$1\":#{base}:\" $11 \":\" $3 \":\" $4 \":\" $5 \":\" $6\"#0/1\" ; print $10}'"
    out  = "> #{output_file}.txt"
    `#{cat} | #{awk1} | #{awk2} #{out}`
  rescue => e
    FileUtils.rm(output_file,     :force=>true)
    throw e
  ensure
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end