$: << File.expand_path(File.dirname(__FILE__) + "../")
require 'constants'
Dir.foreach("#{ALIGNMENTS_FOLDER}/") do |file|
  next unless file =~ /\.sorted\.bam$/
  
  running_file        = running_file(file, "bam_index")
  input_file		      = "#{ALIGNMENTS_FOLDER}/#{file}"
  output_file         = input_file + ".bai"
  next if File.exists? output_file  # The file has been processed in the past
  next if File.exists? running_file # This is being processed
  puts file
  
  Dir.chdir(TMP_FOLDER)
  `touch #{running_file}`
  begin
    puts `samtools index #{input_file}`
  rescue => e
    raise
  ensure
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end