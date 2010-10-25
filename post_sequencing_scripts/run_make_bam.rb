require 'constants'
Dir.foreach("#{ALIGNMENTS_FOLDER}/") do |file|
  next unless file =~ /\.sorted\.sam$/
  running_file  = running_file(file, "make_bam")
  tmp_file      = "#{TMP_FOLDER}/#{file.gsub(".sam",".bam")}"
  output_file		= "#{ALIGNMENTS_FOLDER}/#{file.gsub(".sam",".bam")}"
  input_file		= "#{ALIGNMENTS_FOLDER}/#{file}"
  next if File.exists? output_file # The file has been processed in the past
  next if File.exists? running_file #This is being processed
  puts file
  `touch #{running_file}`
  begin
    `samtools view -S -b -o #{tmp_file} #{input_file}`
    FileUtils.mv(tmp_file, output_file)
  rescue => e
    FileUtils.rm(output_file,     :force=>true)
    throw e
  ensure
    FileUtils.rm(tmp_file,        :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end
