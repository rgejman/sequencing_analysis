$: << File.expand_path(File.dirname(__FILE__) + "../")
require 'constants'
EXTEND_LENGTH = 100

Dir.foreach("#{ALIGNMENTS_FOLDER}/") do |file|
  next unless file =~ /\.sorted.bam$/
  running_file          = running_file(file, "igvtools")
  tmp_file              = "#{TMP_FOLDER}/#{file}.cov.tdf"
  output_file           = "#{IGVTOOLS_OUTPUT_FOLDER}/#{file}.cov.tdf"
  input_path            = "#{ALIGNMENTS_FOLDER}/#{file}"
  next if File.exists? output_file #Already processed
  next if File.exists? running_file #Being processed
  puts file
  `touch #{running_file}`
  begin
    `igvtools count -e #{EXTEND_LENGTH} "#{input_path}" "#{tmp_file}" "#{GENOMES_FOLDER}/#{GENOME}.genome"`
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