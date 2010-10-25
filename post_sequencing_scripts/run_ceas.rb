$: << File.expand_path(File.dirname(__FILE__) + "/../")

require 'constants'

## ASK SCOTT: Does CEAS not extend the input alignments, like igvtools can do?

Dir.foreach("#{ALIGNMENTS_FOLDER}/") do |file|
  next unless file =~ /\.wig$/
  running_file      = running_file(file, "ceas")
  tmp_output_file		= "#{CEAS_FOLDER}/"
  input_file		    = "#{ALIGNMENTS_FOLDERS}/#{file}"
  basename          = file.gsub(".bed","")
  next if !File.exists? wig_file  #there must be a corresponding wig file if we want to run ceas.
  next if File.exists? output_file # The file has been processed in the past
  next if File.exists? running_file #This is being processed
  Dir.chdir(TMP_FOLDER)
  puts file
  `touch #{running_file}`
  begin
    #Do average signaling profile
    `ceas -w #{wig_file} -g #{CEAS_ANNOTATION_TABLES}/mm9/refGene`
  rescue => e
    FileUtils.rm(output_file,     :force=>true)
    throw e
  ensure
    FileUtils.rm(tmp_file,        :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end


# Ceas has 3 different modes:
# => 1.ChIP region annotation: "CEAS estimates the relative enrichment level of ChIP regions in each gene feature with respect to the whole genome."
# => 2.Gene-centered annotation: "To this end, CEAS divides every gene into three equal fractions and, for each fraction, calculates the percentage of the area covered by ChIP regions. It also estimates the percentages of the promoter and downstream of the gene (3kb upstream of TSS and downstream of TTS by default) that are covered by ChIP regions."
# => 3.Average signal profiling within/near important genomic features: uses only ceas [options] -g gdb -w wig