#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'

# Bowtie options
BT_NUM_THREADS		      = 10
BT_OVERRIDE_OFFRATE	    = 9
if ARGV.length > 0
  files = ARGV
else
  files = 
end

Dir.foreach("#{FASTQ_FOLDER}/") do |file|
  next unless file =~ /\.txt$/
  base                    = file.gsub('.txt','.sam')
  running_file            = running_file(base, "alignment")
  tmp_file                = "#{TMP_FOLDER}/#{base}"
  output_file             = "#{ALIGNMENTS_FOLDER}/#{base}"

  next if File.exists? tmp_file
  next if File.exists? output_file
  next if File.exists? running_file

  puts file
  `touch #{running_file}`
  begin
    `bowtie -p #{BT_NUM_THREADS} --best -m 1 #{GENOME} -o #{BT_OVERRIDE_OFFRATE} --sam "#{FASTQ_FOLDER}/#{file}" "#{tmp_file}"`
    FileUtils.mv(tmp_file, output_file)
  rescue => e
    FileUtils.rm(output_file,     :force=>true)
    throw e
  ensure
    FileUtils.rm(tmp_file, :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  break # We break so that other scripts have a chance to execute before we try this one again.
end

# bowtie options:
# -p = # threads
# -n = max # of mistmatches in seed
# --best = hits guaranteed best stratum; ties broken by quality
# -m = Supress all alignments if < <int> exist
# -o = override offrate of index; must be >= index's offrate (prevents HD
# 	thrashing.)