#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'

# Bowtie options
BT_NUM_THREADS		      = 10

files = Dir.glob("**/*_fastq.txt").collect{|p| File.join(FASTQ_CHIP_FOLDER,p)}
for file in files
  base                    = file.split("/").last.gsub('_fastq.txt','.sam')
  running_file            = running_file(base, "alignment")
  input_file              = file
  tmp_file                = "#{TMP_FOLDER}/#{base}"
  user                    = base.split("_").first
  output_file             = "#{ALIGNMENTS_FOLDER}/#{user}/#{base}"

  puts tmp_file
  next if File.exists? tmp_file
  put output_file
  next if File.exists? output_file
  next if File.exists? running_file

  puts file
  `touch #{running_file}`
  begin
    ## Do not align the last base because it has a higher error rate.
    `bowtie --chunkmbs 128 -p #{BT_NUM_THREADS} --best -m 2 #{GENOME} --trim3 1 --sam "#{input_file}" "#{tmp_file}"`
    FileUtils.mv(tmp_file, output_file)
    FileUtils.rm(input_file) # Removes the unsorted SAM file.
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