#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'

# Bowtie options
BT_NUM_THREADS		      = 24

files = Dir.glob("#{FASTQ_CHIP_FOLDER}/**/*_fastq.txt")
for file in files
  base                    = file.split("/").last.gsub('_fastq.txt','.sorted') #.bam is added automatically
  running_file            = running_file(base, "alignment")
  input_file              = file
  tmp_file                = "#{TMP_FOLDER}/#{base}"
  user                    = base.split("_").first
  output_file             = "#{ALIGNMENTS_FOLDER}/#{user}/#{base}.bam"

  next if File.exists? tmp_file
  next if File.exists? output_file
  next if File.exists? running_file

  puts file
  `touch #{running_file}`
  begin
    ## Do not align the last base because it has a higher error rate.
    bt_cmd        = "bowtie --chunkmbs 256 -p #{BT_NUM_THREADS} --best -m 2 #{GENOME} --trim3 1 --sam \"#{input_file}\""
    convert_bam   = "samtools view -hbSu -"
    sort_bam      = "samtools sort - #{tmp_file}"
    `#{bt_cmd} | #{convert_bam} | #{sort_bam}`
    `mkdir -p #{ALIGNMENTS_FOLDER}/#{user}`
    FileUtils.mv(tmp_file + ".bam", output_file)
    puts `samtools index #{output_file}`
  rescue => e
    FileUtils.rm(output_file,     :force=>true)
    throw e
  ensure
    FileUtils.rm(tmp_file,        :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  #break #since this script is memory and CPU intensive, we just let it loop through all the files needing conversion
end

# bowtie options:
# -p = # threads
# -n = max # of mistmatches in seed
# --best = hits guaranteed best stratum; ties broken by quality
# -m = Supress all alignments if < <int> exist
# -o = override offrate of index; must be >= index's offrate (prevents HD
# 	thrashing.)