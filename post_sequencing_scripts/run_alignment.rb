#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
samples_res = conn.query("SELECT * FROM sequencing_samples,sequencing_run where sequencing_run_id = sequencing_run.id and user != 'Control' and user != 'control' and type = 'chip'")

# Bowtie options
BT_NUM_THREADS	= 24

samples_res.each_hash do |sample|
  date          = sample["run_at"][0,10].gsub("-","_")
  lane          = sample["lane"].to_i
  user          = sample["user"].capitalize
  name          = sample["name"]
  type          = sample["type"].downcase
  run_id        = sample["sequencing_run_id"]
  genome        = GENOMES[sample["genome"]]
  
  base_file     = sample_filebase(run_id, date, lane, user, name)
  running_file  = running_file(base_file, "alignment")
  fastq_file    = "#{FASTQ_CHIP_FOLDER}/#{user}/#{base_file}_fastq.txt"
  base          = base_file + ".sorted" #.bam is added automatically
  tmp_file      = "#{TMP_FOLDER}/#{base}"
  output_file   = "#{ALIGNMENTS_FOLDER}/#{user}/#{base}.bam"
  
  next if File.exists? tmp_file
  next if File.exists? output_file
  next if File.exists? running_file
  
  puts fastq_file
  `touch #{running_file}`
  begin
    ## Do not align the last base because it has a higher error rate.
    bt_cmd        = "bowtie --chunkmbs 256 -p #{BT_NUM_THREADS} --best --strata -m 2 #{genome} --trim3 1 --sam \"#{fastq_file}\""
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
  
end