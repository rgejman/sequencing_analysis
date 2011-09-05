#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'

BT_NUM_THREADS	= 18

def run_rmdup(user, base_file)
  running_file  = running_file(base_file, "rmdup")
  bam_in        = "#{ALIGNMENTS_FOLDER}/#{user}/#{base_file}.sorted.bam"
  bam_out       = "#{ALIGNMENTS_FOLDER}/#{user}/#{base_file}.rmdup.sorted.bam"
  tmp_file      = "#{TMP_FOLDER}/#{base_file}.rmdup.sorted.bam"
  
  return unless File.exists? bam_in
  return if File.exists? tmp_file
  return if File.exists? bam_out
  return if File.exists? running_file
  
  puts "#{bam_in}--removing duplicate reads"
  begin
    ## Do not align the last base because it has a higher error rate.
    cmd        = "samtools rmdup #{bam_in} #{tmp_file}"
    puts cmd
    `#{cmd}`
    FileUtils.mv(tmp_file, bam_out)
    puts `samtools index #{bam_out}`
  rescue => e
    FileUtils.rm(bam_out,     :force=>true)
    throw e
  ensure
    FileUtils.rm(tmp_file,        :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
end

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
samples_res = conn.query("SELECT * FROM sequencing_samples,sequencing_run where sequencing_run_id = sequencing_run.id and user != 'Control' and user != 'control' and type = 'chip' and rmdup=1")

# Run through the files in the DB

samples_res.each_hash do |sample|
  date          = sample["run_at"][0,10].gsub("-","_")
  lane          = sample["lane"].to_i
  user          = sample["user"].capitalize
  name          = sample["name"]
  type          = sample["type"].downcase
  run_id        = sample["sequencing_run_id"]
  genome        = GENOMES[sample["genome"]]
  
  base_file     = sample_filebase(run_id, date, lane, user, name)
  next if sample["post_process"].to_i == 0
  run_rmdup(user, base_file)
end