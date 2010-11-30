#!/usr/bin/env ruby -KU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
require 'constants'
require 'mysql'
NUM_THREADS = 24

conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
res = conn.query("SELECT * FROM sequencing_run WHERE run_at IS NOT NULL and illumina_run_id IS NOT NULL ORDER BY run_at desc")
res.each_hash do |sequencing_run|
  run_id                        = sequencing_run["illumina_run_id"]
  running_file                  = running_file(run_id, "BclToQseq")
  next if File.exists? running_file
  intensities_folder            = "#{RAW_FOLDER}/#{run_id}/Data/Intensities"
  base_calls_folder             = "#{intensities_folder}/BaseCalls"
  output_folder                 = base_calls_folder

  `touch running_file`
  begin
    cmd = "setupBclToQseq.py -i #{base_calls_folder} -p #{intensities_folder} -o #{output_folder} --in-place --overwrite"
    puts cmd
    #`#{cmd}`
    Dir.chdir(output_folder)
    cmd = "make -j #{NUM_THREADS}"
    puts cmd
    #`cmd`
  rescue => e
    throw e
  ensure
    FileUtils.rm(running_file,    :force=>true)
    conn.close
  end
  break
end