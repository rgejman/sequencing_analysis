#!/usr/bin/env ruby -wKU
$: << File.expand_path(File.dirname(__FILE__) + "/../")
def alive?(pid)
  pid = Integer("#{ pid }")
  begin
    Process::kill 0, pid
    true
  rescue Errno::ESRCH
    false
  end
end

require 'constants'
files = Dir.glob("#{ALIGNMENTS_FOLDER}/**/*.bam")
forks = []
for path in files
  forks << fork do
    name                          = File.basename(path).gsub(".bam","")
    user                          = name.split("_")[0]
    running_file                  = running_file(name, "fastqc")
    fastqc_tmp_folder_path        = path.gsub(".bam","") + "c"
    fastqc_tmp_zip_path           = "#{path}c.zip"
    fastqc_all_output_folder_path         = "#{FASTQC_FOLDER}/#{user}/#{name}_bamqc"
    fastqc_aligned_output_folder_path     = "#{FASTQC_FOLDER}/#{user}/#{name}_aligned_bamqc"

    next if File.exists? fastqc_output_folder_path # The file has been processed in the past
    next if File.exists? running_file #This is being processed
    puts name
    `touch #{running_file}`
    begin
      cmd_all     = "fastqc -f bam #{path}"
      cmd_aligned = "fastqc -f bam_mapped #{path}"
      `mkdir -p #{FASTQC_FOLDER}/#{user}/`
      
      `#{cmd_all}}`
      FileUtils.mv(fastqc_tmp_folder_path,fastqc_all_output_folder_path)
      `#{cmd_aligned}`
      FileUtils.mv(fastqc_tmp_folder_path,fastqc_aligned_output_folder_path)
    rescue => e
      FileUtils.rm(fastqc_output_folder_path,             :force=>true)
      FileUtils.rm(fastqc_aligned_output_folder_path,     :force=>true)
      FileUtils.rm(fastqc_tmp_zip_path,                   :force=>true)
      FileUtils.rm(fastqc_tmp_folder_path,                :force=>true)
      throw e
    ensure
      FileUtils.rm(fastqc_tmp_folder_path,        :force=>true)
      FileUtils.rm(fastqc_tmp_zip_path,           :force=>true)
      FileUtils.rm(running_file,                    :force=>true)
    end
  end
  while forks.length >= 20
    puts "#{forks.length} running. 20 max."
    Process.wait(forks[0])
    forks.delete_at(0)
  end
  
end
for fork in forks
  puts "P: Waiting for child (#{fork})"
  Process.wait(fork)
end