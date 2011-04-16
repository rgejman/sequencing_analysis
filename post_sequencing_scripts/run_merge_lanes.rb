conn = Mysql::new(MYSQL_HOST, MYSQL_USER, MYSQL_PASS, MYSQL_DB)
rows = conn.query("SELECT * FROM lanes_to_merge")

# Run through the files in the DB

rows.each_hash do |row|
  merged_name     = row["merged_name"]
  user            = row["user"]
  merged_filename = merged_name + ".sorted.bam"
  merged_filepath = "#{ALIGNMENTS_FOLDER}/#{user}/#{merged_filename}"
  tmp_file        = "#{TMP_FOLDER}/#{merged_filename}.sorted.bam"
  running_file    = running_file(merged_name, "merge_lanes")
  
  next if File.exists? merged_filepath # Already been merged
  next if File.exists? running_file #Being processed
  
  lanes = []
  missing_files = false
  for l in 1..10
    lane = row["lane_#{l}"]
    next if lane.nil? || lane.strip == ""
    lane = "#{ALIGNMENTS_FOLDER}/#{user}/#{lane}"
    lanes << lane
    unless File.exists? lane # Already been merged
      puts "#{merged_name} missing: #{lane}"
      missing_files = true
    end
  end
  next if missing_files
  
  `touch #{running_file}`
  puts "Merging #{lanes.length} lanes into #{merged_name}"
  begin
    for lane in lanes
      `samtools merge #{tmp_file} #{lanes.join(" ")}`
      FileUtils.mv(tmp_file, merged_filepath)
      `samtools index #{merged_filepath}`
    end
  rescue => e
    throw e
  ensure
    FileUtils.rm(tmp_file,        :force=>true)
    FileUtils.rm(merged_filepath, :force=>true)
    FileUtils.rm(running_file,    :force=>true)
  end
  
  
  
end