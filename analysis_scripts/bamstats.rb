#!/usr/bin/ruby -wKU
require 'thread'
require "open3"

if ARGV.length < 1
  puts "Usage: bamstats.rb *.bam"
  exit(1)
end
files = ARGV.select {|f| f =~ /\.bam$/}
puts "Processing:"
files.each {|f| puts f}

if files.length == 0
  puts "No BAM files found."
end

threads = []

semaphore = Mutex.new


for file in files
  threads << Thread.new do
    Open3.popen3("bamtools stats -in #{file}") { |stdin, stdout, stderr|
      # stdin = input stream
      # stdout = output stream
      # stderr = stderr stream
      out = nil
      io_threads = []
      io_threads << Thread.new(stderr) do |terr|
        while (line = terr.gets)
          puts "stderr: #{line}"
        end
      end
      io_threads << Thread.new(stdout) do |tout|
        while (line = tout.gets)
          next unless line =~ /(Mapped)/
          out = line
        end
      end
      io_threads.each{|t| t.join()} #in order to cleanup when you're done.
      semaphore.synchronize do
        puts "#{file}: #{out}"
      end
    }
  end
end

threads.each {|t| t.join()}