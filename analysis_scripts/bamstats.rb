#!/usr/bin/ruby -wKU

if ARGV.length != 1
  puts "Usage: bamstats.rb *.bam"
end
s = ARGV[0]
files = Dir[s]
files = files.select {|f| f =~ /\.bam$/}

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
      out = []
      io_threads = []
      io_threads << Thread.new(stderr) do |terr|
        while (line = terr.gets)
          puts "stderr: #{line}"
        end
      end
      io_threads << Thread.new(stdout) do |tout|
        while (line = tout.gets)
          next unless line =~ /(Total)|(Mapped)/
          out << line
        end
      end
      io_threads.each{|t| t.join()} #in order to cleanup when you're done.
      semaphore.synchronize do
        puts file
        out.each {|o| puts o}
      end
    }
  end
end