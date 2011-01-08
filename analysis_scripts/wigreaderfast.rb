# Ron Shlomo Gejman
# ron.gejman@gmail.com
# All rights granted for commercial and non-commercial use
require 'wigreader.rb'

class WigReaderFast < WigReader

  def initialize(filename)
    @data = {}
    chr   = nil
    step  = nil
    print "Getting header locations..."
    # Use grep to identify the locations of the chromosome headers in the wiggle file
    headers = `grep -n v #{filename}`.split("\n").collect{|l| l.split(":")[0].to_i - 1} #subtract one to get base-0 pos.
    puts "done"
    ## MEMORY INTENSIVE STEP
    print "Slurping wiggle file..."
    lines = File.readlines(filename)
    puts "done"

    puts "Parsing wiggle file"
    ## END MEMORY INTENSIVE STEP
    pipes = []
    forks = []
    for header_pos_index in (0...headers.length)
      header_pos = headers[header_pos_index]
      next_header_pos = header_pos_index < (headers.length-1) ? headers[header_pos_index]+1 : nil
      rd, wr = IO.pipe
      forks << fork do
        rd.close
        pre = lines.length
        if next_header_pos.nil?
          lines = lines[header_pos..-1]
        else
          lines = lines[header_pos...next_header_pos] #Trim the array; keep only the lines for my header
        end
        line = lines.shift.chomp
        raise "ERROR: This was supposed to be a header line. Instead got: #{line} for #{header_pos}." if line == nil or line[0,1] != "v"
        tmp, chr, step = line.split(" ").collect{|a| a.split("=")[1]}
        step = step.to_i
        d = {}
        d[chr] = {}
        d[chr] = {:step=>step, :positions=>{}, :chr=>chr}
        last_pos      = nil
        puts "Reading #{chr} with #{step}nt steps"
        for i in (0...lines.length)
          line = lines[i].chomp!
          pos,fpkm = line.split(" ")
          pos = pos.to_i

          # Wiggle files are not required to be contiguous. So the position "1025" may follow the position "900."
          # In this case, we fill in the missing values with "0.0"
          while !last_pos.nil? and last_pos < (pos-step)
            last_pos += step
            d[chr][:positions][last_pos] = 0.0
          end
          d[chr][:positions][pos] = fpkm.to_f
          last_pos = pos
        end
        puts "Done reading #{chr} with #{step}nt steps"
        wr.puts chr 
        wr.puts step
        wr.flush
        for position in d[chr][:positions].keys.sort
          wr.puts "#{position} #{d[chr][:positions][position]}"
          wr.flush
        end
        wr.close
        puts "Done passing #{chr} to pipe. Will exit"
        exit(0)
      end
      pipes << [rd,wr]
    end
    puts "Waiting for processes to finish"
    forks.each {|id| Process.wait(id); puts "Process #{id} done"}
    puts "Finished waiting for processes"
    for rd,wr in pipes
      puts "Reading from pipe"
      wr.close
      d = rd.read.split("\n")
      rd.close
      chr = d.shift
      step = d.shift
      positions = {}
      for l in d #iterate over the array of positions
        pos,fpkm = l.split(" ")
        positions[pos.to_i] = fpkm.to_f
      end
      next if positions.keys.length < 5 # this is the continuation of a chr that did not fit evently into the wig blocks
      raise "Unexpected error" if @data.has_key? chr
      @data[chr] = {:step=>step,:positions=>positions}
    end

    lines_skipped = 0

    # After we have read in the entire dataset, we go back and find the first and last positions on the chromosome.
    # This could be a tad faster if performed while reading... but it clutters up the code too much to do it that way

    for chr in @data.keys
      positions = @data[chr][:positions].keys.sort
      @data[chr][:start]  = positions.first
      @data[chr][:end]    = positions.last
    end

    # Set the "lines" array to nil and hope that the GC notices it.
    lines = nil

    # Print how many lines were skipped (hopefully at the end of the chromosomes)
    # Should not be more lines than the # of chromosomes.
    puts "Skipped #{lines_skipped} lines."
    #pp @data
  end
end
