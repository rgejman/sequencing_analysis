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
    headers = `grep -n v #{filename}`.split("\n").collect{|l| l.split(":")[0].to_i} #subtract one to get base-0 pos.
    puts "done"
    
    pipes = []
    forks = []
    for header_pos_index in (0...headers.length)
      header_pos = headers[header_pos_index]
      next_header_pos = header_pos_index < (headers.length-1) ? headers[header_pos_index+1] : nil
      rd, wr = IO.pipe
      forks << fork do
        rd.close
        if next_header_pos.nil?
          lines = `tail -n +#{header_pos} #{filename}`.split("\n")
        else
          len = next_header_pos-header_pos
          lines = `tail -n +#{header_pos} #{filename} | head -n #{len}`.split("\n")
        end
        line = lines.shift.chomp
        raise "ERROR: This was supposed to be a header line. Instead got: #{line} for #{header_pos}." if line == nil or line[0,1] != "v"
        raise "ERROR: Last line should be data, not header or nil. #{next_header_pos}. #{lines.last} | #{lines[lines.length-2]}" if lines.last == nil or lines.last[0,1] == "v"
        tmp, chr, step = line.split(" ").collect{|a| a.split("=")[1]}
        step = step.to_i
        d = {}
        d[chr] = {}
        d[chr] = {:step=>step, :positions=>{}, :chr=>chr}
        last_pos      = nil
        
        # this is a partial chr. No need to grab it.
        if lines.length < 5
          wr.close
          exit(0)
        end
        puts "Reading #{chr} with #{step}nt steps and #{lines.length} lines}"
        for i in (0...lines.length)
          line = lines[i].chomp
          break if line[0,1] == "v"
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
        puts "Done reading #{chr} with #{step}nt steps and #{lines.length} lines"
        for position in d[chr][:positions].keys.sort
          wr.puts "#{chr} #{step} #{position} #{d[chr][:positions][position]}"
        end
        wr.close
        puts "Done passing #{chr} to pipe. Will exit"
        exit(0)
      end
      wr.close
      pipes << rd
    end
    
    #If there are pipes whose contents have not been emptied:
    while pipes.any? {|rd| !rd.eof? } do
      for rd in pipes
        next if rd.eof?
        tokens = rd.readline.chomp.split(" ")
        chr = tokens[0]
        step = tokens[1].to_i
        pos = tokens[2].to_i
        val = tokens[3].to_f
        @data[chr] ||= {:step=>step}
        @data[chr][:positions] ||= {}
        @data[chr][:positions][pos] = val
      end
    end
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

if __FILE__ == $0
  f = "/media/bigdisk/sequencing/wig/Eugene/Eugene_WT_H3K9me2_CD4_all.sorted.wig"
  reader = WigReaderFast.new(f)
end