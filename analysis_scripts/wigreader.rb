# Ron Shlomo Gejman
# ron.gejman@gmail.com
# All rights granted for commercial and non-commercial use

class WigReader
  
  
  # WigReader is a class for loading and using variableStep Wiggle files
  # WARNING: The entire contents of a wiggle file are loaded into memory TWICE.
    
  # Store all the wiggle entries in @data.
  attr_accessor :data

  def initialize(filename)
    @data = {}
    chr   = nil
    step  = nil
    
    ## MEMORY INTENSIVE STEP
    
    lines = File.readlines(filename)
    
    ## END MEMORY INTENSIVE STEP
    
    skip_until_header = false 
    lines_skipped = 0
    last_pos      = nil
    
    # Digest Wiggle file, line by line.
    # When a header line is encountered we initialize a new hash (inside @data) for the chromosome.
    
    while line = lines.shift
      line.chomp!
      if line[0,1] == "v" #e.g. variableStep chrom=chr11 span=25
        skip_until_header = false
        tmp, chr, step = line.split(" ").collect{|a| a.split("=")[1]}
        step = step.to_i
        if @data.has_key? chr
          # Ends of chromosomes may not line up with the "step" size. If this happens, the Wiggle file will have a new
          # header with a new stepsize. If we encounter the same header twice, we skip all the entries for the second header
          # since it is probably just 1 line at the end of the chromosome
          # This is not optimal or always safe, but allows for easy summarization w/o missing much data.
          # The "lines skipped" variable will print at the end to tell you how many lines you missed this way.
          # It should not be more than the # of chromosomes!
          skip_until_header = true
          next
        end
        puts "Reading #{chr} with #{step}nt steps"
        @data[chr] = {:step=>step, :positions=>{}}
        next
      end
      if skip_until_header
        lines_skipped += 1
      end
      pos,fpkm = line.split(" ")
      pos = pos.to_i
      
      # Wiggle files are not required to be contiguous. So the position "1025" may follow the position "900."
      # In this case, we fill in the missing values with "0.0"
      while !last_pos.nil? and last_pos < (pos-step)
        last_pos += step
        @data[chr][:positions][last_pos] = 0.0
      end
      @data[chr][:positions][pos] = fpkm.to_f
      last_pos = pos
    end
    
    # Set the "lines" array to nil and hope that the GC notices it.
    lines = nil
    GC.start
    
    # After we have read in the entire dataset, we go back and find the first and last positions on the chromosome.
    # This could be a tad faster if performed while reading... but it clutters up the code too much to do it that way
    
    for chr in @data.keys
      positions = @data[chr][:positions].keys.sort
      @data[chr][:start]  = positions.first
      @data[chr][:end]    = positions.last
    end
    
    # Print how many lines were skipped (hopefully at the end of the chromosomes)
    # Should not be more lines than the # of chromosomes.
    puts "Skipped #{lines_skipped} lines."
    puts "Chrs: #{@data.keys.join(',')}"
    #pp @data
  end
  
  #Return the FPKM value for the region between coordinates "a" and "b"
  def fpkm(chr,a,b)
    if !@data.has_key? chr
      if chr == "MT"
        new_key = "chrM"
      else
        new_key = "chr" + chr
      end
      if !@data.has_key? new_key
        raise "Neither #{chr} nor #{new_key} are in the wig file"
      else
        chr = new_key
      end
    end
    raise "Start and End cannot be equal" if a == b
    start = @data[chr][:start]
    last  = @data[chr][:end]
    step  = @data[chr][:step]
    if a < b
      raise "Start (#{a}) is before first pos on this chromosome (#{start})" if a < start
      raise "End (#{b}) is after last pos on this chromosome (#{last})" if b > last
      return calculate_fpkm(chr,a,b,start,step)
    else #reverse strand
      raise "Start (#{b}) is before first pos on this chromosome (#{start})" if b < start
      raise "End (#{b}) is after last pos on this chromosome (#{last})" if a > last
      return calculate_fpkm(chr,b,a,start,step)
    end
  end
  
  private
  def calculate_fpkm(chr,a,b,start,step)
    numerator   = 0.0
    denominator = 0.0
    
    # The "pos" variable holds the current position along the chromosome.
    # We calculate the initial position by finding the first wiggle block containing "a"
    pos = ((((a.to_f - start.to_f)/step.to_f).floor * step) + start).to_i
    
    # If the first base of the wiggle block is not equal to the first base we are interested in
    # we must take a fraction of that wiggle block (equal to the # of bases overlapping with the region of interest)
    if pos < a
      bases = step - (start - pos)
      raise "pos<a: # of bases (#{bases}) is < 0!: #{step} #{start} #{pos}" if bases < 0
      numerator += @data[chr][:positions][pos] * (bases) / step.to_f # This is the fraction of the first block
      denominator += (bases / step.to_f)
      pos += step
    end
    # for the rest of the wiggle blocks we just take the entire block value
    until pos > (b-step)
      numerator += @data[chr][:positions][pos]
      denominator += 1
      pos += step
    end
    # if the region spills over into another partial wiggle block then we also have to take a fraction of that block
    if pos <= b #we have to add the partial last fragment
      bases = b - pos + 1
      numerator += @data[chr][:positions][pos] * bases / step.to_f
      denominator += (bases / step.to_f)
    end
    return numerator.to_f / denominator
  end
end

if __FILE__ == $0
  #reader = WigReader.new("smalltest.wig")
  #reader = nil
  reader = WigReader.new("bigtest.wig")
end