function readInstance(file::String)

    #Open file, extract data, close file
    o = open(file,"r")
    data = readlines(o)
    close(o)
  
    #Get simple numeric data
    n = parse(Int64,data[1])
    d = parse(Int64,data[2])
    f = parse(Int64,data[3])
    Amin = parse(Int64,data[4])
    Nr = parse(Int64,data[5])
    R = parse(Int64,data[9])
  
    #Get region definitions
    regions = Dict()
    for i in 1:Nr
      regions[i] = []
    end
  
    region_data = split(data[7]," ")
    for i in 1:n
      k = parse(Int64,region_data[i])
      if (k != 0)
        append!(regions[k],i)
      end
    end
  
    #Get airport coordinates
    coords = Matrix{Int64}(zeros(n,2))
    for i in 1:n
      line = split(data[10+i]," ")
      coords[i,1] = parse(Int64,line[1])
      coords[i,2] = parse(Int64,line[2])
    end
  
    return n,d,f,Amin,Nr,R,regions,coords
  end