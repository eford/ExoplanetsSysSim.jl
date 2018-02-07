if length(ARGS)<2
   warn("# usage:  convert_csv_to_jld.jl input.csv")
   warn("# WARNING:  THIS DOESN'T OUTPUT A JLD FILE YET!")
   exit(255)
end

csv_file = convert(String,ARGS[1])
if !isfile(csv_file)
   warn(string("Can't read ",csv_file))
   exit(255)
end
jld_file = replace(csv_file,r".csv$",".jld")


filename = csv_file
(data,header) = readcsv(filename,header=true)

ncols = length(header)
data_dict = Dict{String,Any}()
for i in 1:ncols
  data_vec = data[:,i]
  try 
    data_vec = convert(Vector{Float64},data[:,i])
    #dat_dict[header[i]] = data_vec
    println("# Converted col:", header[i])
  catch
    #println("# Couldn't convert colm: ",header[i])
  end
end


