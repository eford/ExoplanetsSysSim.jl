filename = "q1_q17_dr25_koi.csv"
#filename = "q1_q17_dr25_stellar.csv"
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



