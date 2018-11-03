using Serialization
using HTTP
using JSON


url = "https://gist.githubusercontent.com/tkluck/e1cd1746c69aa17e4a37114d22649627/raw/7fbe9763fae27f14924262ad03606f1c3af4400e/CPImport.txt"
ResDataType = Dict{Tuple{Int,Int}, Vector{Int}}
result = ResDataType()

data = String(HTTP.request("GET", url).body)
data = replace(data, r"allConwayPolynomials := " => "")
data = replace(data, r"];"                       => "]")

list = JSON.parse(data)
# Frank Lübeck's data format signifies the end of the array with a null marker
pop!(list)

result = ResDataType((p,n) => coeffs for (p,n,coeffs) in list)

open("conwaypolynomials.data", "w") do f
    serialize(f, result)
end
