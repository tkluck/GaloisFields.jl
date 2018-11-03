using Serialization
using JSON

ResDataType = Dict{Tuple{Int,Int}, Vector{Int}}
result = ResDataType()

try
    data = read(`curl -L 'http://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/CPimport.txt'`)
    data = transcode(String, data)

    data = replace(data, r"allConwayPolynomials := " => "")
    data = replace(data, r"];"                       => "]")

    list = JSON.parse(data)
    # Frank Lübeck's data format signifies the end of the array with a null marker
    pop!(list)

    result = ResDataType((p,n) => coeffs for (p,n,coeffs) in list)
catch exception
    @warn """
        Issue downloading Conway polynomials database: $exception
        Please retry the build step again later.
    """
end

open("conwaypolynomials.data", "w") do f
    serialize(f, result)
end
