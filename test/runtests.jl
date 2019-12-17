using GaloisFields # make sure pre-compilation happens in a single process

tasks = mapreduce(&, ["using-bitintegers.jl", "not-using-bitintegers.jl"]) do script
    fullpath = joinpath(@__DIR__, script)
    `$(Base.julia_cmd()) $fullpath`
end

run(tasks)
