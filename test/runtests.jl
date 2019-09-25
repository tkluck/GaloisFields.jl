tasks = map(["using-bitintegers.jl", "not-using-bitintegers.jl"]) do script
    run(`$(julia_cmd()) $(joinpath(@__DIR__, script)))`)
end

for t in tasks
    wait(t)
end

success = maximum(t -> t.exitcode, t in tasks)
exit(success)
