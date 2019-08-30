using Revise, BenchmarkTools, GaloisFields, Nemo
using Primes
using BenchmarkTools: BenchmarkGroup, run, median, ratio, judge

galoisfields = BenchmarkGroup()
nemo         = BenchmarkGroup()

function benchmarks!(grp, x)
    grp["square"]                 = @benchmarkable $x^2
    grp["sum"]                    = @benchmarkable $x + $x
    grp["polynomial"]             = @benchmarkable $x^10 + 5*$x^8 + 10*$x^6 + 28*$x^4 + $x + 3
    grp["integer multiplication"] = @benchmarkable 5 * $x
    grp["large power"]            = @benchmarkable $x^1000
    grp
end

function compare!(p, n)
    D, α = GaloisField(p, n, :α)
    E, β = FiniteField(p, n, "β")
    galoisfields["q=$p^$n"] = benchmarks!(BenchmarkGroup(), α)
    nemo["q=$p^$n"] = benchmarks!(BenchmarkGroup(), β)
end

compare!(2, 1)
compare!(2, 10)
compare!(2, 64)
compare!(3, 1)
compare!(3, 2)
compare!(3, 8)
compare!(29, 1)
compare!(29, 3)
compare!(31, 1)
compare!(prevprime(typemax(Int128)), 1)

# we don't have the Conway polynomial for 29^42 in our database,
# so we provide it manually.
const F = @GaloisField! 29 γ^42 + -17γ^37 - 2
const G, δ = FiniteField(29, 42, "δ")
galoisfields["q=29^42"] = benchmarks!(BenchmarkGroup(), γ)
nemo["q=29^42"] = benchmarks!(BenchmarkGroup(), δ)

println(judge(ratio(median(run(galoisfields)), median(run(nemo)))))
