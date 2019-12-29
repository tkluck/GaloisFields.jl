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
    if n == 1
        D = GaloisField(p, n, :α)
        α = rand(D)
        E, _ = FiniteField(p, n, "β")
        β = rand(E)
    else
        D, α = GaloisField(p, n, :α)
        E, β = FiniteField(p, n, "β")
    end
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

# add a comparison between a 'naive' broadcast and the broadcast machinery
# (This is not a comparison against Nemo but still nice to have in the same place.)
let p = prevprime(typemax(Int16))
    G = GaloisField(p, 1)
    x = rand(0:p-1, 1_000_000)
    y = rand(0:p-1, 1_000_000)
    z = rand(0:p-1)
    galoisfields["broadcast"] = BenchmarkGroup()
    nemo["broadcast"] = BenchmarkGroup()

    galoisfields["broadcast"]["add"] = @benchmarkable $(G.(x)) .+ $(G.(y))
    nemo["broadcast"]["add"] = @benchmarkable mod.($x .+ $y, $p)

    galoisfields["broadcast"]["muladd"] = @benchmarkable $(G.(x)) .+ $(G(z)) .* $(G.(y))
    nemo["broadcast"]["muladd"] = @benchmarkable mod.($x .+ $z .* $y, $p)
end

println(judge(ratio(median(run(galoisfields)), median(run(nemo)))))
