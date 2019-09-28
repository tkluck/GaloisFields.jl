using Primes
using GaloisFields

using BitIntegers

const G = @GaloisField! 𝔽₂₉ α^2 - 2
const H = @GaloisField! G   β^3 + 2β + 1
const J = @GaloisField! H   γ^7 - 2

const TestFields = [
    @GaloisField ℤ/2ℤ
    @GaloisField ℤ/3ℤ
    @GaloisField ℤ/5ℤ
    @GaloisField ℤ/7ℤ
    @GaloisField ℤ/67ℤ

    @GaloisField! 𝔽₄ α
    @GaloisField! 𝔽₉ β

    @GaloisField! 𝔽₃ α^2 + 1

    G
    H
    J
    [GaloisField(prevprime(typemax(I)))
     for I in [Int8, Int16, Int32, Int64, Int128]]

    @GaloisField! 𝔽₂ α^2 + α + 1
    @GaloisField! 𝔽₅ α^2 - 2
    @GaloisField! 2^2 α
    @GaloisField! 5^2 α
    @GaloisField! 5^6 α

    # BitIntegers
    GaloisField(nextprime(Int256(typemax(Int128))+1))
    GaloisField(prevprime(typemax(Int256)))
]

include("arithmetic.jl")
