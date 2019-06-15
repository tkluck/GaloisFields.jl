
iterate(K::Type{<:PrimeField}) = (K(0), 1)
iterate(K::Type{<:PrimeField}, n) = n < char(K) ? (K(n), n + 1) : nothing

iter(K::Type{<:ExtensionField}) = Iterators.product([basefield(K) for _ = 1:n(K)]...)
function iterate(K::Type{<:ExtensionField})
    it = iterate(iter(K))
    it == nothing && return nothing
    return K(it[1]), it[2]
end
function iterate(K::Type{<:ExtensionField}, state)
    it = iterate(iter(K), state)
    it == nothing && return nothing
    return K(it[1]), it[2]
end

iterate(K::Type{<:BinaryField}) = K(Bits(), 0), 1
iterate(K::Type{<:BinaryField}, state) = state < 2^n(K) ? (K(Bits(), state), state + 1) : nothing

"""
    Yield q^n in the smallest Signed that fits it
    and is at least Int.
"""
function safepow(q::Signed, n)
    q = q isa Base.SmallSigned ? Int(q) : q
    bits = leading_zeros(zero(q)) -  leading_zeros(q)
    m = n * bits
    if m <= leading_zeros(zero(Int))
        return q^n
    elseif m <= leading_zeros(zero(Int64))
        return Int64(q)^n
    elseif m <= leading_zeros(zero(Int128))
        return Int128(q)^n
    else
        return BigInt(q)^n
    end
end

Base.IteratorSize(::Type{<:AbstractGaloisField}) = Base.HasLength()
Base.length(F::Type{<:PrimeField}) = char(F)
@generated Base.length(::Type{F}) where F <: AbstractGaloisField = safepow(length(basefield(F)), n(F))
