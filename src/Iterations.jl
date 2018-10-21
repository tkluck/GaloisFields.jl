
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
