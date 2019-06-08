"""
    F = ExtensionField{F <: AbstractGaloisField, N, α, MinPoly}

Algebraic extension of a finite field ``F`` of degree ``N``.
"""
struct ExtensionField{F <: AbstractGaloisField, N, α, MinPoly, Conway} <: AbstractExtensionField
    n::NTuple{N, F}
    ExtensionField{F, N, α, MinPoly, Conway}(n::NTuple{N, F}) where
        {F <: AbstractGaloisField, N, α, MinPoly, Conway} = new(n)
end

basefield(::Type{ExtensionField{F, N, α, MinPoly, Conway}}) where {F, N, α, MinPoly, Conway} = F
char(::Type{ExtensionField{F, N, α, MinPoly, Conway}})      where {F, N, α, MinPoly, Conway} = char(F)
n(::Type{ExtensionField{F, N, α, MinPoly, Conway}})         where {F, N, α, MinPoly, Conway} = N
genname(::Type{ExtensionField{F, N, α, MinPoly, Conway}})   where {F, N, α, MinPoly, Conway} = α
minpoly(::Type{ExtensionField{F, N, α, MinPoly, Conway}})   where {F, N, α, MinPoly, Conway} = MinPoly
isconway(::Type{ExtensionField{F, N, α, MinPoly, Conway}})  where {F, N, α, MinPoly, Conway} = Conway

expansion(a::ExtensionField) = a.n

# -----------------------------------------------------------------------------
#
# Addition and substraction
#
# -----------------------------------------------------------------------------
zero(T::Type{<:ExtensionField}) =
    T(ntuple(i -> zero(basefield(T)), n(T)))
one( T::Type{<:ExtensionField}) =
    T(ntuple(i -> i == 1 ? one(basefield(T)) : zero(basefield(T)), n(T)))
gen( T::Type{<:ExtensionField}) =
    T(ntuple(i -> i == 2 ? one(basefield(T)) : zero(basefield(T)), n(T)))

+(a::F, b::F) where F<:ExtensionField = F(a.n .+ b.n)
-(a::F, b::F) where F<:ExtensionField = F(a.n .- b.n)

+(a::ExtensionField) = copy(a)
-(a::ExtensionField) = typeof(a)(.-a.n)

iszero(a::ExtensionField) = all(iszero, a.n)

# -----------------------------------------------------------------------------
#
# The interesting extension field operations: multiplication and division
#
# -----------------------------------------------------------------------------
function _rem(a::AbstractVector{C}, b::AbstractVector{C}) where C
    a = copy(a)
    len_a = findlast(!iszero, a)
    len_b = findlast(!iszero, b)
    while len_a !== nothing && len_a >= len_b
        q = a[len_a] // b[len_b]
        a[len_a-len_b+1:len_a] .-= q .* b[1:len_b]
        len_a = findprev(!iszero, a, len_a - 1)
    end
    a
end

function _gcdx(a::AbstractVector{C}, b::AbstractVector{C}) where C
    a = copy(a)
    b = copy(b)
    len_a = findlast(!iszero, a)
    len_b = findlast(!iszero, b)
    len_a === nothing && (len_a = 0)
    len_b === nothing && (len_b = 0)
    m = max(len_a, len_b)
    s0, s1 = zeros(C, m), zeros(C, m)
    t0, t1 = zeros(C, m), zeros(C, m)
    s0[1] = one(C)
    t1[1] = one(C)
    # The loop invariant is: s0*a0 + t0*b0 == a
    while len_b != 0
        s2 = copy(s0)
        t2 = copy(t0)
        while len_a >= len_b
            q = a[len_a] // b[len_b]
            a[len_a-len_b+1:len_a] .-= q .* b[1:len_b]

            deg_diff = len_a - len_b
            s2[deg_diff+1:end] .-= q .* s1[1:end-deg_diff]
            t2[deg_diff+1:end] .-= q .* t1[1:end-deg_diff]

            len_a = findprev(!iszero, a, len_a-1)
            len_a == nothing && (len_a = 0)
        end
        a, b = b, a
        len_a, len_b = len_b, len_a
        s0, s1 = s1, s2
        t0, t1 = t1, t2
    end
    return (a, s0, t0)
end

function *(::Direct, a::F, b::F) where F <: ExtensionField
    N = n(F)
    coeffs = zeros(basefield(F), 2N - 1)
    for (i, a_i) in enumerate(a.n)
        for (j, b_j) in enumerate(b.n)
            coeffs[i+j-1] += a_i * b_j
        end
    end
    coeffs = _rem(coeffs, collect(minpoly(F)))
    return F(ntuple(i -> coeffs[i], n(F)))
end

function inv(::Direct, a::F) where F <: ExtensionField
    iszero(a) && throw(DivideError())
    N = n(F)
    coeffs = collect(a.n)
    d, u, v = _gcdx(coeffs, collect(minpoly(F)))
    @assert !iszero(d[1]) && all(iszero, @view d[2:end])
    u ./= d[1]
    return F(ntuple(i -> u[i], n(F)))
end

/(::Direct, a::F, b::F) where F <: ExtensionField = a * inv(b)
//(::Direct, a::F, b::F) where F <: ExtensionField = a * inv(b)
^(::Direct, a::F, n::Integer) where F <: ExtensionField = Base.power_by_squaring(a, n)

*(a::F, b::F)       where F <: ExtensionField = zech_op(F, *, a, b)
/(a::F, b::F)       where F <: ExtensionField = zech_op(F, /, a, b)
//(a::F, b::F)      where F <: ExtensionField = zech_op(F, //, a, b)
^(a::F, n::Integer) where F <: ExtensionField = zech_op(F, ^, a, n)
inv(a::F)           where F <: ExtensionField = zech_op(F, inv, a)

# -----------------------------------------------------------------------------
#
# Constructors and promotion
#
# -----------------------------------------------------------------------------
promote_rule(F::Type{<:ExtensionField}, ::Type{<:Integer}) = F
function convert(F::Type{<:ExtensionField}, i::Integer)
    convert(F, convert(basefield(F), i))
end

(::Type{F})(n::F) where F<:ExtensionField = F(n.n)

# -----------------------------------------------------------------------------
#
# Random number
#
# -----------------------------------------------------------------------------
function rand(rng::AbstractRNG, ::SamplerType{F}) where F <: ExtensionField
    F(ntuple(i -> rand(rng, basefield(F)), n(F)))
end
