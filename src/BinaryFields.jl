# we use the Bits singleton to avoid ambiguity between the constructor that
# takes a bit pattern and the constructor that takes an Integer value (to be
# interpreted mod 2).
struct Bits end
struct BinaryField{I <: Unsigned, N, α, MinPolyMask, Conway} <: AbstractExtensionField
    n::I
    BinaryField{I, N, α, MinPolyMask, Conway}(::Bits, n) where
        {I <: Unsigned, N, α, MinPolyMask, Conway} = new(n)
end

basefield(::Type{BinaryField{I, N, α, MinPolyMask, Conway}})   where {I, N, α, MinPolyMask, Conway} = GaloisField(2)
char(::Type{BinaryField{I, N, α, MinPolyMask, Conway}})        where {I, N, α, MinPolyMask, Conway} = 2
n(::Type{BinaryField{I, N, α, MinPolyMask, Conway}})           where {I, N, α, MinPolyMask, Conway} = N
genname(::Type{BinaryField{I, N, α, MinPolyMask, Conway}})     where {I, N, α, MinPolyMask, Conway} = α
minpolymask(::Type{BinaryField{I, N, α, MinPolyMask, Conway}}) where {I, N, α, MinPolyMask, Conway} = MinPolyMask
isconway(::Type{BinaryField{I, N, α, MinPolyMask, Conway}})    where {I, N, α, MinPolyMask, Conway} = Conway

struct BitsIter{I <: Integer}
    n::I
end
expansion(a::BinaryField) = BitsIter(a.n)

Base.iterate(it::BitsIter) = iszero(it.n) ? nothing : (GaloisField(2)(it.n & 1), 1)
function Base.iterate(it::BitsIter, state)
    if it.n < big"1" << state
        return nothing
    else
        GaloisField(2)((it.n & (1 << state)) >> state), state + 1
    end
end
Base.length(it::BitsIter) = iszero(it.n) ? 0 : trailing_zeros(prevpow(2, it.n)) + 1
Base.eltype(it::BitsIter) = GaloisField(2)

function minpoly(F::Type{<:BinaryField})
    minpoly = big"1" << n(F) | minpolymask(F)
    return BitsIter(minpoly)
end

# -----------------------------------------------------------------------------
#
# Addition and substraction
#
# -----------------------------------------------------------------------------
zero(T::Type{<:BinaryField{I}}) where I <: Unsigned = T(Bits(), I(0))
one( T::Type{<:BinaryField{I}}) where I <: Unsigned = T(Bits(), I(1))
gen( T::Type{<:BinaryField{I}}) where I <: Unsigned = T(Bits(), I(2))

+(a::F, b::F) where F<:BinaryField = F(Bits(), a.n ⊻ b.n)
-(a::F, b::F) where F<:BinaryField = F(Bits(), a.n ⊻ b.n)

+(a::BinaryField) = copy(a)
-(a::BinaryField) = copy(a)

iszero(a::BinaryField) = iszero(a.n)

# -----------------------------------------------------------------------------
#
# The interesting extension field operations: multiplication and division
#
# -----------------------------------------------------------------------------
_widen(::Type{UInt8}) = UInt16
_widen(::Type{UInt16}) = UInt32
_widen(::Type{UInt32}) = UInt64
_widen(::Type{UInt64}) = UInt128
include("BitReversal.jl")
function carrylessmul(a::I, b::I) where I <: Integer
    J = _widen(I)
    res = zero(J)
    a = J(a)
    b = bitswap(J(b))
    for _ = 1:8sizeof(J)
        res <<= 1
        res |= (count_ones(b & a) & 1) % J
        b >>= 1
    end
    return res
end

function *(::Direct, a::F, b::F) where F <: BinaryField{I} where I
    N = n(F)
    c = carrylessmul(a.n, b.n)
    J = typeof(c)

    minpoly = (J(1) << N) | minpolymask(F)
    while (lz = leading_zeros(c)) < (8sizeof(J) - N)
        to_shift = 8sizeof(J) - lz - N - 1
        c ⊻= minpoly << to_shift
    end
    return F(Bits(), c % I)
end

function _gcdx(::Bits, a::I, b::I) where I <: Unsigned
    a = copy(a)
    b = copy(b)
    len_a = 8sizeof(I) - leading_zeros(a)
    len_b = 8sizeof(I) - leading_zeros(b)
    m = max(len_a, len_b)
    s0, s1 = one(I), zero(I)
    t0, t1 = zero(I), one(I)
    # The loop invariant is: s0*a0 + t0*b0 == a
    while len_b != 0
        s2 = copy(s0)
        t2 = copy(t0)
        while len_a >= len_b
            deg_diff = len_a - len_b
            a ⊻= b << deg_diff
            s2 ⊻= s1 << deg_diff
            t2 ⊻= t1 << deg_diff

            len_a = 8sizeof(I) - leading_zeros(a)
        end
        a, b = b, a
        len_a, len_b = len_b, len_a
        s0, s1 = s1, s2
        t0, t1 = t1, t2
    end
    return (a, s0, t0)
end

function inv(::Direct, a::F) where F <: BinaryField{I} where I
    iszero(a) && throw(DivideError())
    N = n(F)
    # the minimum polynomial is (1 << N) | minpolymask(F)
    # and if N is maximal, this takes us outside of
    # the allocated bit size. So we do the first step
    # of the Euclidean algorithm "manually" instead of
    # letting _gcdx do it.
    len_a = (8sizeof(I) - leading_zeros(a.n))
    len_diff = (N + 1) - len_a
    minpoly = (I(1) << N) | minpolymask(F)
    b = minpoly ⊻ (a.n << len_diff)
    d, u, v = _gcdx(Bits(), a.n, b)
    # apply the result of the "manual" first step
    u ⊻= carrylessmul(I(1) << len_diff, v)
    @assert d == I(1)
    return F(Bits(), u % I)
end

/(::Direct, a::F, b::F)  where F <: BinaryField = a * inv(b)
//(::Direct, a::F, b::F) where F <: BinaryField = a * inv(b)
^(::Direct, a::F, n::Integer) where F <: BinaryField = Base.power_by_squaring(a, n)

*(a::F, b::F)       where F <: BinaryField = zech_op(F, *, a, b)
/(a::F, b::F)       where F <: BinaryField = zech_op(F, /, a, b)
//(a::F, b::F)      where F <: BinaryField = zech_op(F, //, a, b)
^(a::F, n::Integer) where F <: BinaryField = zech_op(F, ^, a, n)
inv(a::F)           where F <: BinaryField = zech_op(F, inv, a)

# -----------------------------------------------------------------------------
#
# Constructors and promotion
#
# -----------------------------------------------------------------------------
promote_rule(F::Type{<:BinaryField}, ::Type{<:Integer}) = F
function convert(F::Type{<:BinaryField{I}}, i::Integer) where I
    F(Bits(), (i & 1) % I)
end

(::Type{F})(n::F) where F<:BinaryField = F(Bits(), n.n)
