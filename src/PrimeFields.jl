import Base.Checked: add_with_overflow, sub_with_overflow

import .Util: hilo_mul, widen_bits

"""
    posmod(x, y)

Return mod(x, y) under the assumption that y is positive. In general,
returns

    rem(x, y) + ifelse(rem(x, y) < 0, y, 0)

This is slightly faster than the mod operation.
"""
@inline function posmod(x, y)
    z = rem(x, y)
    z + ifelse(signbit(z), y, zero(y))
end

"""
    PrimeField{I<:Integer, p}

A type representing an element in ℤ/pℤ.
"""
struct PrimeField{I<:Integer, p} <: AbstractGaloisField
    n::I
    PrimeField{I, p}(n::Integer) where {I, p} = new(posmod(n, p))
    PrimeField{I, p}(::NonNegative, n::Integer) where {I, p} = new(rem(n, p))
    PrimeField{I, p}(::Reduced, n::Integer) where {I, p} = new(n)
end

char(::Type{PrimeField{I,p}})    where {I, p} = p
n(::Type{PrimeField{I,p}})       where {I, p} = 1
inttype(::Type{PrimeField{I,p}}) where {I, p} = I

char(a::AbstractGaloisField)    = char(typeof(a))
n(a::AbstractGaloisField)       = n(typeof(a))
inttype(a::AbstractGaloisField) = inttype(typeof(a))

function inttype(p::Integer)
    T = Int8
    while isbitstype(T)
        p <= typemax(T) && return T
        T = widen_bits(T)
    end
    error("Integer type of prime ($p) must be a bitstype.")
end

# -----------------------------------------------------------------------------
#
# Arithmetic
#
# -----------------------------------------------------------------------------
zero(F::Type{<:PrimeField}) = F(Reduced(), zero(inttype(F)))
one( F::Type{<:PrimeField}) = F(Reduced(),  one(inttype(F)))

@Base.pure function _power_of_two(n, p)
    return rem(big"2"^n, p) % typeof(p)
end

@Base.pure _overflow_value(F, I) = _power_of_two(leading_zeros(zero(I)), char(F))
@Base.pure _overflow_value(F) = _overflow_value(F, inttype(F))

function +(a::F, b::F) where F<:PrimeField
    iszero(a) && return +b
    iszero(b) && return +a
    if char(F) > div(typemax(inttype(F)), 2)
        m, overflowed = add_with_overflow(a.n, b.n)
        if overflowed
            # the effect of the overflow is that typemax(I) + 1
            # is replaced by typemin(I). In other words, we
            # basically substracted
            #     typemax(I) + 1 - typemin(I)
            # We need to add it back to make up for this. And
            # because we're about to do a mod p reduction,
            # we can add back the mod p residue class.
            m += 2rem(typemax(inttype(F)), char(F)) + 2
            # TODO: can this last operation overflow?
        end
        F(m) # note: not necessarily non-negative in this case!
    else
        F(NonNegative(), a.n + b.n)
    end
end

function -(a::F, b::F) where F<:PrimeField
    iszero(a) && return -b
    iszero(b) && return +a
    if char(F) > div(typemax(inttype(F)), 2)
        m, overflowed = sub_with_overflow(a.n, b.n)
        if overflowed
            # see above for an explanation of this magic
            # number.
            m -= 2rem(typemax(inttype(F)), char(F)) + 2
        end
        F(m) # note: not necessarily non-negative in this case!
    else
        F(NonNegative(), char(F) + a.n - b.n)
    end
end

+(a::PrimeField) = a
-(a::PrimeField) = iszero(a.n) ? a : typeof(a)(Reduced(), char(typeof(a)) - a.n)

function *(a::F, b::F) where F<:PrimeField
    iszero(a) && return zero(F)
    iszero(b) && return zero(F)
    isone(a) && return +b
    isone(b) && return +a
    F(NonNegative(), Base.widemul(a.n, b.n))
end

*(a::F, b::F) where F <: PrimeField{Int128} = a * b.n # use allocation-free function below

# a version that avoids widemul (which allocates a BigInt for Int128)
function *(a::F, b::Union{Int8, Int16, Int32, Int64, Int128}) where F <: PrimeField{Int128}
    iszero(a) && return zero(F)
    iszero(b) && return zero(F)
    isone(a) && return F(b)
    isone(b) && return +a

    u = a.n
    v = Int128(b)
    res = zero(UInt128)

    # TODO: prove that this terminates
    while true
        hi, lo = hilo_mul(u, v)

        res, overflow = add_with_overflow(res, lo)
        hi += overflow

        iszero(hi) && return F(res)

        u = hi
        v = _overflow_value(F, UInt128)
    end
end

function *(a::PrimeField, b::Integer)
    iszero(a) && return zero(a)
    iszero(b) && return zero(a)
    isone(a) && return oftype(a, b)
    isone(b) && return +a
    oftype(a, Base.widemul(a.n, b))
end

*(a::Integer, b::PrimeField) = b * a

"""
    _invmod(n, m)

Multiplicative inverse of n (mod m), assuming 0 <= n < m and m prime.
"""
function _invmod(n::T, m::T) where T<:Integer
    g, x, y = gcdx(n, m)
    r = x + ifelse(signbit(x), m, zero(m))
    r
end
_invmod(n::Integer, m::Integer) = _invmod(promote(n, m)...)

inv(a::F)     where F<:PrimeField = iszero(a) ? throw(DivideError()) : F(Reduced(), _invmod(a.n, char(F)))
/(a::F,b::F)  where F<:PrimeField = a * inv(b)
//(a::F,b::F) where F<:PrimeField = a * inv(b)

iszero(a::PrimeField) = iszero(a.n)

# -----------------------------------------------------------------------------
#
# Constructors and promotion
#
# -----------------------------------------------------------------------------
promote_rule(F::Type{<:PrimeField}, ::Type{<:Integer}) = F
convert(F::Type{<:PrimeField}, i::Integer) = F(i)

(::Type{F})(n::F) where F<:PrimeField = F(Reduced(), n.n)
convert(::Type{F}, n::F) where F<:PrimeField = n

# -----------------------------------------------------------------------------
#
# Random number
#
# -----------------------------------------------------------------------------
function rand(rng::AbstractRNG, ::SamplerType{F}) where F <: PrimeField
    F(rand(rng, 0 : char(F)-1))
end
