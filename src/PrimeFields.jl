import Base.Checked: add_with_overflow, sub_with_overflow

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
    for I in [Int8, Int16, Int32, Int64, Int128]
        if p <= typemax(I)
            return I
        end
    end
    throw("Primes greater than Int128 are currently unsupported")
end

# -----------------------------------------------------------------------------
#
# Arithmetic
#
# -----------------------------------------------------------------------------
zero(F::Type{<:PrimeField}) = F(Reduced(), zero(inttype(F)))
one( F::Type{<:PrimeField}) = F(Reduced(),  one(inttype(F)))

function +(a::F, b::F) where F<:PrimeField
    iszero(a) && return +b
    iszero(b) && return +a
    if char(F) > div(typemax(inttype(F)), 2)
        m, overflowed = add_with_overflow(a.n, b.n)
        if overflowed
            m += 2rem(typemax(inttype(F)), char(F)) + 2
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

inv(a::F)     where F<:PrimeField = F(Reduced(), invmod(a.n, char(F)))
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
