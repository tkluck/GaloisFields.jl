import Base.Checked: add_with_overflow, sub_with_overflow

"""
    PrimeField{I<:Integer, p}

A type representing an element in ℤ/pℤ.
"""
struct PrimeField{I<:Integer, p} <: AbstractGaloisField
    n::I
    PrimeField{I, p}(n::Integer) where {I, p} = new(mod(n, p))
    PrimeField{I, p}(::Reduced, n::Integer) where {I, p} = new(n)
end

char(::Type{PrimeField{I,p}})    where {I, p} = p
n(::Type{PrimeField{I,p}})       where {I, p} = 1
inttype(::Type{PrimeField{I,p}}) where {I, p} = I

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
    if char(F) > div(typemax(inttype(F)), 2)
        m, overflowed = add_with_overflow(a.n, b.n)
        if overflowed
            m += 2rem(typemax(inttype(F)), char(F)) + 2
        end
        F(m)
    else
        F(a.n + b.n)
    end
end

function -(a::F, b::F) where F<:PrimeField
    if char(F) > div(typemax(inttype(F)), 2)
        m, overflowed = sub_with_overflow(a.n, b.n)
        if overflowed
            m -= 2rem(typemax(inttype(F)), char(F)) + 2
        end
        F(m)
    else
        F(a.n - b.n)
    end
end

+(a::PrimeField) = a
-(a::PrimeField) = a.n == 0 ? a : typeof(a)(Reduced(), char(typeof(a)) - a.n)

*(a::F, b::F) where F<:PrimeField = F(Base.widemul(a.n, b.n))

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

(::Type{F})(n::F) where F<:PrimeField = F(n.n)
convert(::Type{F}, n::F) where F<:PrimeField = n
