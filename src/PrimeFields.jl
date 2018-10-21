
"""
    PrimeField{I<:Integer, p}

A type representing an element in ℤ/pℤ.
"""
struct PrimeField{I<:Integer, p} <: AbstractGaloisField
    n::I
    PrimeField{I, p}(n::I) where {I, p} = new(mod(n, p))
    PrimeField{I, p}(::Reduced, n::I) where {I, p} = new(n)
end

char(::Type{PrimeField{I,p}}) where {I, p} = p
eltype(::Type{PrimeField{I,p}}) where {I, p} = I

# -----------------------------------------------------------------------------
#
# Arithmetic
#
# -----------------------------------------------------------------------------
zero(F::Type{<:PrimeField{I}}) where I = F(zero(I))
one( F::Type{<:PrimeField{I}}) where I = F(one(I))

+(a::F, b::F) where F<:PrimeField = F(a.n + b.n)
-(a::F, b::F) where F<:PrimeField = F(a.n - b.n)

+(a::PrimeField) = a
-(a::PrimeField) = typeof(a)(Reduced(), char(typeof(a)) - a.n)

*(a::F, b::F) where F<:PrimeField = F(a.n * b.n)

inv(a::F)     where F<:PrimeField = F(Reduced(), invmod(a.n, char(F)))
/(a::F,b::F)  where F<:PrimeField = a * inv(b)
//(a::F,b::F) where F<:PrimeField = a * inv(b)

iszero(a::PrimeField) = iszero(a.n)

show(io::IO, a::PrimeField) = show(io, a.n)
function show(io::IO, ::Type{PrimeField{I,p}}) where {I,p}
    number = replace("$p", r"[0-9]" => x->['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉'][parse(Int,x) + 1])
    write(io, "𝔽$number")
end

promote_rule(F::Type{<:PrimeField}, ::Type{<:Integer}) = F
convert(F::Type{PrimeField{I,p}}, i::Integer) where {I,p} = F(Reduced(), I(mod(i, p)))

(::Type{F})(n::F) where F<:PrimeField = F(n.n)
convert(::Type{F}, n::F) where F<:PrimeField = n
