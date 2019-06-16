@doc read(open(joinpath(@__DIR__, "..", "README.md")), String)
module GaloisFields

import LinearAlgebra: norm, tr
import Random: AbstractRNG, SamplerType
import Serialization: deserialize

import Polynomials: Poly, coeffs
import Primes: factor, Factorization

# imports for overloading
import Base: zero, one, +, -, *, /, //, ^, inv, iszero
import Base: show
import Base: convert, promote_rule, promote_type, eltype
import Base: iterate
import Base: rand

"""
    abstract type AbstractGaloisField <: Number end

A type representing finite fields.
"""
abstract type AbstractGaloisField <: Number end

"""
    abstract type AbstractExtensionField <: AbstractGaloisField end

A type representing a finite extension of an underlying finite field.
"""
abstract type AbstractExtensionField <: AbstractGaloisField end

"""
    Reduced()

A helper singleton used for asserting that an input value
has already been reduced mod p.
"""
struct Reduced end

"""
    NonNegative()

A helper singleton used for asserting that an input value
is non-negative.

This means we can use `rem(...)` instead of `mod(...)`, which
is significantly faster.
"""
struct NonNegative end

"""
    p = char(GaloisField(3)) # returns 3

Return the characteristic of a finite field, or 0 for <:Integer or <:Rational{<Integer}.
"""
char(::Type{<:Rational{<:Integer}}) = 0
char(::Type{<:Integer}) = 0

include("BoundedIntegers.jl")
include("PrimeFields.jl")
include("ZechLog.jl")
include("ExtensionFields.jl")
include("BinaryFields.jl")
include("Conversions.jl")
include("Iterations.jl")
include("Reinterpret.jl")
include("Broadcast.jl")
include("Display.jl")
include("LinearAlgebra.jl")

"""
    const F = GaloisField(p)
    const F,α = GaloisField(p, :β => [1, 0, 1])
    const F,α = GaloisField(p, n, :β)

Return a type representing a finite field.

The single-argument signature returns the finite field ``ℤ/pℤ``.

The two-arguments signature returns an algebraic extension of that field,
with minimum polynomial given by the second argument: a dense representation
of the univariate, monic polynomial, with ascending degree.

The three-arguments signature returns an algebraic extension of that field,
with minimum polynomial equal to the
[Conway polynomial](https://en.wikipedia.org/wiki/Conway_polynomial_(finite_fields))
 for ``(p,n)``. The `GaloisFields` package ships with a database of Conway
 polynomials and will raise an error if it does not contain an entry for
 ``(p,n)``.

Note that in the latter two cases, the variable name (e.g. β above) is part of the
type. This lets you define identifications between isomorphic (sub)fields. For
example, with the following definition

    const F = @GaloisField! 𝔽₂ β^2 + β + 1
    const G = @GaloisField! 𝔽₂ γ^2 + γ + 1

the fields ``F`` and ``G`` are isomorphic, but not canonically. We might
define

    @GaloisFields.identify β => γ + 1
    @GaloisFields.identify γ => β + 1

to allow for conversions like

    G(β)
    convert(F, γ + 1)

In the Conway case, you do not have to define your own identifications, as
the Conway polynomials satisfy compatibility relations that allow us to use
certain distinguished inclusions between them.
"""
function GaloisField end

GaloisField(p::Integer, minpoly::Poly) = GaloisField(GaloisField(p, 1), minpoly)
GaloisField(p::Integer, minpoly::Pair) = GaloisField(GaloisField(p, 1), minpoly)
GaloisField(F::Type{<:AbstractGaloisField}, minpoly::Poly) = GaloisField(F, minpoly.var => coeffs(minpoly))
function GaloisField(F::Type{<:AbstractGaloisField}, minpoly::Pair{Symbol, <:AbstractVector{<:Number}}, conway=false)
    sym, coeffs = minpoly
    mp = tuple(map(F, coeffs)...)
    N = length(coeffs) - 1
    if char(F) == 2 && F <: PrimeField
        I = if N <= 8
            UInt8
        elseif N <= 16
            UInt16
        elseif N <= 32
            UInt32
        elseif N <= 64
            UInt64
        else
            # fall through to ExtensionField
            nothing
        end
        if I !== nothing
            minpolymask = zero(I)
            for (i, c) in enumerate(mp)
                minpolymask |= (c.n % I) << (i - 1)
            end
            BF = BinaryField{I, N, sym, minpolymask, conway}
            return BF, gen(BF)
        end
    end
    EF = ExtensionField{F, N, sym, mp, conway}
    return EF, gen(EF)
end

const DATAFILE = joinpath(@__DIR__, "..", "deps", "conwaypolynomials.data")
_conwaypolynomials = nothing
function conwaypolynomial(p::Integer, n::Integer)
    global _conwaypolynomials
    if _conwaypolynomials == nothing
        _conwaypolynomials = deserialize(open(DATAFILE))
    end
    return _conwaypolynomials[p, n]
end

GaloisField(q::Integer) = GaloisField(factor(q))
GaloisField(q::Integer, sym::Symbol) = GaloisField(factor(q), sym)
GaloisField(p::Integer, n::Integer) = GaloisField(p, n, gensym())
GaloisField(factors::Factorization) = GaloisField(factors, gensym())

function GaloisField(p::Integer, n::Integer, sym::Symbol)
    I = inttype(p)
    # standardize on what type of integer we use in the type
    # parameter. This allows us to just write e.g. PrimeField{I, 2} where I
    # instead of PrimeField{I, Int8(2)}.
    J = p <= typemax(Int) ? Int : I
    𝔽ₚ = PrimeField{inttype(p), J(p)}
    if n == 1
        return 𝔽ₚ
    else
        coeffs = conwaypolynomial(p, n)
        return GaloisField(𝔽ₚ, sym => map(𝔽ₚ, coeffs), true)
    end
end

function GaloisField(factors::Factorization, sym::Symbol)
    if length(factors) != 1
        error("There is no finite field of order $(prod(f))")
    end
    (p, n), = factors
    return GaloisField(p, n, sym)
end


"""
    F = @GaloisField 3
    F = @GaloisField ℤ/3ℤ
    F = @GaloisField 𝔽₃
    F,α = @GaloisField 3^2

Different ways of defining a finite field.

For defining a Galois field of prime power order with ``n>1``, consider using
`@GaloisField!` instead, which allows specifying the display name of the
generator.
"""
macro GaloisField(expr)
    res = _parse_declaration(expr)
    if res === nothing
        error("Not implemented: @GaloisField $expr")
    end
    return GaloisField(res)
end

function _factorization(p::Integer, n::Integer)
    return Factorization{Int}(Dict(p => n))
end

function _parse_declaration(expr)
    # @GaloisField q
    if expr isa Integer
        return expr
    elseif expr isa Expr
        # @GaloisField ℤ/pℤ
        if expr.head == :call && expr.args[1] == :/ &&
            expr.args[2] == :ℤ && expr.args[3].head == :call &&
            expr.args[3].args[1] == :* && expr.args[3].args[3] == :ℤ
            p = expr.args[3].args[2]
            factors = _factorization(p, 1)
            return factors
        elseif expr.head == :call && expr.args[1] == :^ &&
            expr.args[2] isa Integer && expr.args[3] isa Integer
            p, n = expr.args[2:end]
            factors = _factorization(p, n)
            return factors
        end
    # @GaloisField 𝔽₃₇
    elseif expr isa Symbol
        str = collect(string(expr))
        if str[1] == '𝔽'
            s = ['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉']
            indices = indexin(str[2:end], s)
            q = 0
            for ix in indices
                q = 10q + ix - 1
            end
            return q
        end
    end
    return nothing
end

_parsepoly(x) = x
_parsepoly(x::Symbol) = :( $(Poly([0, 1], x) ) )
function _parsepoly(expr::Expr)
    Expr(expr.head, expr.args[1], map(_parsepoly, expr.args[2:end])...)
end

"""
    const G = @GaloisField! 3 β^2 + 1
    const G = @GaloisField! 𝔽₃ β^2 + 1
    const G = @GaloisField! 3^2 β
    const K = GaloisField(3)
    const G = @GaloisField! K β^2 + 1

Define a finite field `G` and inject a variable for its
primitive element into the current scope.

Note that the variable name (e.g. β above) is part of the type. This lets you
define identifications between isomorphic (sub)fields. For example, with the
following definition

    const F = @GaloisField! 𝔽₂ β^2 + β + 1
    const G = @GaloisField! 𝔽₂ γ^2 + γ + 1

the fields ``F`` and ``G`` are isomorphic, but not canonically. We might
define

    @GaloisFields.identify β => γ + 1
    @GaloisFields.identify γ => β + 1

to allow for conversions like

    G(β)
    convert(F, γ + 1)
"""
macro GaloisField!(expr, minpoly)
    if minpoly isa Symbol
        poly = QuoteNode(minpoly)
        sym = minpoly
    else
        poly = @eval $(_parsepoly(minpoly))
        sym = poly.var
    end
    decl = _parse_declaration(expr)
    F = something(decl, esc(expr))
    quote
        EF, $(esc(sym)) = $GaloisField($F, $poly)
        EF
    end
end

if VERSION < v"1.3-"
    # https://github.com/JuliaLang/julia/pull/31822
    Base.deepcopy(x::AbstractGaloisField) = x
end

export GaloisField, @GaloisField, @GaloisField!, char
export norm, tr

end
