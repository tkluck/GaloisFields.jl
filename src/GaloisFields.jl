"""
    using GaloisFields

A module for finite fields. Synopsis:

    using GaloisFields
    # a finite field of prime order
    F = GaloisField(3)
    # a finite field of prime power order, plus a primitive element
    F,β = GaloisField(3, 2)
    # friendlier syntax for the same thing
    F,β = @GaloisField 3^2
    # supply your own minimum polynomial for a primitive element
    F,β = GaloisField(3, :β => [2, 1, 1])
    # friendlier syntax for the same thing; inject β
    F = @GaloisField! 3 β^2 + β + 2

    # friendly syntax
    F = @GaloisField ℤ/3ℤ
    F = @GaloisField 𝔽₃

    # it is also possible to declare a Galois field as an extension
    # of a field you defined earlier:

    # a finite field of degree 2 over F
    # (γ generates it over F, not necessarily over the prime field!)
    G,γ = GaloisField(F, 2)
    # friendlier syntax
    G,γ = @GaloisField F^2
    # supply your own minimum polynomial for a primitive element
    G,γ = GaloisField(F, :γ => [2, 1, 1])
    # friendlier syntax for the same thing; inject γ
    F = @GaloisField! F γ^2 + γ + 2

    # friendly syntax
    F = @GaloisField! ℤ/3ℤ γ^2 + 1
    F = @GaloisField! 𝔽₃    γ^2 + 1

In all cases, the variable name (e.g. β or γ) is part of the type; this
lets you define identifications between isomorphic (sub)fields. For example,
with the following definition

    F = @GaloisField! 𝔽₂ β^2 + β + 1
    G = @GaloisField! 𝔽₂ γ^2 + γ + 1

the fields ``F`` and ``G`` are isomorphic, but not canonically. We might
define

    @GaloisFields.identify β => γ + 1
    @GaloisFields.identify γ => β + 1

to allow for conversions like

    G(β)
    convert(F, γ + 1)

This module has a special case for efficient binary representation of
power-of-two fields.
"""
module GaloisFields

using Polynomials: Poly, coeffs

# imports for overloading
import Base: zero, one, +, -, *, /, //, inv
import Base: show
import Base: convert, promote_rule, promote_type, eltype

"""
    abstract type AbstractGaloisField <: Number end

A type representing finite fields.
"""
abstract type AbstractGaloisField <: Number end

"""
    Reduced()

A helper singleton used for asserting that an input value
has already been reduced mod p.
"""
struct Reduced end

"""
    p = char(GaloisField(3)) # returns 3

Return the characteristic of a finite field, or 0 for <:Integer or <:Rational{<Integer}.
"""
char(x) = char(typeof(x))
char(::Type{<:Rational{<:Integer}}) = 0
char(::Type{<:Integer}) = 0

include("PrimeFields.jl")
include("ExtensionFields.jl")
include("Conversions.jl")

"""
    F = GaloisField(p)
    F,α = GaloisField(p, [1, 0, 1])

Return a type representing a finite field.

The single-argument signature returns the finite field ``ℤ/pℤ``.

The two-arguments signature returns an algebraic extension of that field,
with minimum polynomial given by the second argument: a dense representation
of the univariate, monic polynomial, with ascending degree.
"""
GaloisField(p::Integer) = PrimeField{typeof(p), p}
GaloisField(p::Integer, args...) = GaloisField(GaloisField(p), args...)
GaloisField(F::Type{<:PrimeField}, minpoly::Poly) = GaloisField(F, minpoly.var => coeffs(minpoly))
function GaloisField(F::Type{<:PrimeField}, minpoly::Pair{Symbol, <:AbstractVector{<:Number}})
    sym, coeffs = minpoly
    mp = tuple(map(F, coeffs)...)
    N = length(coeffs) - 1
    EF = ExtensionField{F, N, sym, mp}
    return EF, gen(EF)
end

macro GaloisField(expr)
    # @GaloisField p
    if expr isa Integer
        return :( $GaloisField($expr) )
    elseif expr isa Expr
        # @GaloisField p^n
        if expr.head == :call && expr.args[1] == :^
            p, n = expr.args[2:end]
            return :( $GaloisField($p, $n) )
        # @GaloisField ℤ/pℤ
        elseif expr.head == :call && expr.args[1] == :/ &&
            expr.args[2] == :ℤ && expr.args[3].head == :call &&
            expr.args[3].args[1] == :* && expr.args[3].args[3] == :ℤ
            p = expr.args[3].args[2]
            return :( $GaloisField($p) )
        end
    # @GaloisField 𝔽₃₇
    elseif expr isa Symbol
        str = collect(string(expr))
        if str[1] == '𝔽'
            s = ['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉']
            indices = indexin(str[2:end], s)
            p = 0
            for ix in indices
                p = 10p + ix - 1
            end
            return :( $GaloisField($p) )
        end
    end
    throw("Not implemented: @GaloisField $expr")
end

parsepoly(x) = x
parsepoly(x::Symbol) = :( $(Poly([0, 1], x) ) )
function parsepoly(expr::Expr)
    Expr(expr.head, expr.args[1], map(parsepoly, expr.args[2:end])...)
end

macro GaloisField!(expr, minpoly)
    poly = @eval $(parsepoly(minpoly))
    quote
        F = @GaloisField $expr
        EF, $(esc(poly.var)) = $GaloisField(F, $poly)
        EF
    end
end


export GaloisField, @GaloisField, @GaloisField!, char

end
