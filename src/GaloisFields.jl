"""
    using GaloisFields

A module for finite fields. Synopsis:

    using GaloisFields
    F = GaloisField(3)
    F = @GaloisField ℤ/3ℤ
    F = @GaloisField 𝔽₃

    F, β = GaloisField(3, :β => [2, 1, 1])
    F = @GaloisField! 𝔽₃ β^2 + β + 2

See the docstrings for `GaloisField`, `@GaloisField`, and `@GaloisField!` for details.
"""
module GaloisFields

using Polynomials: Poly, coeffs

# imports for overloading
import Base: zero, one, +, -, *, /, //, inv, iszero
import Base: show
import Base: convert, promote_rule, promote_type, eltype
import Base: iterate

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
char(::Type{<:Rational{<:Integer}}) = 0
char(::Type{<:Integer}) = 0

include("PrimeFields.jl")
include("ExtensionFields.jl")
include("Conversions.jl")
include("Iterations.jl")
include("Display.jl")

"""
    F = GaloisField(p)
    F,α = GaloisField(p, :β => [1, 0, 1])

Return a type representing a finite field.

The single-argument signature returns the finite field ``ℤ/pℤ``.

The two-arguments signature returns an algebraic extension of that field,
with minimum polynomial given by the second argument: a dense representation
of the univariate, monic polynomial, with ascending degree.

Note that in the latter case, the variable name (e.g. β above) is part of the
type. This lets you define identifications between isomorphic (sub)fields. For
example, with the following definition

    F = @GaloisField! 𝔽₂ β^2 + β + 1
    G = @GaloisField! 𝔽₂ γ^2 + γ + 1

the fields ``F`` and ``G`` are isomorphic, but not canonically. We might
define

    @GaloisFields.identify β => γ + 1
    @GaloisFields.identify γ => β + 1

to allow for conversions like

    G(β)
    convert(F, γ + 1)
"""
GaloisField(p::Integer) = PrimeField{typeof(p), p}
GaloisField(p::Integer, args...) = GaloisField(GaloisField(p), args...)
GaloisField(F::Type{<:AbstractGaloisField}, minpoly::Poly) = GaloisField(F, minpoly.var => coeffs(minpoly))
function GaloisField(F::Type{<:AbstractGaloisField}, minpoly::Pair{Symbol, <:AbstractVector{<:Number}})
    sym, coeffs = minpoly
    mp = tuple(map(F, coeffs)...)
    N = length(coeffs) - 1
    EF = ExtensionField{F, N, sym, mp}
    return EF, gen(EF)
end

"""
    F = @GaloisField 3
    F = @GaloisField ℤ/3ℤ
    F = @GaloisField 𝔽₃

Different ways of defining a finite field of a given order.
"""
macro GaloisField(expr)
    res = _parse_declaration(expr)
    if res === nothing
        throw("Not implemented: @GaloisField $expr")
    end
    return res
end

function _parse_declaration(expr)
    # @GaloisField p
    if expr isa Integer
        return :( $GaloisField($expr) )
    elseif expr isa Expr
        # @GaloisField ℤ/pℤ
        if expr.head == :call && expr.args[1] == :/ &&
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
    return nothing
end

_parsepoly(x) = x
_parsepoly(x::Symbol) = :( $(Poly([0, 1], x) ) )
function _parsepoly(expr::Expr)
    Expr(expr.head, expr.args[1], map(_parsepoly, expr.args[2:end])...)
end

"""
    G = @GaloisField! 3 β^2 + 1
    G = @GaloisField! 𝔽₃ β^2 + 1
    K = GaloisField(3)
    G = @GaloisField! K β^2 + 1

Define a finite field `G` and inject a variable for its
primitive element into the current scope.

Note that the variable name (e.g. β above) is part of the type. This lets you
define identifications between isomorphic (sub)fields. For example, with the
following definition

    F = @GaloisField! 𝔽₂ β^2 + β + 1
    G = @GaloisField! 𝔽₂ γ^2 + γ + 1

the fields ``F`` and ``G`` are isomorphic, but not canonically. We might
define

    @GaloisFields.identify β => γ + 1
    @GaloisFields.identify γ => β + 1

to allow for conversions like

    G(β)
    convert(F, γ + 1)
"""
macro GaloisField!(expr, minpoly)
    poly = @eval $(_parsepoly(minpoly))
    decl = _parse_declaration(expr)
    F = something(decl, esc(expr))
    quote
        EF, $(esc(poly.var)) = $GaloisField($F, $poly)
        EF
    end
end

export GaloisField, @GaloisField, @GaloisField!, char

end
