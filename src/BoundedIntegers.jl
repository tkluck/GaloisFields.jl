module BoundedIntegers

import Base: +, -, *, ^, div, rem, divrem, mod
import Base: typemin, typemax, promote_rule, convert
import Base: zero, one, iszero, isone

import ..Util: widen_bits

"""
    BoundedInteger{Bounds, I <: Integer}

An integer type whose values are of type I and guaranteed to lie inside Bounds.
This is useful for making compile-time decisions about widening I.
"""
struct BoundedInteger{Bounds, I <: Integer} <: Integer
    n::I
    BoundedInteger{Bounds, I}(n::Integer) where {Bounds, I <: Integer} = new(n) #(@assert(n in Bounds); new(n))
end

bounds(::Type{BoundedInteger{Bounds, I}}) where {Bounds, I} = Bounds
inttype(::Type{BoundedInteger{Bounds, I}}) where {Bounds, I} = I
typemin(T::Type{<:BoundedInteger}) = first(bounds(T))
typemax(T::Type{<:BoundedInteger}) = last(bounds(T))
bounds(n::BoundedInteger) = bounds(typeof(n))
inttype(n::BoundedInteger) = inttype(typeof(n))
val(n::BoundedInteger) = n.n

BoundedInteger(n::Integer) = BoundedInteger{n:n, typeof(n)}(n)
BoundedInteger{Bounds}(n::Integer) where Bounds = BoundedInteger{Bounds, typeof(n)}(n)

zero(::Type{BoundedInteger{Bounds, I}}) where {Bounds, I} = (@assert(zero(I) in Bounds); BoundedInteger{Bounds, I}(zero(I)))
one(::Type{BoundedInteger{Bounds, I}}) where {Bounds, I} = (@assert(one(I) in Bounds); BoundedInteger{Bounds, I}(one(I)))
iszero(a::BoundedInteger) = iszero(val(a))
isone(a::BoundedInteger) = isone(val(a))

_joinbounds(op, b) = b
_joinbounds(op, b, c...) = _joinbounds(op, b, _joinbounds(op, c...))
function _joinbounds(op, b, c)
    extrema = (
        op(widen_bits(first(b)), widen_bits(first(c))),
        op(widen_bits(first(b)), widen_bits(last(c) )),
        op(widen_bits(last(b)),  widen_bits(first(c))),
        op(widen_bits(last(b)),  widen_bits(last(c) )),
    )
    lo, hi = min(extrema...), max(extrema...)
    I = _mintype(lo : hi)
    I(lo) : I(hi)
end

@inline function _mintype(bounds)
    T = Int8
    while isbitstype(T)
        first(bounds) >= typemin(T) && last(bounds) <= typemax(T) && return T
        T = widen_bits(T)
    end
    error("BoundedInteger out of bounds: $bounds")
end

@inline function _apply(op, x::BoundedInteger...)
    resultbounds = _joinbounds(op, map(bounds, x)...)
    I = promote_type(map(inttype, x)...)
    J = _mintype(resultbounds)
    K = promote_type(I, J)
    T = BoundedInteger{resultbounds, K}
    T(op(map(K, map(val, x))...))
end

for convexop in [:+, :-, :*, :div]
    @eval $convexop(x::BoundedInteger...) = _apply($convexop, x...)
end

@inline ^(a::BoundedInteger, b::Integer) = _apply(^, a, BoundedInteger(b))
@inline div(a::BoundedInteger, b::Integer) = _apply(div, a, BoundedInteger(b))
@inline rem(a::BoundedInteger, b::Integer) = BoundedInteger{-abs(b-1):abs(b-1)}(rem(val(a), b))
@inline mod(a::BoundedInteger, b::Integer) = b > 0 ?
                                    BoundedInteger{0:(b-1)}(mod(val(a), b)) :
                                    BoundedInteger{(b+1):0}(mod(val(a), b))
@inline divrem(a::BoundedInteger, b::Integer) = div(a, b), rem(a, b)

@inline rem(a::BoundedInteger, I::Type{<:Integer}) = val(a) % I

promote_rule(T::Type{<:BoundedInteger}, S::Type) = promote_type(inttype(T), S)
convert(T::Type{<:Number}, a::BoundedInteger) = convert(T, val(a))

function Base.show(io::IO, a::BoundedInteger)
    print(io, "BoundedInteger{$(bounds(a))}($(val(a)))")
end

end # module
