"""
    GaloisFields.Broadcast

This module optimizes operations for Vector{<:PrimeField} by doing only a single
mod(..., char(F)) operation for a fused broadcasted operation (e.g. `mod(x.n +
y.n - z.n, char(F))` instead of `mod(mod(x.n + y.n, char(F)) - z.n, char(F))`.

This is applied for broadcasted combinations of the `+, -, *` operations over
(vectors of) `PrimeFields`.

It is possible we need to widen to a larger integer type for the intermediate
computation. The utility function `bounds` can compute the required type
at compile time.

Implementation note: it is very important that the functions are inferable.
Two patterns that we use to ensure that, are:

1) Instead of `map(a -> widenleaves(I, a), args)`, we use

    widenleavestuple(I) = ()
    widenleavestuple(I, arg, args...) = (widenleaves(I, arg), widenleavestuple(I, args...)...)

This ensures that `I` does not need to be captured into a closure, and also
that the result type is completely inferable.

2) We user `broadcast((a_i, I) -> a_i % I, a, I)` instead of `map(a_i -> a_i %I, a)`
also to prevent having to capture in a closure.
"""
module Broadcast

import Base: eltype
import Base.Broadcast: Broadcasted, BroadcastStyle, broadcasted
import Base.Broadcast: Style, DefaultArrayStyle, AbstractArrayStyle
import Base.Broadcast: result_style, instantiate

import ..Util: joinbounds
import GaloisFields: inttype, char, PrimeField, posmod

# -----------------------------------------------------------------------------
#
# Define FusedModStyle broadcasting style and declare when it should be used.
#
# -----------------------------------------------------------------------------
struct FusedModStyle{F, Style} <: BroadcastStyle end
eltype(::FusedModStyle{F, Style}) where {F, Style} = F
style(::FusedModStyle{F, Style})  where {F, Style} = Style()

function BroadcastStyle(T::Type{<:F}) where F <: PrimeField
    st = DefaultArrayStyle{0}()
    FusedModStyle{F, typeof(st)}()
end

function BroadcastStyle(T::Type{<:AbstractArray{F}}) where F <: PrimeField
    st = DefaultArrayStyle{ndims(T)}()
    FusedModStyle{F, typeof(st)}()
end
function BroadcastStyle(T::Type{<:NTuple{N, F}}) where {N, F <: PrimeField}
    st = Style{Tuple}()
    FusedModStyle{F, typeof(st)}()
end
function BroadcastStyle(s::FusedModStyle{F}, t::FusedModStyle{F}) where F <: PrimeField
    st = result_style(style(s), style(t))
    FusedModStyle{F, typeof(st)}()
end
function BroadcastStyle(s::FusedModStyle, t::AbstractArrayStyle{0})
    st = result_style(style(s), t)
    FusedModStyle{eltype(s), typeof(st)}()
end
function BroadcastStyle(s::AbstractArrayStyle{0}, t::FusedModStyle)
    st = result_style(s, style(t))
    FusedModStyle{eltype(t), typeof(st)}()
end
BroadcastStyle(s::FusedModStyle, t::FusedModStyle) = result_style(style(s), style(t))
BroadcastStyle(s::FusedModStyle, t::BroadcastStyle) = t

const FusedModBroadcasted{F, InnerStyle} = Broadcasted{<:FusedModStyle{F, InnerStyle}}

# -----------------------------------------------------------------------------
#
# Override instantiate to use lazy reduction operations
#
# -----------------------------------------------------------------------------

instantiate(bc::Broadcasted{<:FusedModStyle{F}}) where F <: PrimeField = instantiate(fieldvals(F, unreducedbroadcast(bc)))

# -----------------------------------------------------------------------------
#
# Helper types and methods to compute the fused mod broadcast
#
# -----------------------------------------------------------------------------

struct UnreducedBroadcast{F <: PrimeField, BC <: Broadcasted, Bounds}
    bc :: BC
end

@inline UnreducedBroadcast(F, bounds, bc) = UnreducedBroadcast{F, typeof(bc), bounds}(bc)

const TupleOf{F, N} = NTuple{N, F}
# the kinds of arguments we may find for a mergeable broadcast
const BroadcastableWithBounds{F} = Union{UnreducedBroadcast{F}, AbstractArray{F}, TupleOf{F}, F}

bounds(::Type{<:UnreducedBroadcast{F, BC, Bounds}}) where {F, BC, Bounds} = Bounds
bounds(a::Type{<:AbstractArray{<:PrimeField{I}}}) where I = I(0) : I(char(eltype(a)) - 1)
bounds(::Type{<:TupleOf{F}}) where F <: PrimeField{I} where I = I(0) : I(char(F) - 1)
bounds(a::Type{<:PrimeField{I}}) where I = I(0) : I(char(a) - 1)

widenleaves(I, ubc::UnreducedBroadcast) = widenleaves(I, ubc.bc)
widenleaves(I, a::AbstractArray{<:PrimeField}) = broadcasted(widenleaves, I, a)
widenleaves(I, a::AbstractArray{<:Integer}) = broadcasted(widenleaves, I, a)
widenleaves(I, a::TupleOf{<:PrimeField}) = map(a_i -> widenleaves(I, a_i), a)
widenleaves(I, a::TupleOf{<:Integer}) = map(a_i -> widenleaves(I, a_i), a)
widenleaves(I, a::PrimeField) = a.n % I
widenleaves(I, a::Integer) = a % I
function widenleaves(I, bc::Broadcasted)
    if bc.f == widenleaves
        return broadcasted(widenleaves, I, bc.args[2])
    else
        extendedleaves = widenleavestuple(I, bc.args...)
        return Broadcasted{typeof(BroadcastStyle(typeof(bc)))}(bc.f, extendedleaves, bc.axes)
    end
end

widenleavestuple(I) = ()
widenleavestuple(I, arg, args...) = (widenleaves(I, arg), widenleavestuple(I, args...)...)

intvals(ubc::UnreducedBroadcast) = ubc.bc
intvals(a::AbstractArray{<:PrimeField}) = reinterpret(inttype(eltype(a)), a)
intvals(a::AbstractArray{<:Integer}) = a
intvals(a::TupleOf{<:PrimeField}) = map(a -> a.n, a)
intvals(a::TupleOf{<:Integer}) = a
intvals(a::PrimeField) = a.n
intvals(a::Integer) = a

intvalstuple() = ()
intvalstuple(arg, args...) = (intvals(arg), intvalstuple(args...)...)

reducedvals(F, ubc::UnreducedBroadcast) = broadcasted(posmod, intvals(ubc), char(F))
reducedvals(F, a::AbstractArray{<:PrimeField}) = broadcasted(a_i -> a_i.n, a)
reducedvals(F, a::AbstractArray{<:Integer}) = a
reducedvals(F, a::TupleOf{<:PrimeField}) = map(a_i -> a_i.n, a)
reducedvals(F, a::TupleOf{<:Integer}) = a
reducedvals(F, a::PrimeField) = a.n
reducedvals(F, a::Integer) = posmod(a, char(F))

reducedvalstuple(F) = ()
reducedvalstuple(F, arg, args...) = (reducedvals(F, arg), reducedvalstuple(F, args...)...)

fieldvals(F, bc::Broadcasted) = broadcasted(F, bc)
fieldvals(F, ubc::UnreducedBroadcast) = broadcasted(F, intvals(ubc))
fieldvals(F, a::AbstractArray{<:PrimeField}) = a
fieldvals(F, a::AbstractArray{<:Integer}) = broadcasted(F, a)
fieldvals(F, a::TupleOf{<:PrimeField}) = a
fieldvals(F, a::TupleOf{<:Integer}) = map(F, a)
fieldvals(F, a::PrimeField) = a
fieldvals(F, a::Integer) = F(a)

fieldvalstuple(F) = ()
fieldvalstuple(F, arg, args...) = (fieldvals(F, arg), fieldvalstuple(F, args...)...)

unreducedbroadcast(a) = a
unreducedbroadcast(bc::FusedModBroadcasted{F}) where F <: PrimeField = unreducedbroadcast(F, bc.f, style(BroadcastStyle(typeof(bc))), bc.axes, map(unreducedbroadcast, bc.args)...)

function unreducedbroadcast(F, f, innerstyle, axes, args...)
    fieldargs = fieldvalstuple(F, args...)
    bc = Broadcasted{typeof(innerstyle)}(f, fieldargs)
    return UnreducedBroadcast(F, bounds(F), bc)
end

const FusableOps = Union{typeof(+), typeof(-), typeof(*), typeof(^), typeof(posmod)}

function unreducedbroadcast(::Type{F}, f::FusableOps, innerstyle, axes, args::BroadcastableWithBounds{F}...) where F <: PrimeField
    resultbounds = joinbounds(f, map(bounds ∘ typeof, args)...)
    I = eltype(resultbounds)
    J = promote_type(map(eltype ∘ bounds ∘ typeof, args)...)
    if promote_type(I, J) == J
        intargs = intvalstuple(args...)
        return UnreducedBroadcast(F, resultbounds, Broadcasted{typeof(innerstyle)}(f, intargs, axes))
    elseif promote_type(I, Int) == Int
        extargs = widenleavestuple(I, args...)
        return UnreducedBroadcast(F, resultbounds, Broadcasted{typeof(innerstyle)}(f, extargs, axes))
    elseif J != inttype(F)
        redargs = reducedvalstuple(F, args...)
        return unreducedbroadcast(f, innerstyle, axes, redargs...)
    else
        fieldargs = fieldvalstuple(F, args...)
        bci = Broadcasted{typeof(innerstyle)}((a -> a.n)∘f, fieldargs, axes)
        return UnreducedBroadcast(F, bounds(F), bci)
    end
end

# -----------------------------------------------------------------------------
#
# Transform / and // to a version that also works when reinterpreting as integers
#
# -----------------------------------------------------------------------------
const DivOp = Union{typeof(/), typeof(//)}
function unreducedbroadcast(::Type{F}, f::DivOp, innerstyle, axes, a::BroadcastableWithBounds{F}, b::BroadcastableWithBounds{F}) where F <: PrimeField
    inverse = unreducedbroadcast(broadcasted(invmod, intvals(b), char(F)))
    return unreducedbroadcast(F, *, innerstyle, axes, a, inverse)
end

end # module
