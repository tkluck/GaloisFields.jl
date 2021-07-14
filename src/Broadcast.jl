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
BroadcastStyle(s::FusedModStyle, t::BroadcastStyle) = result_style(style(s), t)

const FusedModBroadcasted{F, InnerStyle} = Broadcasted{<:FusedModStyle{F, InnerStyle}}

# -----------------------------------------------------------------------------
#
# Override instantiate to use deferred reduction operations
#
# -----------------------------------------------------------------------------

instantiate(bc::FusedModBroadcasted) = instantiate(maybe_unreducedbroadcast(bc))

# -----------------------------------------------------------------------------
#
# Helper types and methods to compute the fused mod broadcast
#
# -----------------------------------------------------------------------------

"""
    struct UnreducedBroadcast

Represents a broadcast whose eltype should be F, but which has delayed
the mod operation and so for now has eltype an Integer.
The instantiate(...) function can be used to convert it to a Broadcast
object with eltype equal to F.
"""
struct UnreducedBroadcast{F <: PrimeField, BC <: Broadcasted, Bounds}
    bc :: BC
end

eltype(::UnreducedBroadcast{F}) where F <: PrimeField = F
@inline UnreducedBroadcast(F, bounds, bc) = UnreducedBroadcast{F, typeof(bc), bounds}(bc)

const TupleOf{F, N} = NTuple{N, F}

bounds(::Type{<:UnreducedBroadcast{F, BC, Bounds}}) where {F, BC, Bounds} = Bounds
bounds(a::Type{<:AbstractArray{<:PrimeField{I}}}) where I = I(0) : I(char(eltype(a)) - 1)
bounds(::Type{<:TupleOf{F}}) where F <: PrimeField{I} where I = I(0) : I(char(F) - 1)
bounds(a::Type{<:PrimeField{I}}) where I = I(0) : I(char(a) - 1)

"""
    struct Widener{I <: Integer}

A callable object that traverses a tree of Broadcasted operations
and converts the arguments to the integer type I.
"""
struct Widener{I <: Integer}
end

eltype(::Widener{I}) where I = I

(w::Widener)(a::PrimeField) = a.n % eltype(w)
(w::Widener)(a::Integer) = a % eltype(w)
(w::Widener)(a::AbstractArray{<:PrimeField}) = broadcasted(w, a)
(w::Widener)(a::AbstractArray{<:Integer}) = broadcasted(w, a)
(w::Widener)(a::TupleOf{<:PrimeField}) = map(w, a)
(w::Widener)(a::TupleOf{<:Integer}) = map(w, a)
(w::Widener)(ubc::UnreducedBroadcast) = UnreducedBroadcast(eltype(ubc), bounds(typeof(ubc)), w(ubc.bc))
(w::Widener)(bc::Broadcasted) = begin
    B = Broadcasted{typeof(BroadcastStyle(typeof(bc)))}
    if bc.f isa FusableOps
        args = map(w, bc.args)
        return B(bc.f, args, bc.axes)
    elseif bc.f isa Widener
        w′ = Widener{promote_type(eltype(w), eltype(bc.f))}()
        return B(w′, bc.args, bc.axes)
    else
        return B(w, (bc,), nothing)
    end
end

widenleavestuple(I) = ()
widenleavestuple(I, arg, args...) = (Widener{I}()(arg), widenleavestuple(I, args...)...)

intvals(a::Integer) = a
intvals(a::PrimeField) = a.n
intvals(a::AbstractArray{<:PrimeField}) = reinterpret(inttype(eltype(a)), a)
intvals(a::AbstractArray{<:Integer}) = a
intvals(a::TupleOf{<:PrimeField}) = map(intvals, a)
intvals(a::TupleOf{<:Integer}) = a
intvals(ubc::UnreducedBroadcast) = ubc.bc

intvalstuple() = ()
intvalstuple(arg, args...) = (intvals(arg), intvalstuple(args...)...)

maybe_unreducedbroadcast(a) = a
maybe_unreducedbroadcast(bc::FusedModBroadcasted{F}) where F <: PrimeField = _maybe_unreducedbroadcast(F, bc.f, style(BroadcastStyle(typeof(bc))), bc.axes, map(maybe_unreducedbroadcast, bc.args)...)

_maybe_unreducedbroadcast(F, f, innerstyle, axes, args...) = Broadcasted{typeof(innerstyle)}(f, maybe_reduce.(args), axes)

notnarrower(I, J) = promote_type(I, J) == I
iswider(I, J) = !notnarrower(J, I)

# the kinds of arguments we may find for a mergeable broadcast
const FusableOps = Union{typeof(+), typeof(-), typeof(*), typeof(^)}
const BroadcastableWithBounds{F} = Union{UnreducedBroadcast{F}, AbstractArray{F}, TupleOf{F}, F}

function _maybe_unreducedbroadcast(::Type{F}, f::FusableOps, innerstyle, axes, args::BroadcastableWithBounds{F}...) where F <: PrimeField
    resultbounds = joinbounds(f, map(bounds ∘ typeof, args)...)
    I = eltype(resultbounds)
    J = promote_type(map(eltype ∘ bounds ∘ typeof, args)...)
    # The integer operation won't overflow because J is at least as big as I
    if notnarrower(J, I)
        intargs = intvalstuple(args...)
        bci = Broadcasted{typeof(innerstyle)}(f, intargs, axes)
        return UnreducedBroadcast(F, resultbounds, bci)
    # It might overflow, but in any case it stays a bits type. We widen the arguments
    # before doing the operation to prevent the overflow.
    elseif notnarrower(Int, I)
        w = Widener{I}()
        intargs = intvalstuple(args...)
        bci = Broadcasted{typeof(innerstyle)}(f, intargs, axes)
        wbci = w(bci)
        return UnreducedBroadcast(F, resultbounds, wbci)
    else
        return Broadcasted{typeof(innerstyle)}(f, maybe_reduce.(args), axes)
    end
end

instantiate(ubc::UnreducedBroadcast) = broadcasted(eltype(ubc), intvals(ubc))
maybe_reduce(x) = x
maybe_reduce(ubc::UnreducedBroadcast) = instantiate(ubc)

# -----------------------------------------------------------------------------
#
# Transform / and // to a version that also works when reinterpreting as integers
#
# -----------------------------------------------------------------------------
const DivOp = Union{typeof(/), typeof(//)}
function _maybe_unreducedbroadcast(::Type{F}, f::DivOp, innerstyle, axes, a::BroadcastableWithBounds{F}, b::BroadcastableWithBounds{F}) where F <: PrimeField
    b_inverse = _maybe_unreducedbroadcast(F, inv, innerstyle, axes, b)
    return _maybe_unreducedbroadcast(F, *, innerstyle, axes, a, b_inverse)
end

function _maybe_unreducedbroadcast(::Type{F}, f::typeof(inv), innerstyle, axes, a::BroadcastableWithBounds{F}) where F <: PrimeField
    bc_inverse = broadcasted(invmod, intvals(a), char(F))
    return UnreducedBroadcast(F, bounds(F), bc_inverse)
end

# -----------------------------------------------------------------------------
#
# Transform == to a version that does only a single mod reduction
#
# -----------------------------------------------------------------------------
function _maybe_unreducedbroadcast(::Type{F}, f::typeof(==), innerstyle, axes, a::BroadcastableWithBounds{F}, b::BroadcastableWithBounds{F}) where F <: PrimeField
    bc_diff = _maybe_unreducedbroadcast(F, -, innerstyle, axes, a, b)
    return broadcasted(iszero, maybe_reduce(bc_diff))
end


# -----------------------------------------------------------------------------
#
# Testability
#
# -----------------------------------------------------------------------------
_is_integer_broadcast(leaf) = eltype(leaf) <: Integer
_is_integer_broadcast(bc::Broadcasted) = all(_is_integer_broadcast, bc.args) || bc.f isa Widener
_is_integer_broadcast(ubc::UnreducedBroadcast) = _is_integer_broadcast(ubc.bc)
isfused(bc) = false
isfused(bc::FusedModBroadcasted) = _is_integer_broadcast(maybe_unreducedbroadcast(bc))

broadcasted_calls(x) = x
broadcasted_calls(expr::Expr) = if expr.head == :call
    Expr(:call, broadcasted, broadcasted_calls.(expr.args)...)
elseif expr.head == :$
    expr.args[1]
else
    Expr(expr.head, broadcasted_calls.(expr.args)...)
end
macro broadcasted(expr)
    expr = broadcasted_calls(expr)
    esc(expr)
end

end # module
