"""
    GaloisFields.Broadcast

Currently, this module contains optimizations for `Vector{<:PrimeField}`.

This module optimizes operations for Vector{<:PrimeField} by doing only a single
mod(..., char(F)) operation for a fused broadcasted operation (e.g. `mod(x.n +
y.n - z.n, char(F))` instead of `mod(mod(x.n + y.n, char(F)) - z.n, char(F))`.

This is applied for broadcasted combinations of the `+, -, *` operations over
vectors of `PrimeFields`.

We use a specially crafted `BoundedInteger` type to widen the intermediate
results where necessary. See `GaloisFields.BoundedIntegers`.
"""
module Broadcast

import Base: copyto!, similar, eltype
import Base.Broadcast: Broadcasted, BroadcastStyle, broadcasted
import Base.Broadcast: DefaultArrayStyle, AbstractArrayStyle
import Base.Broadcast: flatten, copyto!

import ..BoundedIntegers: BoundedInteger
import GaloisFields: inttype, char, PrimeField, posmod

struct FusedModBroadcast{F} <: BroadcastStyle end
eltype(::FusedModBroadcast{F}) where F = F

BroadcastStyle(::Type{<:AbstractArray{F}}) where F <: PrimeField = FusedModBroadcast{F}()
BroadcastStyle(::T, ::T) where T <: FusedModBroadcast = T()
BroadcastStyle(s::FusedModBroadcast, t::FusedModBroadcast) = AbstractArrayStyle{1}
BroadcastStyle(s::FusedModBroadcast, t::AbstractArrayStyle{0}) = s
BroadcastStyle(s::AbstractArrayStyle{0}, t::FusedModBroadcast) = t
BroadcastStyle(s::FusedModBroadcast, t::BroadcastStyle) = t

similar(bc::Broadcasted{FusedModBroadcast{F}}, ::Type{F}) where F <: PrimeField = similar(Array{F}, axes(bc))

# transform / and // to a version that also works when reinterpreting as integers
broadcasted(st::FusedModBroadcast, ::Union{typeof(/), typeof(//)}, arg1, arg2) = broadcasted(st, *, arg1, (inv ∘ eltype(st)).(arg2))
FusableOps = Union{typeof(+), typeof(-), typeof(*), typeof(^)}
broadcasted(st::FusedModBroadcast, f::FusableOps, args...) = Broadcasted{typeof(st)}(f, args)
broadcasted(st::FusedModBroadcast, f::Function, args...) = broadcasted(DefaultArrayStyle{1}(), f, args...)

_boundedtype(T::Type{<:PrimeField}) = BoundedInteger{0:char(T)-1, inttype(T)}
_reinterpret(F, x) = x
_reinterpret(F, x::Integer) = _reinterpret(F, F(x))
_reinterpret(F, x::PrimeField) = reinterpret(_boundedtype(typeof(x)), x.n)
_reinterpret(F, x::AbstractArray{<:PrimeField}) = reinterpret(_boundedtype(eltype(x)), x)
function copyto!(dest::AbstractArray{F}, bc::Broadcasted{FusedModBroadcast{F}}) where F <: PrimeField
    bcf = flatten(bc)
    args = map(x -> _reinterpret(F, x), bcf.args)
    red(x) = posmod(x, char(F)) % inttype(F)
    copyto!(reinterpret(inttype(F), dest), broadcasted(red ∘ bcf.f, args...))
    dest
end

end # module
