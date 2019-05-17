module Broadcast

import Base: copyto!, similar
import Base.Broadcast: Broadcasted, BroadcastStyle, broadcasted
import Base.Broadcast: DefaultArrayStyle, AbstractArrayStyle
import Base.Broadcast: flatten, copyto!

import SIMD: VecRange, vload, vstore, Vec, vifelse

import GaloisFields: inttype, char, PrimeField

function vecmod(x, y)
    T = typeof(x)
    z = rem(x, y)
    z += vifelse(z < 0, T(y), T(0))
    z
end

_DenseVector{T} = Union{
    DenseVector{T},
    # TODO: should check for unitrange having step size +1
    SubArray{T,1,Vector{T},Tuple{UnitRange{Int64}},true},
}
_SIMDable = PrimeField{<:Union{Int8, Int16, Int32, Int64}}

struct SIMDBroadcast{F} <: BroadcastStyle end
eltype(::SIMDBroadcast{F}) where F = F

BroadcastStyle(::Type{<:_DenseVector{F}}) where F <: _SIMDable = SIMDBroadcast{F}()
BroadcastStyle(::T, ::T) where T <: SIMDBroadcast = T()
BroadcastStyle(s::SIMDBroadcast, t::SIMDBroadcast) = AbstractArrayStyle{1}
BroadcastStyle(s::SIMDBroadcast, t::AbstractArrayStyle{0}) = s
BroadcastStyle(s::AbstractArrayStyle{0}, t::SIMDBroadcast) = t
BroadcastStyle(s::SIMDBroadcast, t::BroadcastStyle) = t

similar(bc::Broadcasted{SIMDBroadcast{F}}, ::Type{F}) where F <: PrimeField = similar(Array{F}, axes(bc))

_getslice(x, ix) = x
_getslice(x::PrimeField, ix) = x.n
_getslice(x::_DenseVector{F}, vr::VecRange{N}) where F <: PrimeField where N = vload(Vec{N, inttype(F)}, pointer(reinterpret(inttype(F), x), vr.i))
_setslice!(x::_DenseVector{F}, val::Vec{N, I}, vr::VecRange{N}) where F <: PrimeField{I} where {N, I} = vstore(val, pointer(reinterpret(inttype(F), x), vr.i))

# can't do invmod(...) on SIMD, so transform here
# TODO: there's many more operations we should similarly exclude, e.g. F.(x)
broadcasted(st::SIMDBroadcast, ::Union{typeof(/), typeof(//)}, arg1, arg2) = broadcasted(st, *, arg1, broadcasted(st, inv, arg2))
broadcasted(st::SIMDBroadcast, ::typeof(inv), arg::AbstractVector) = map(inv ∘ eltype(st), arg)
broadcasted(st::SIMDBroadcast, ::typeof(inv), arg::Number) = inv(eltype(st)(arg))

function copyto!(dest::_DenseVector{F}, bc::Broadcasted{SIMDBroadcast{F}}) where F <: PrimeField
    bcf = flatten(bc)
    f = bcf.f
    args = bcf.args

    m = length(dest)
    N = 8
    lane = VecRange{N}(1)
    for i = 0 : div(m, N) - 1
        j = i * N
        _setslice!(
            dest,
            # We're fusing the operation before doing the mod operation.
            # This is way faster, but it may overflow. To reduce the odds, we construct
            # PrimeField{I, p} always with a bigger I than necessary. See PrimeField.jl.
            # (TODO: A safer solution is to be implemented.)
            vecmod(f(map(a -> _getslice(a, lane + j), args)...), char(F)),
            lane + j,
        )
    end
    r = rem(m, N)
    for j = m - r + 1 : m
        dest[j] = bcf[j]
    end

    dest
end

end # module
