"""
Reinterpretation between the bits types IntX and PrimeField. These have the same
underlying representation and we want to be able to performantly convert between
them.

This would be _almost_ done if we just implemented reinterpret. The only issue is
that `_getindex_ra` and `_setindex_ra!`` use `fieldcount(T) == 0` to determine whether
T is suitable for re-interpretation. In our cases, `fieldcount(T) == 1`, so
we have to override these two methods.
"""
module Reinterpret

import Base: reinterpret, ReinterpretArray, _getindex_ra, _setindex_ra!
import Base: @propagate_inbounds

import GaloisFields: PrimeField, Reduced

@inline reinterpret(F::Type{<:PrimeField{I}}, x::I) where I = F(Reduced(), x)
@inline reinterpret(T::Type{I}, x::PrimeField{I})   where I = x.n

const CONVERSIONS = [
    (:(PrimeField{I}), :I, :I),
    (:I, :(PrimeField{I}), :I),
]

for (From, To, Where) in CONVERSIONS
    @eval begin
        @inline @propagate_inbounds _getindex_ra(a::ReinterpretArray{<:$To, N, <:$From}, i1::Int, tailinds::TT) where {N, TT} where $Where = reinterpret(eltype(a), a.parent[i1, tailinds...])
        @inline @propagate_inbounds function _setindex_ra!(a::ReinterpretArray{<:$To, N, <:$From}, v, i1::Int, tailinds::TT) where {N, TT} where $Where
           v = convert(eltype(a), v)::eltype(a)
           return setindex!(a.parent, reinterpret(eltype(a.parent), v), i1, tailinds...)
         end
    end
end

end
