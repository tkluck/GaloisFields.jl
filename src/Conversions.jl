struct IdentifiedVar{sym} end

identifications = Dict()

macro identify(expr)
    if expr.head == :call && expr.args[1] == :(=>)
        sym, target = expr.args[2:end]
        quote
            let tgt = $(esc(target))
                $identifications[$(QuoteNode(sym)), typeof(tgt)] = tgt
            end
        end
    else
        throw("Usage: @identify <symbol> => <target value>")
    end
end

function promote_rule(::Type{K}, ::Type{L}) where K <: PrimeField where L <: PrimeField
    if char(K) == char(L)
        return PrimeField{promote_type(eltype(K), eltype(L)), char(L)}
    else
        return Union{}
    end
end

function promote_rule(::Type{K}, ::Type{L}) where K <: ExtensionField where L <: PrimeField
    if char(K) == char(L)
        # TODO: should possibly extend the integer type, as well
        return K
    else
        return Union{}
    end
end

_upgrade(F::Type{<:ExtensionField}, i::AbstractGaloisField) = _upgrade(F, _upgrade(basefield(F), i))
_upgrade(F::Type{<:ExtensionField{K}}, i::K) where
    K <: AbstractGaloisField =
    F(ntuple(j -> j == 1 ? i : zero(K), n(F)))


function convert(K::Type{<:AbstractGaloisField}, x::PrimeField)
    if char(K) == char(typeof(x))
        return _upgrade(K, x)
    else
        throw(InexactError("Cannot convert $x to $K"))
    end
end

@generated function promote_rule(::Type{K}, ::Type{L}) where K <: ExtensionField where L <: ExtensionField
    if K == L
        return quote
            return K
        end
    elseif promote_rule(basefield(K), L) == basefield(K)
        return quote
            return K
        end
    elseif char(K) == char(L)
        sym = genname(L)
        if (sym, K) in keys(identifications)
            return quote
                return K
            end
        else
            return quote
                return Union{}
            end
        end
    else
        return quote
            return Union{}
        end
    end
end

@generated function convert(::Type{K}, x::L) where K <: ExtensionField where L <: ExtensionField
    if K == L
        return quote
            return x
        end
    elseif promote_rule(basefield(K), L) == basefield(K)
        return quote
            return $_upgrade(K, x)
        end
    elseif char(K) == char(L)
        sym = genname(L)
        if (sym, K) in keys(identifications)
            target = identifications[sym, K]
            powers = map(i -> target^(i-1), 1:n(L))
            return quote
                # cannot use closures in an @generated function body.
                # This is why we don't just have
                #     mapreduce((c, p) -> c * p, +, zip(expansion(x), $powers))
                res = zero(K)
                for (c, p) in zip(expansion(x), $powers)
                    res += c * p
                end
                return res
            end
        else
            throw("Cannot convert $x to $K; use @GaloisFields.identify $sym => <target> to define the conversion")
        end
    else
        throw("Cannot convert $x to $K")
    end
end

(E::Type{ExtensionField{F, N, α, MinPoly}})(n::Integer) where {F, N, α, MinPoly} = convert(E, n)
(E::Type{PrimeField{I,p}})(n::PrimeField) where {I, p} = convert(E, n)
(E::Type{ExtensionField{F, N, α, MinPoly}})(n::PrimeField) where {F, N, α, MinPoly} = convert(E, n)
(E::Type{PrimeField{I,p}})(n::ExtensionField) where {I, p} = convert(E, n)
(E::Type{ExtensionField{F, N, α, MinPoly}})(n::ExtensionField) where {F, N, α, MinPoly} = convert(E, n)
