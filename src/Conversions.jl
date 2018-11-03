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

function promote_rule(::Type{K}, ::Type{L}) where K <: PrimeField{I} where L <: PrimeField{I} where I
    if char(K) == char(L)
        return PrimeField{I, char(L)}
    else
        return Union{}
    end
end

function promote_rule(::Type{K}, ::Type{L}) where K <: AbstractExtensionField where L <: PrimeField
    if char(K) == char(L)
        # TODO: should possibly extend the integer type, as well
        return K
    else
        return Union{}
    end
end

_upgrade(F::Type{<:AbstractExtensionField}, i::AbstractGaloisField) = _upgrade(F, _upgrade(basefield(F), i))
_upgrade(F::Type{<:ExtensionField{K}}, i::K) where
    K <: AbstractGaloisField =
    F(ntuple(j -> j == 1 ? i : zero(K), n(F)))
_upgrade(F::Type{<:BinaryField}, i::PrimeField{I, 2}) where I =
    F(Bits(), i.n)

function convert(K::Type{<:AbstractGaloisField}, x::PrimeField)
    if char(K) == char(typeof(x))
        return _upgrade(K, x)
    else
        throw(InexactError("Cannot convert $x to $K"))
    end
end

@generated function promote_rule(::Type{K}, ::Type{L}) where K <: AbstractExtensionField where L <: AbstractExtensionField
    if K == L
        return quote
            return K
        end
    elseif promote_rule(basefield(K), L) == basefield(K)
        return quote
            return K
        end
    elseif char(K) == char(L)
        if isconway(K) && isconway(L)
            if n(K) % n(L) == 0
                return quote
                    return K
                end
            else
                return quote
                    return Union{}
                end
            end
        end
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

@generated function convert(::Type{K}, x::L) where K <: AbstractExtensionField where L <: AbstractExtensionField
    if K == L
        return quote
            return x
        end
    elseif promote_rule(basefield(K), L) == basefield(K)
        return quote
            return $_upgrade(K, x)
        end
    elseif char(K) == char(L)
        if isconway(K) && isconway(L)
            if n(K) % n(L) == 0
                p = char(K)
                m = (p^n(K) - 1) ÷ (p^n(L) - 1)
                target = gen(K)^m
            else
                throw("There is no inclusion from $L to $K")
            end
        else
            sym = genname(L)
            if (sym, K) in keys(identifications)
                target = identifications[sym, K]
            else
                throw("Cannot convert $x to $K; use @GaloisFields.identify $sym => <target> to define the conversion")
            end
        end
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
        throw("Cannot convert $x to $K")
    end
end

for to in [BinaryField, ExtensionField]
    @eval (E::Type{<:$to})(n::Integer) = convert(E, n)
end

for to in [PrimeField, BinaryField, ExtensionField]
    for from in [PrimeField, BinaryField, ExtensionField]
        @eval (E::Type{<:$to})(n::$from) = convert(E, n)
    end
end
