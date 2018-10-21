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
        throw("Not implemented: @identify $expr")
    end
end

function promote_rule(F::Type{<:AbstractGaloisField}, G::Type{<:PrimeField})
    if char(F) == char(G)
        return F
    else
        return Union{}
    end
end

function convert(T::Type{<:AbstractGaloisField}, x::PrimeField)
    if char(T) == char(typeof(x))
        return convert(T, x.n)
    else
        throw(InexactError("Cannot convert $x to $T"))
    end
end

@generated function promote_rule(::Type{AF}, ::Type{EF}) where AF <: AbstractGaloisField where EF <: ExtensionField
    if char(AF) == char(EF)
        sym = genname(EF)
        if nothing !== (target = get(identifications, (sym, AF), nothing))
            return quote
                return $AF
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

@generated function convert(::Type{AF}, x::EF) where AF <: AbstractGaloisField where EF <: ExtensionField
    if char(AF) == char(EF)
        sym = genname(EF)
        if nothing !== (target = get(identifications, (sym, AF), nothing))
            powers = map(i -> target^(i-1), 1:n(EF))
            return quote
                # cannot use closures in an @generated function body.
                # This is why we don't just have
                #     mapreduce((c, p) -> c * p, +, zip(x.n, $powers))
                res = zero($AF)
                for (c, p) in zip(x.n, $powers)
                    res += c * p
                end
                return res
            end
        else
            throw(InexactError("Cannot convert $x to $T; use @GaloisFields.identify $sym => <target> to define the conversion"))
        end
    else
        throw(InexactError("Cannot convert $x to $T"))
    end
end
