"""
    Direct()

A helper singleton for arithmetic operations that benefit from a log-table.
Field implementations implement e.g.

    *(::Direct, a, b)

and then define

    *(a::F, b::F) = zech_op(F, *, a, b)

This will use the `::Direct` implementation for computing a log table
where appropriate, and will either use the log table or the `::Direct`
implementation for computing results.
"""
struct Direct end

logtables = IdDict()
function cycliclogtable(cyclic_generator, limit)
    global logtables
    # get!(f, ...) doesn't have an implementation for IdDict
    if cyclic_generator in keys(logtables)
        return logtables[cyclic_generator]
    end
    table = [(one(cyclic_generator), 0)]
    x, n = cyclic_generator, 1
    while x != one(x)
        push!(table, (x, n))
        x = *(Direct(), x, cyclic_generator)
        n += 1
        @assert n <= limit
    end

    result = (Dict(x => n for (x, n) in table), Dict(n => x for (x, n) in table), n)
    logtables[cyclic_generator] = result
    result
end

function zech_op(logtable, exptable, n, ::typeof(*), a, b)
    iszero(a) && return zero(a)
    iszero(b) && return zero(b)
    exptable[mod(logtable[a] + logtable[b], n)]
end

function zech_op(logtable, exptable, n, ::Union{typeof(/), typeof(//)}, a, b)
    iszero(b) && throw(DivideError())
    exptable[mod(logtable[a] - logtable[b], n)]
end

function zech_op(logtable, exptable, n, ::typeof(inv), a)
    iszero(a) && throw(DivideError())
    exptable[mod(-logtable[a], n)]
end

function zech_op(logtable, exptable, n, ::typeof(^), a, N)
    iszero(a) && throw(DivideError())
    exptable[mod(N * logtable[a], n)]
end

overrides = Dict()

function enable_zech_multiplication(F::Type{<:AbstractGaloisField})
    overrides[F] = true
end

function disable_zech_multiplication(F::Type{<:AbstractGaloisField})
    overrides[F] = false
end

@generated function zech_op(::Type{F}, op, args...) where F <: AbstractGaloisField
    p = char(F)
    q = p^n(F)
    default_condition = p != 2 && q < 2^16
    if get(overrides, F, default_condition)
        logtable, exptable, order = cycliclogtable(gen(F), q - 1)
        if order == q - 1
            return quote
                zech_op($logtable, $exptable, $order, op, args...)
            end
        else
            if F in keys(overrides)
                # the user explicitly requested Zech logarithms, but we
                # can't provide them. Throw a warning.
                @warn "Primitive element for $F is not a multiplicative generator. We cannot use Zech logarithms and therefore operations may be slow"
           end
            # fall through
        end
    end
    return quote
        op(Direct(), args...)
    end
end
