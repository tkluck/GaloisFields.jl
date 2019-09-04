import Base: @pure, power_by_squaring, literal_pow, widemul

"""
    F = ExtensionField{F <: AbstractGaloisField, N, α, MinPoly}

Algebraic extension of a finite field ``F`` of degree ``N``.
"""
struct ExtensionField{F <: AbstractGaloisField, N, α, MinPoly, Conway} <: AbstractExtensionField
    coeffs :: NTuple{N, F}
    ExtensionField{F, N, α, MinPoly, Conway}(coeffs::NTuple{N, F}) where
        {F <: AbstractGaloisField, N, α, MinPoly, Conway} = new(coeffs)
end

basefield(::Type{ExtensionField{F, N, α, MinPoly, Conway}}) where {F, N, α, MinPoly, Conway} = F
char(::Type{ExtensionField{F, N, α, MinPoly, Conway}})      where {F, N, α, MinPoly, Conway} = char(F)
n(::Type{ExtensionField{F, N, α, MinPoly, Conway}})         where {F, N, α, MinPoly, Conway} = N
genname(::Type{ExtensionField{F, N, α, MinPoly, Conway}})   where {F, N, α, MinPoly, Conway} = α
minpoly(::Type{ExtensionField{F, N, α, MinPoly, Conway}})   where {F, N, α, MinPoly, Conway} = MinPoly
isconway(::Type{ExtensionField{F, N, α, MinPoly, Conway}})  where {F, N, α, MinPoly, Conway} = Conway

expansion(a::ExtensionField) = a.coeffs

# -----------------------------------------------------------------------------
#
# Addition and substraction
#
# -----------------------------------------------------------------------------
zero(T::Type{<:ExtensionField}) =
    T(ntuple(i -> zero(basefield(T)), n(T)))
one( T::Type{<:ExtensionField}) =
    T(ntuple(i -> i == 1 ? one(basefield(T)) : zero(basefield(T)), n(T)))
gen( T::Type{<:ExtensionField}) =
    T(ntuple(i -> i == 2 ? one(basefield(T)) : zero(basefield(T)), n(T)))

+(a::F, b::F) where F<:ExtensionField = F(a.coeffs .+ b.coeffs)
-(a::F, b::F) where F<:ExtensionField = F(a.coeffs .- b.coeffs)

+(a::ExtensionField) = copy(a)
-(a::ExtensionField) = typeof(a)(.-a.coeffs)

iszero(a::ExtensionField) = all(iszero, a.coeffs)

# -----------------------------------------------------------------------------
#
# The interesting extension field operations: multiplication and division
#
# -----------------------------------------------------------------------------
function _rem(a::AbstractVector{C}, b::AbstractVector{C}) where C
    a = copy(a)
    len_a = findlast(!iszero, a)
    len_b = findlast(!iszero, b)
    while len_a !== nothing && len_a >= len_b
        q = a[len_a] // b[len_b]
        a[len_a-len_b+1:len_a] .-= q .* b[1:len_b]
        len_a = findprev(!iszero, a, len_a - 1)
    end
    a
end

function _gcdx(a::AbstractVector{C}, b::AbstractVector{C}) where C
    a = copy(a)
    b = copy(b)
    len_a = something(findlast(!iszero, a), 0)
    len_b = something(findlast(!iszero, b), 0)
    m = max(len_a, len_b)
    s0, s1 = zeros(C, m), zeros(C, m)
    t0, t1 = zeros(C, m), zeros(C, m)
    s0[1] = one(C)
    t1[1] = one(C)
    # The loop invariant is: s0*a0 + t0*b0 == a
    while len_b != 0
        s2 = copy(s0)
        t2 = copy(t0)
        while len_a >= len_b
            q = a[len_a] // b[len_b]
            a[len_a-len_b+1:len_a] .-= q .* b[1:len_b]

            deg_diff = len_a - len_b
            s2[deg_diff+1:end] .-= q .* s1[1:end-deg_diff]
            t2[deg_diff+1:end] .-= q .* t1[1:end-deg_diff]

            len_a = something(findprev(!iszero, a, len_a-1), 0)
        end
        a, b = b, a
        len_a, len_b = len_b, len_a
        s0, s1 = s1, s2
        t0, t1 = t1, t2
    end
    return (a, s0, t0)
end

@pure @generated function _pow_of_generator_rem(::Type{F}, m) where F <: ExtensionField
    N = n(F)
    pows = map(N : 2N - 2) do i
        pow_of_generator = zeros(basefield(F), 2N - 1)
        pow_of_generator[i + 1] = one(basefield(F))
        pow_of_generator_rem = _rem(pow_of_generator, collect(minpoly(F)))
        (i, F(tuple(pow_of_generator_rem[1:N]...)))
    end
    quote
        $(
        map(pows) do p
            :( m == $(p[1]) && return $(p[2]) )
        end...
        )
        error("_pow_of_generator_rem: m (=$m) should be between $(n(F)) and $(2n(F) - 1)")
    end
end

*(::Direct, a::F, b::F) where F <: ExtensionField = _mul_generated(F, a.coeffs, b.coeffs)

literal_pow(::typeof(^), a::F, ::Val{2}) where F <: ExtensionField = _mul_generated(F, a.coeffs)


function inv(::Direct, a::F) where F <: ExtensionField
    iszero(a) && throw(DivideError())
    N = n(F)
    coeffs = collect(a.coeffs)
    d, u, v = _gcdx(coeffs, collect(minpoly(F)))
    @assert !iszero(d[1]) && all(iszero, @view d[2:end])
    u ./= d[1]
    return F(ntuple(i -> u[i], n(F)))
end

/(::Direct, a::F, b::F) where F <: ExtensionField = a * inv(b)
//(::Direct, a::F, b::F) where F <: ExtensionField = a * inv(b)
^(::Direct, a::F, n::Integer) where F <: ExtensionField = power_by_squaring(a, n)

*(a::F, b::F)       where F <: ExtensionField = zech_op(F, *, a, b)
/(a::F, b::F)       where F <: ExtensionField = zech_op(F, /, a, b)
//(a::F, b::F)      where F <: ExtensionField = zech_op(F, //, a, b)
^(a::F, n::Integer) where F <: ExtensionField = zech_op(F, ^, a, n)
inv(a::F)           where F <: ExtensionField = zech_op(F, inv, a)

function power_by_squaring(a::F, n::Integer) where F <: ExtensionField
    iszero(a) && return iszero(n) ? one(a) : zero(a)
    # TODO: when length(F) is a BigInt, this allocates BigInt(n),
    # which is a wasteful allocation when n is small.
    n = mod(n, length(F) - 1)
    if n >= div(length(F) - 1, 2)
        n -= length(F) - 1
    end
    if n < 0
        a, n = inv(a), -n
    end
    n == 0 && return one(a)

    t = trailing_zeros(n) + 1
    n >>= t
    while (t -= 1) > 0
        a = a^2
    end
    b = a
    while n > 0
        t = trailing_zeros(n) + 1
        n >>= t
        while (t -= 1) >= 0
            a = a^2
        end
        b *= a
    end
    return b
end

@inline function _convolution(F, N, lo, hi, vec1, vec2)
    res = sum(vec1[i] * vec2[N - i + 1] for i in lo:hi)
    return F(res)
end

@inline _convolution(F, N, lo, hi, vec) = _convolution(F, N, lo, hi, vec, vec)

const _ApplicableInteger = Union{Int8, Int16, Int32}
@inline function _convolution(::Type{F}, N, lo, hi, vec1::Vec, vec2::Vec) where Vec <: NTuple{M, F} where F <: PrimeField{<:_ApplicableInteger} where M
    p = char(F)
    STEP = leading_zeros(zero(Int64)) - 2(leading_zeros(zero(p)) - leading_zeros(p))
    @assert STEP >= 1

    res = zero(Int64)
    for i in lo:STEP:hi
        for k in i:min(i + STEP - 1, hi)
            res += widemul(vec1[k].n, vec2[N - k + 1].n)
        end
        res = rem(res, p)
    end
    res = F(Reduced(), res)
    res
end

@inline function _convolution(::Type{F}, N, lo, hi, vec::Vec) where Vec <: NTuple{M, F} where F <: PrimeField{<:_ApplicableInteger} where M
    p = char(F)
    STEP = leading_zeros(zero(Int64)) - 2(leading_zeros(zero(p)) - leading_zeros(p) + 1)
    @assert STEP >= 1

    mid = div(lo + hi, 2)
    mid_counts_once = isodd(hi - lo + 1)

    res = zero(Int64)
    for i in lo:STEP:mid
        for k in i:min(i + STEP - 1, mid)
            factor = ifelse(mid_counts_once && k == mid, 1, 2)
            res += factor * widemul(vec[k].n, vec[N - k + 1].n)
        end
        res = rem(res, p)
    end
    res = F(Reduced(), res)
    res
end

# in case of length(operands) == 1: compute its square
# in case of length(operands) == 2: compute their product
# (ugly interface, but reduces code duplication)
@generated function _mul_generated(::Type{F}, operands...) where F <: ExtensionField
    N = n(F)
    coeffs = [Symbol(:coeff, i) for i in 1:N]
    code = quote
    end
    for j in 1:N
        push!(code.args, :( $(coeffs[j]) = _convolution(basefield(F), $j, 1, $j, operands...)))
    end
    for i in N + 1 : 2N - 1
        c = :( if (q = _convolution(basefield(F), $i, $(i + 1 - N), $N, operands...)) |> !iszero
        end )
        coeffs_to_add = _pow_of_generator_rem(F, i - 1).coeffs
        for j in 1 : N
            if !iszero(coeffs_to_add[j])
                push!(c.args[2].args, :( $(coeffs[j]) += q * $(coeffs_to_add[j]) ))
            end
        end
        push!(code.args, c)
    end

    res_tuple_expr = :( tuple($(coeffs...)) )
    push!(code.args, :( return F($res_tuple_expr) ))

    code
end

# -----------------------------------------------------------------------------
#
# Integer and basefield operations
#
# -----------------------------------------------------------------------------
function *(a::F, b::ExtensionField{F}) where F <: AbstractGaloisField
    typeof(b)(a .* b.coeffs)
end

*(a::ExtensionField{F}, b::F) where F <: AbstractGaloisField = b * a

*(a::Integer, b::ExtensionField) = typeof(b)(a .* b.coeffs)
*(a::ExtensionField, b::Integer) = b * a

# -----------------------------------------------------------------------------
#
# Constructors and promotion
#
# -----------------------------------------------------------------------------
promote_rule(F::Type{<:ExtensionField}, ::Type{<:Integer}) = F
function convert(F::Type{<:ExtensionField}, i::Integer)
    convert(F, convert(basefield(F), i))
end

(::Type{F})(n::F) where F<:ExtensionField = F(n.coeffs)

# -----------------------------------------------------------------------------
#
# Random number
#
# -----------------------------------------------------------------------------
function rand(rng::AbstractRNG, ::SamplerType{F}) where F <: ExtensionField
    F(ntuple(i -> rand(rng, basefield(F)), n(F)))
end
