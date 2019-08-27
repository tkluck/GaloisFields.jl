using Random
using Primes

# Some properties of primitive roots we're using here:
#
# Definition 1: A value x is an n-root of unity in 𝔽 if x^n = 1.
# Definition 2: A value x ia a *primitive* n-th root of unity in 𝔽, if in
#               addition to being an n-th root of unity, we have x^k ≠ 1 for
#               any $k < 1$.
#
# Lemma 1: Primitive n-th roots of unity generate the group of n-th roots of unity.
#
# Proof:
# For any $i,j < n$, $x^i$ and $x^j$ must be distinct. If not, we must have
# $x^(i-j)=1$ contradicting the primitive-ness assumption. Then, since degree
# $n$ polynomial $x^n - 1$ can have at most $n$, solutions, the $x^i$ are all
# the n-th roots of unity.
#
# Lemma 2: If x is a primitive n-th roof of unity, then so is x^i for any i,
#          co-prime to n.
# Proof:
# Suppose not, then x^{ij} = 1 for some $j < n$. Then, since `x` is
# primitive, we must have $ij = an$ for some $a$, but since `i`, `n` are
# co-prime, $j$ must divide $n$ (which it can't since $j < n$). QED.
#
# Lemma 3: x is a primitive n-th root of unity if and only of x^(n/p_i) ≠ 1 for
#          all prime factors p_i of $n$.
#
# The forward definition follows by definition. For the reverse direction, let
# y be some primitive root. Then we must have $x = y^k$ for some $k$. Then,
# the reverse condition implies that $n$ and $k$ are co-prime (If $n$ and $k$
# share some prime factor $p_i$, then $n*k/p_i$ is an multiple of $n$, thus we
# would have x^(n/p_i) = 1). Thus by Lemma 2, x is a primitive n-th root of
# unity.

"""
    _rand_root(𝔽, n)

Obtain a random n-th root of unity (not necessarily primitive or ≢ 1)
"""
function _rand_root(rng::AbstractRNG, 𝔽::Type{<:PrimeField}, n)
    # Assumes gcd(char(𝔽) - 1, n) == n
    r = 𝔽(rand(rng, 1:(char(𝔽) - 1)))
    r^div(char(𝔽) - 1, n)
end

"""
    is_primitive_root(𝔽, x, n)

Determine whether `x` is a primitive `n`-th root of unity in 𝔽.
"""
is_primitive_root(𝔽::Type{<:PrimeField}, x, n; n_factors = factor(n)) =
    all(i->!isone(x^(n ÷ i)), keys(n_factors)) # Lemma 3

"""
    any_primitive_root([rng,] 𝔽, n = char(𝔽) - 1)

Obtain some primitive n-th root of unity (not necessarily minimal).

Note that for a primitive n-th root to exist `n` must divide `char(𝔽) - 1`.
Callers should guarantee this invariant.

If `n`, is not-specified it defaults to `n-1`. The resulting `p-1`-st primitive
root of unity is also called a "primitive root of 𝔽", "primitive element of 𝔽",
"primitive root mod p" or "generator (of the multiplicative group 𝔽⋆)".
"""
function any_primitive_root(rng::AbstractRNG, 𝔽::Type{<:PrimeField}, n = char(𝔽) - 1)
    q = char(𝔽) - 1
    @assert gcd(q, n) == n
    n_factors = factor(n)
    # Some implementation use the below algorithm to first find a primitive
    # generator of 𝔽 and then raise it to the n-th power. However, since
    # totient(n) is multiplicative and totient(x)/x <= 1 \forall(x), and since
    # n divides q, we know that totient(q)/q < totient(n)/n, so this method will
    # always be faster.
    while true
        root = _rand_root(rng, 𝔽, n)
        # This accepts with probability totient(n)/n
        is_primitive_root(𝔽, root, n; n_factors=n_factors) && return root
    end
end
any_primitive_root(𝔽::Type{<:PrimeField}, n) =
    any_primitive_root(Random.GLOBAL_RNG, 𝔽, n)

"""
    minimal_primitive_root(𝔽, n = char(𝔽) - 1)

Obtain the (unique) minimal primitive n-th root of unity of the field 𝔽.
Minimality is taken under the order of the canonical embedding into ℤ (i.e.
the one used by `reinterpret` and the display functions).

Note that for a primitive n-th root to exist `n` must divide `char(𝔽) - 1`.
Callers should guarantee this invariant.

If `n`, is not-specified it defaults to `n-1`. The resulting `p-1`-st primitive
root of unity is also called a "primitive root of 𝔽", "primitive element of 𝔽",
"primitive root mod p" or "generator (of the multiplicative group 𝔽⋆)".
"""
function minimal_primitive_root(𝔽::Type{<:PrimeField{T}}, n) where {T}
    root = any_primitive_root(𝔽, n)

    # Iterate over all of them to find the minimal one. (By Lemma 2, iterating
    # over x^i for i co-prime is primitive. By Lemma 1, this is all of them).
    𝔽(minimum(reinterpret(T, root^i) for i in filter(i->isone(gcd(i,n)), 1:n-1)))
end
