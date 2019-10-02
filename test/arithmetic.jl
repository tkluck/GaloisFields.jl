using Test
using LinearAlgebra: norm, tr

const MAXITERATIONS = 100
const MAXITERATIONS2 = round(Int, sqrt(MAXITERATIONS))
const MAXITERATIONS3 = round(Int, cbrt(MAXITERATIONS))

@testset "GaloisFields" begin
    @testset "Constructors" begin
        F = @GaloisField 3
        @test char(F) == 3
        F, α = @GaloisField 9
        @test char(F) == 3

        F = @GaloisField ℤ/3ℤ
        @test char(F) == 3

        F = @GaloisField ℤ/170141183460469231731687303715884105727ℤ
        @test char(F) == 170141183460469231731687303715884105727
        F = @GaloisField 𝔽₁₇₀₁₄₁₁₈₃₄₆₀₄₆₉₂₃₁₇₃₁₆₈₇₃₀₃₇₁₅₈₈₄₁₀₅₇₂₇
        @test char(F) == 170141183460469231731687303715884105727

        p = 29
        F = @GaloisField ℤ/(p*ℤ)
        @test char(F) == p

        n = 10
        F = @GaloisField! p^n δ
        @test char(F) == p
        @test δ^(Int64(p)^n) == δ
    end

    @testset "Arithmetic in $F" for F in TestFields
        @test startswith(repr(F), "𝔽")

        @test F(1) + F(-1) == 0
        @test F(1) + F(1) == F(2)
        @test F(char(F) - 1) + F(1) == 0

        iszero(F(35)) || @test F(35) / F(35) == 1
        iszero(F(34)) || @test F(34) // F(34) == 1
        iszero(F(16)) || @test F(34) // F(16) * F(16) == 34
        @test_throws DivideError F(1) // F(0)
        @test_throws DivideError F(0) // F(0)

        @test zero(F) + one(F) == 1
        @test iszero(zero(F))
        @test iszero(char(F) * one(F))

        @test iszero(-F(0))

        @test F(-1) * F(-1) == 1
        @test F(-1) + F(-1) == -2
        @test F(0) - F(-1) == 1
    end

    @testset "Integer promotion with $F" for F in TestFields
        @test F(2) + 4 == F(6)
        @test F(2) * 4 == F(8)
        @test 2 + F(4) == F(6)
        @test 2 * F(4) == F(8)
        if !iszero(F(41))
            @test F(2) / 41 == F(2) / F(41)
            @test F(2) // 41 == F(2) // F(41)
            @test 2 / F(41) == F(2) / F(41)
            @test 2 // F(41) == F(2) // F(41)
        end
        @testset "Overflow of $I" for I in [Int8, Int16, Int32, Int64, Int128]
            for i in (typemin(I), typemax(I))
                @test F(2) + i == F(big"2" + i)
                @test F(2) - i == F(big"2" - i)
                @test F(2) * i == F(big"2" * i)
            end
        end
    end

    @testset "Extensions of 𝔽₃" begin
        G = @GaloisField! 𝔽₃ α^2 + 1
        H = @GaloisField! 𝔽₃ β^2 + 1
        GaloisFields.identify(α => -β)
        @test char(G) == 3
        @test repr(G) == "𝔽₉"

        @test G(1) + G(-1) == 0

        @test α^2 + 1 == 0
        @test β^2 + 1 == 0

        @test (1 + α) // (1 + α) == 1
        @test (1 - α) // (1 + α) == 2α

        @test α + β == 0
        @test H(α) + β == 0

        # β + 1 doesn't satisfy minimum polynomial
        @test_throws GaloisFields.InclusionError GaloisFields.identify(α => β + 1)
    end

    @testset "Extensions of 𝔽₂" begin
        G = @GaloisField! 𝔽₂ α^2 + α + 1
        H = @GaloisField! 𝔽₂ β^2 + β + 1
        GaloisFields.identify(α => β + 1)
        @test char(G) == 2
        @test repr(G) == "𝔽₄"

        @test G(1) + G(-1) == 0

        @test α^2 + α + 1 == 0
        @test β^2 + β + 1 == 0

        @test (1 + α) // (1 + α) == 1
        @test (1 + α) // α == α

        @test α - β == 1
        @test H(α) - β == 1

        # test for correct handling of integer overflow
        for n in [7:9; 15:17; 31:33; 63:65]
            F, α = GaloisField(2, n)
            @test  α^(n - 1) // α^(n - 1) == 1
            @test  α^(n + 0) // α^(n + 0) == 1
            @test  α^(n + 1) // α^(n + 1) == 1
        end
    end

    @testset "Nested extension of 𝔽₂₉" begin
        G = @GaloisField! 𝔽₂₉ α^2 - 2
        H = @GaloisField! G   β^3 + 2β + 1
        K = @GaloisField! H   γ^7 - 2

        @test H(1) + H(-1) == 0

        @test H(α)^2 == 2
        @test K(α)^2 == 2
        @test β^3 + 2β + 1 == 0
        @test K(β)^3 + 2K(β) + 1 == 0
        @test γ^7 == 2

        @test α + β == β + α
        @test α + β + γ == γ + β + α

        @test G(H(α)) == α
        @test G(K(α)) == α

        @test norm(G, γ) isa G
        @test tr(G, γ) isa G
        @test norm(H, γ) isa H
        @test tr(H, γ) isa H

        @test_throws GaloisFields.InclusionError G(β)
    end

    @testset "Full axiom tests for $F" for F in TestFields
        if length(F) < MAXITERATIONS
            elements = F
        else
            elements = rand(F, MAXITERATIONS)
        end
        @test all(+x == x                       for x in elements)
        @test all(zero(x) + x == x              for x in elements)
        @test all(one(F) * x == x               for x in elements)
        @test all(x + -x == 0 == -x + x         for x in elements)
        @test all(x^0 == 1                      for x in elements)
        @test all(x^1 == x                      for x in elements)
        @test all(x^2 == x * x                  for x in elements)
        @test all(x^3 == x * x * x              for x in elements)
        @test all(x^4 == x * x * x * x          for x in elements)
        @test all(x^5 == x * x * x * x * x      for x in elements)
        @test all(x^length(F) == x              for x in elements)
        @test all(x * inv(x) == 1 == inv(x) * x for x in elements if !iszero(x))
        @test all(inv(x) == x^(-1)              for x in elements if !iszero(x))

        if length(F) < MAXITERATIONS2
            pairs = [(x, y) for x in F for y in F]
        else
            pairs = [(rand(F), rand(F)) for _ in 1:MAXITERATIONS]
        end
        @test all( x + y == y + x      for (x, y) in pairs)
        @test all( x * y == y * x      for (x, y) in pairs)
        @test all( x / y == x * inv(y) for (x, y) in pairs if !iszero(y))

        if length(F) < MAXITERATIONS3
            triples = [(x, y, z) for x in F for y in F for z in F]
        else
            triples = [(rand(F), rand(F), rand(F)) for _ in 1:MAXITERATIONS]
        end
        @test all((x + y) + z == x + (y + z)   for (x, y, z) in triples)
        @test all((x * y) * z == x * (y * z)   for (x, y, z) in triples)
        @test all((x + y) * z == x * z + y * z for (x, y, z) in triples)
        @test all(x * (y + z) == x * y + x * z for (x, y, z) in triples)
    end

    @testset "Conway polynomial database" begin
        K = @GaloisField! 29^4 α
        @test α^(29^4) == α
        L = @GaloisField! 29^2 β

        # Conway polynomials' compatibility conditions give a commutative
        # diagram of inclusions between them
        @test β == α^((29^4 - 1)÷(29^2 - 1))

        # same tests, but now employ Primes to factorize q
        # at construction time
        K,α = GaloisField(29^4)
        @test α^(29^4) == α
        L,β = GaloisField(29^2)
        @test β == α^((29^4 - 1)÷(29^2 - 1))

        M = @GaloisField! 81 γ
        N = @GaloisField! 9 δ
        @test γ^10 == δ

        @test (2γ)^10 == 2^10 * δ

        # Conway identification even with different variable names
        @GaloisField! 17^2 x
        @GaloisField! 17^2 y
        x^3 + x == y^3 + y

        # Rather big fields (make Int64 explicit for 32-bit platforms)
        @GaloisField! 2^50 z
        @test z^(Int64(2)^50) == z

        # Include two different fields in the smallest one that contains both
        @GaloisField! 2^4 x
        @GaloisField! 2^6 z
        @test x^((2^12 - 1)÷(2^6 - 1)) == z^((2^12 - 1)÷(2^4 - 1))

    end

    @testset "Zech logarithms" begin
        F = @GaloisField! 3^8 x
        G = @GaloisField! 3^8 y

        GaloisFields.enable_zech_multiplication(F)
        GaloisFields.disable_zech_multiplication(G)

        # double-check that x and y are identified because we use Conway
        # polynomials for both F and G....
        @test x == y
        # ...which allows us to convert F to G
        Fpairs = rand(F, 100, 2)
        Gpairs = map(G, Fpairs)

        # at which indices can we divide by F[:, 2] ?
        nz = findall(!iszero, Fpairs[:, 2])

        # the actual tests
        @test all(Fpairs[:,  1] .* Fpairs[:,  2] == Gpairs[:,  1] .* Gpairs[:,  2])
        @test all(Fpairs[nz, 1] ./ Fpairs[nz, 2] == Gpairs[nz, 1] ./ Gpairs[nz, 2])
        @test all(Fpairs[:,  1] .^ 2             == Gpairs[:,  1] .^ 2)

        H = @GaloisField! 2^20 z
        K = @GaloisField! 2^20 w

        GaloisFields.enable_zech_multiplication(H)
        GaloisFields.disable_zech_multiplication(K)

        @test z^100 + z + 1 == w^100 + w + 1
        @test (z^100 + z) / (z + 1) == (w^100 + w) / (w + 1)
    end

    @testset "Display" begin
        I = @GaloisField ℤ/2ℤ
        @test repr(I(0)) == "0"
        @test repr(I(1)) == "1"

        F = @GaloisField! 𝔽₂ α^2 + α + 1
        @test repr(F(0)) == "0"
        @test repr(α) == "α"
        @test repr(α + 1) == "α + 1"
        G = @GaloisField! 𝔽₅ α^2 - 2
        H = @GaloisField! G   β^3 + β + 1
        @test repr(G(0)) == "0"
        @test repr(H(0)) == "0"
        @test repr(α + β) == "β + α"
        @test repr(α * β) == "α * β"
        @test repr(α * β + β) == "(α + 1) * β"
        K = @GaloisField! 2^2 α
        @test repr(α^2) == "α + 1"
        L = @GaloisField! 5^2 α
        @test repr(α - 2) == "α + 3"
        M = @GaloisField! 5^6 α
        @test repr(3α^3 - 2) == "3 * α^3 + 3"
    end

    @testset "Broadcast" begin
        F = @GaloisField 𝔽₂₉

        x = rand(1:char(F), 100)
        y = rand(1:char(F)-1, 100)

        @test F[x;] .+ F[y;] == F.(x .+ y)
        @test F[x;] .* F[y;] == F.(x .* y)
        @test F[x;] .- F[y;] == F.(x .- y)
        @test F[x;] ./ F[y;] == F.(x .* invmod.(y, char(F)))

        @test tuple(F.(x)...) .+ tuple(F.(y)...) == tuple(F.(x .+ y)...)
        @test tuple(F.(x)...) .* tuple(F.(y)...) == tuple(F.(x .* y)...)
        @test tuple(F.(x)...) .- tuple(F.(y)...) == tuple(F.(x .- y)...)
        @test tuple(F.(x)...) ./ tuple(F.(y)...) == tuple(F.(x .* invmod.(y, char(F)))...)

        @test F(x[1]) .+ F[y;] == F.(x[1] .+ y)
        @test x[1]    .+ F[y;] == F.(x[1] .+ y)

        @test F(x[1]) ./ F[y;]   == F.(x[1] .* invmod.(y, char(F)))
        @test x[1]    ./ F[y;]   == F.(x[1] .* invmod.(y, char(F)))
        @test F[x;]   ./ F(y[1]) == F.(x    .* invmod(y[1], char(F)))
        @test F[x;]   ./ y[1]    == F.(x    .* invmod(y[1], char(F)))

        @test F(x[1]) .// F[y;]   == F.(x[1] .* invmod.(y, char(F)))
        @test x[1]    .// F[y;]   == F.(x[1] .* invmod.(y, char(F)))
        @test F[x;]   .// F(y[1]) == F.(x    .* invmod(y[1], char(F)))
        @test F[x;]   .// y[1]    == F.(x    .* invmod(y[1], char(F)))

        @test 3F[x;] == 3 .* F[x;] == F(3) * F[x;] == F(3) .* F[x;]
        @test Int8(3) * F[x;] == Int8(3) .* F[x;] == F(3) * F[x;] == F(3) .* F[x;]

        @test F.(x) == F[x;]
        @test convert.(F, x) == F[x;]

        # corner case: fuse operations with intermediate results bigger than integer type
        @test F[x;] .* F[x;] .* F[x;] .* F[x;] == map(x -> x^4, F[x;])
    end

    @testset "Random selection" begin
        F = @GaloisField 𝔽₂₉
        G = @GaloisField! 𝔽₅ α^2 - 2
        H = @GaloisField! G   β^3 + β + 1
        K = @GaloisField! 𝔽₆₄ γ

        x = rand(F, 100)
        y = rand(G, 100)
        z = rand(H, 100)
        w = rand(K, 100)

        @test x .+ x .* x == map(a -> a + a * a, x)
        @test y .+ y .* y == map(a -> a + a * a, y)
        @test z .+ z .* z == map(a -> a + a * a, z)
        @test w .+ w .* w == map(a -> a + a * a, w)
    end

    @testset "Primitive roots of unity" begin
        let 𝔽₁₀₃₁ = GaloisField(1031), n = 103
            # The really naive way to check for primitive roots of unity
            # 1031 is small, so this is fast enough.
            naive_roots_of_unity = filter(1:1030) do x
                e = 𝔽₁₀₃₁(x)
                # Is this a root of unity?
                isone(e^n) || return false
                # Is it primitive?
                for i = 1:n-1
                    isone(e^i) && return false
                end
                return true
            end

            # Generate a bunch of primitive roots of unity and do some basic
            # sanity checks.
            let random_roots_of_unity = [GaloisFields.any_primitive_root(𝔽₁₀₃₁, n) for _ = 1:1000]
                @test all(x->x in naive_roots_of_unity, random_roots_of_unity)
                # Make sure they're not all the same
                @test any(x->x != random_roots_of_unity[1], random_roots_of_unity)
            end

            # Check to make sure that we're getting the correct minimum root of
            # unity.
            @test minimum(naive_roots_of_unity) ==
                GaloisFields.minimal_primitive_root(𝔽₁₀₃₁, n)
        end
    end

    @test_throws ErrorException GaloisField(10)
    @test_throws ErrorException GaloisField(10, 1)
end
