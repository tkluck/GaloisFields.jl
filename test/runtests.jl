using Test
using Primes
using GaloisFields

@testset "GaloisFields" begin
    @testset "Prime field arithmetic" begin
        F = @GaloisField ℤ/37ℤ
        @test char(F) == 37
        @test repr(F) == "𝔽₃₇"

        @test F(1) + F(-1) == 0
        @test F(1) + F(1) == F(2)
        @test F(36) + F(1) == 0

        @test F(35) / F(35) == 1
        @test F(34) // F(34) == 1
        @test F(34) // F(16) * F(16) == 34

        @test zero(F) + one(F) == 1
        @test iszero(zero(F))
        @test iszero(char(F) * one(F))

        @test iszero(-F(0))

        # test for correct handling of integer overflow
        for I in [Int8, Int16, Int32, Int64, Int128]
            p = prevprime(typemax(I))
            G = GaloisField(p)
            @test G(-1) * G(-1) == 1
            @test G(-1) + G(-1) == -2
            @test G(0) - G(-1) == 1
        end
    end

    @testset "Integer promotions" begin
        F = @GaloisField ℤ/37ℤ

        @test F(2) + 4 == F(6)
        @test F(2) * 4 == F(8)
        @test F(2) / 4 == F(2) / F(4)
        @test F(2) // 4 == F(2) // F(4)
        @test 2 + F(4) == F(6)
        @test 2 * F(4) == F(8)
        @test 2 / F(4) == F(2) / F(4)
        @test 2 // F(4) == F(2) / F(4)
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

        @test_throws GaloisFields.InclusionError G(β)
    end

    @testset "Iterations" begin
        I = @GaloisField ℤ/2ℤ
        J = @GaloisField ℤ/3ℤ
        F = @GaloisField! 𝔽₂ α^2 + α + 1
        G = @GaloisField! 𝔽₅ α^2 - 2
        H = @GaloisField! G   β^3 + β + 1
        K = @GaloisField! 2^2 α
        L = @GaloisField! 5^2 α
        M = @GaloisField! 5^6 α
        for Q in [I, J, F, G, H, K, L, M]
            @test all(+x == x for x in Q)
            @test all(-x == 0 - x for x in Q)
            @test all(x^0 == 1 for x in Q)
            @test all(x^length(Q) == x for x in Q)
            @test all(x * inv(x) == 1 for x in Q if !iszero(x))
        end

        for Q in [I, J, F, G, K, L] # the smaller ones
            @test all( x + y == y + x for x in Q for y in Q)
            @test all( +x == x for x in Q for y in Q)
            @test all( x * y == y * x for x in Q for y in Q)
            @test all( x / y == x * inv(y) for x in Q for y in Q if !iszero(y))
        end
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

        @test x^3 + x + 1 == y^3 + y + 1
        @test (x^3 + x) / (x + 1) == (y^3 + y) / (y + 1)

        H = @GaloisField! 2^20 z
        K = @GaloisField! 2^20 w

        GaloisFields.enable_zech_multiplication(H)
        GaloisFields.disable_zech_multiplication(K)

        @test z^100 + z + 1 == w^100 + w + 1
        @test (z^100 + z) / (z + 1) == (w^100 + w) / (w + 1)

        @test_throws GaloisFields.InclusionError H(x)
        @test_throws GaloisFields.InclusionError F(z)
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
end
