using Test
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
    end

    @testset "Extensions of 𝔽₃" begin
        G = @GaloisField! 𝔽₃ α^2 + 1
        H = @GaloisField! 𝔽₃ β^2 + 1
        @GaloisFields.identify α => -β
        @test char(G) == 3
        @test repr(G) == "𝔽₉"

        @test α^2 + 1 == 0
        @test β^2 + 1 == 0

        @test (1 + α) // (1 + α) == 1
        @test (1 - α) // (1 + α) == 2α

        @test α + β == 0
        @test H(α) + β == 0
    end

    @testset "Extensions of 𝔽₂" begin
        G = @GaloisField! 𝔽₂ α^2 + α + 1
        H = @GaloisField! 𝔽₂ β^2 + β + 1
        @GaloisFields.identify α => β + 1
        @test char(G) == 2
        @test repr(G) == "𝔽₄"

        @test α^2 + α + 1 == 0
        @test β^2 + β + 1 == 0

        @test (1 + α) // (1 + α) == 1
        @test (1 + α) // α == α

        @test α - β == 1
        @test H(α) - β == 1
    end
end
