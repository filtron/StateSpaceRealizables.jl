function test_state_space()

    @testset "StateSpace" begin

        b = 5.0
        k = 101
        imaxis = collect(im * LinRange(-b, b, k))
        ts = collect(LinRange(0.0, b, k))

        A1 = -tril(ones(2, 2))
        B1 = ones(2, 2)
        C1 = [0.0 1.0]
        D1 = ones(1, 2)

        A2 = diagm(-1 => ones(2)) - diagm(0 => ones(3))
        B2 = vcat(1.0, zeros(2, 1))
        C2 = [1.0 -1.0 0.0; 1.0 -1.0 0.0]
        D2 = zeros(2, 1)

        A3 = randn(2, 2)
        B3 = randn(2, 2)
        C3 = randn(1, 2)
        D3 = randn(1, 2)

        for E in (ContinuousTE, DiscreteTE)
            S1 = StateSpace{E}(A1, B1, C1, D1)
            negS1 = -S1

            S2 = ProperStateSpace{E}(A2, B2, C2)
            negS2 = -S2

            S12 = S1 ∘ S2

            S3 = StateSpace{E}(A3, B3, C3, D3)

            S13 = S1 + S3

            S4 = ProperStateSpace{E}(A3, B3, C3)


            @test ninputs(S1) == 2
            @test noutputs(S1) == 1
            @test nstates(S1) == 2
            @test all(ssparams(S1) .== (A1, B1, C1, D1))
            @test all(ssparams(S2) .== (A2, B2, C2, nothing))
            @test ssrealize(S1) == S1
            @test ssrealize(S2) == S2
            @test ssrealize(A1, B1, C1, D1; te = E()) == S1
            @test isproper(S1) == IsNotProper()
            @test isproper(S2) == IsProper()

            @test_nowarn poles(S1)
            @test_nowarn S3 + S4
            @test_nowarn S4 + S3
            @test_nowarn S3 - S4
            @test_nowarn S4 - S3

            @test negS1.A == A1
            @test negS1.B == B1
            @test negS1.C == -C1
            @test negS1.D == -D1

            @test negS2.A == A2
            @test negS2.B == B2
            @test negS2.C == -C2

            @test ninputs(S12) == 1
            @test noutputs(S12) == 1
            @test nstates(S12) == nstates(S1) + nstates(S2)

            @test ninputs(S13) == ninputs(S1) == ninputs(S3)
            @test noutputs(S13) == noutputs(S1) == noutputs(S3)
            @test nstates(S13) == nstates(S1) + nstates(S3)

            E == ContinuousTE && @test_nowarn sqr_magnitude_response(S1, imaxis)
            E == ContinuousTE && @test_nowarn sqr_magnitude_response(S2, imaxis)

            E == ContinuousTE && @test_nowarn impulse_response(S1, ts)
            E == ContinuousTE && @test_nowarn impulse_response(S2, ts)
        end
    end






end
