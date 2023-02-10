function test_expsmooth()



    @testset "Expsmoother" begin
        νs = 0:2
        λs = [0.5, 1.5, 2.5]

        hs = map(Base.splat(ExpSmoother{ContinuousTE}), zip(λs, νs))

        k = 11
        b = 5.0
        imaxis = im * collect(LinRange(-b, b, k))
        ts = collect(LinRange(0.0, b, k))

        for (i, h) in enumerate(hs)

            A, B, C = ssparams(h)

            @test ninputs(h) == 1
            @test noutputs(h) == 1
            @test nstates(h) == νs[i] + 1
            @test isproper(h) == true
            @test all( x-> x == -λs[i], poles(h))
            @test length(poles(h)) == nstates(h)

            # this assumes balanced realization ...
            @test A + A' ≈ - B*B'
            @test LinearAlgebra.norm_sqr(C) ≈ 1.0
            @test ssrealize(h) == ssrealize(A, B, C; te = ContinuousTE())
            @test_nowarn sqr_magnitude_response(h, imaxis)
            @test_nowarn impulse_response(h, ts)

        end

    end

end