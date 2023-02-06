function test_inner_functions()

    abstol = 1e-12

    k = 101
    b = 5.0
    imaxis = im * collect(LinRange(-b, b, k))

    λs = [1.0, 2.0, 3.0]
    Ls = LaguerreInner.(λs)

    αs = [2.0, 4.0, 6.0]
    βs = αs .^ 2 + 2π * [0.25, 0.5, 0.75]
    Ks = map(x -> KautzInner(x...), zip(αs, βs))

    P1 = CompositeInner(Ls, Ks)
    P2 = CompositeInner(Ls)
    P3 = CompositeInner(Ks)

    Ps = (P1, P2, P3)

    @testset "InnerFunction" begin

        @testset "LaguerreInner" begin

            for L in Ls
                ninputs(L) == 1
                noutputs(L) == 1
                isproper(L) == false
                S = ssrealize(L)
                (; A, B, C, D) = S
                A2, B2, C2, D2 = ssparams(L)
                @test all((A, B, C, D) .== (A2, B2, C2, D2))
                @test A + A' + B * B' ≈ zero(A) atol = abstol
                @test B ≈ -C' atol = abstol
                @test D * D' ≈ diagm(one.(diag(D)))
                @test all(abs2.(L.(imaxis)) .≈ 1.0)
            end

        end

        @testset "KautzInner" begin

            for K in Ks
                ninputs(K) == 1
                noutputs(K) == 1
                isproper(K) == false
                S = ssrealize(K)
                (; A, B, C, D) = S
                A2, B2, C2, D2 = ssparams(K)
                @test all((A, B, C, D) .== (A2, B2, C2, D2))
                @test A + A' + B * B' ≈ zero(A) atol = abstol
                @test B ≈ -C' atol = abstol
                @test D * D' ≈ diagm(one.(diag(D)))
                @test all(abs2.(K.(imaxis)) .≈ 1.0)
            end

        end

        @testset "CompositeInner" begin

            for P in Ps
                ninputs(P) == 1
                noutputs(P) == 1
                isproper(P) == false
                S = ssrealize(P)
                (; A, B, C, D) = S
                A2, B2, C2, D2 = ssparams(P)
                @test all((A, B, C, D) .== (A2, B2, C2, D2))
                @test A + A' + B * B' ≈ zero(A) atol = abstol
                @test B ≈ -C' atol = abstol
                @test D * D' ≈ diagm(one.(diag(D)))
                @test all(abs2.(P.(imaxis)) .≈ 1.0)
            end

            @test all(ssparams(P1) .≈ ssparams(reduce(*, Ls) * reduce(*, Ks)))
            @test all(ssparams(P2) .≈ ssparams(reduce(*, Ls)))
            @test all(ssparams(P3) .≈ ssparams(reduce(*, Ks)))

        end

    end

end