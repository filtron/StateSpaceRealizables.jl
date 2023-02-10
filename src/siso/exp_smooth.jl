
"""
    ExpSmoother{E,T<:Real} <: StateSpaceRealizable{E}

Type for representing exponential smoothers, i.e. the following frequency response

H(z) = ν! / √(2ν!) (2λ)^{ν + 1/2} / (z + λ)^{ν + 1},
"""
struct ExpSmoother{E,T<:Real} <: StateSpaceRealizable{E}
    λ::T
    ν::Integer
end

function ExpSmoother{E}(λ::T, ν::Integer) where {E<:TimeEvolution,T<:Real}
    return ExpSmoother{E,T}(λ, ν)
end

ninputs(h::ExpSmoother) = 1
nstates(h::ExpSmoother) = h.ν + 1
noutputs(h::ExpSmoother) = 1
isproper(h::ExpSmoother) = true
poles(h::ExpSmoother)    = fill(-h.λ, nstates(h))

# this is internally balanced:
function ssparams(h::ExpSmoother{<:ContinuousTE,T}) where {T}
    A = h.λ * (I - T(2) * tril(ones(nstates(h), nstates(h))))
    B = sqrt(T(2) * h.λ) * ones(nstates(h), 1)
    C =
        exp.(
            -logfactorial.(collect(0:nstates(h)-1)) -
            logfactorial.(collect(nstates(h)-1:-1:0)) .+
            2 * logfactorial(nstates(h) - 1) .-
            logfactorial(2 * nstates(h) - 2) / 2
        ) .* (-1) .^ collect(0:nstates(h)-1)
    return A, B, permutedims(C)
end

function sqr_magnitude_response!(mag2::AbstractVector{<:AbstractMatrix}, h::ExpSmoother{<:ContinuousTE,T}, zs) where {T}
    for (i, z) in enumerate(zs)
        mag2[i][1,1] = exp( T(2) * logfactorial(h.ν) - logfactorial(2 * h.ν)  + (T(2) * h.ν + one(T)) * log(2 * h.λ)) / (abs2(z) +  h.λ^2)^(h.ν + 1)
    end
end


function impulse_response!(hs::AbstractVector{<:AbstractMatrix}, h::ExpSmoother{<:ContinuousTE,T}, ts) where {T}
    for (i, t) in enumerate(ts)
        hs[i][1, 1] = exp( (h.ν + T(0.5)) * log(2 * h.λ) - T(0.5) * logfactorial(2 * h.ν) ) * t^h.ν * exp( -h.λ *  t)
    end
end