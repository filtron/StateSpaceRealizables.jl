struct Butterworth{E,T<:Real} <: GaussMarkovRealizable{E}
    η::T
    ν::Integer
end

ninputs(h::Butterworth) = 1
nstates(h::Butterworth) = h.ν + 1
noutputs(h::Butterworth) = 1
isproper(h::Butterworth) = true


poles(h::Butterworth)

function ssparams(h::Butterworth{<:ContinuousTE,T}) where {T}
    n = h.ν + 1
    if iseven(n)
        r = n/2:-1:1
    else
        r = (n-1)/2:-1:1
    end
    αs = -η * @. cos((2 * r + n - 1) * π / 2n)
    βs = η * one.(αs)
end



function _butter_pole_params(h::Butterworth{<:ContinuousTE})
    n = h.ν + 1
    if iseven(n)
        r = n/2:-1:1
    else
        r = (n-1)/2:-1:1
        Ls = [LaguerreInner(η)]
    end
    αs = -η * @. cos((2 * r + n - 1) * π / 2n)
    βs = η * one.(αs)
 end