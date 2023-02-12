struct Oscillator{E,T<:Real} <: StateSpaceRealizable{E}
    α::T 
    β::T 
end

Oscillator{E}(α::T, β::T) where {E<:TimeEvolution,T<:Real} = 
    Oscillator{E,T}(α, β)

ninputs(h::Oscillator) = 1
nstates(h::Oscillator) = 2
noutputs(h::Oscillator) = 1
isproper(h::Oscillator) = true

function poles(h::Oscillator)
    if h.α > h.β
        φ = one(complex(T)) * sqrt((h.α + h.β) * (h.α - h.β))
    elseif K.α < K.β
        φ = im * sqrt((h.β + h.α) * (h.β - h.α))
    elseif K.α == K.β
        φ = one(complex(T)) * zero(h.α)
    end
    return [-h.α - φ, -h.α + φ]
end

function ssparams(h::Oscillator{E,T}) where {E,T}
    A = [-T(2)*h.α -h.β; h.β zero(T)] 
    B = [sqrt(T(4) * h.α); zero(T);;]
    C = [one(T) zero(T)] 
    return A,B, C    
end 

