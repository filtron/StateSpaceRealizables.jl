"""
    AbstractInnerFunction{E<:TimeEvolution}

Abstract type for representing inner functions.
"""
abstract type AbstractInnerFunction{E} <: StateSpaceRealizable{E} end


isproper(::AbstractInnerFunction) = false

"""
    AbstractInnerFunction{E<:TimeEvolution}

Abstract type for representing inner functions.
"""
abstract type UnivariateInnerFunction{E,T<:Number} <: AbstractInnerFunction{E} end


"""
    LaguerreInner{E,T} <: UnivariateInnerFunction{E,T}

Type for representing Laguerre inner functions.
"""
struct LaguerreInner{E,T} <: UnivariateInnerFunction{E,T}
    λ::T
end

LaguerreInner{E}(λ::Real) where {E<:ContinuousTE} =
    λ > zero(λ) ? LaguerreInner{E,typeof(λ)}(λ) :
    throw(DomainError(λ, "λ must be positive."))

LaguerreInner(λ; te::TimeEvolution = ContinuousTE()) = LaguerreInner{typeof(te)}(λ)

(L::LaguerreInner)(z::Number) = (z - L.λ) / (z + L.λ)

"""
    LaguerreInner{E,T} <: UnivariateInnerFunction{E,T}

Type for representing Kautz inner functions.
"""
struct KautzInner{E,T} <: UnivariateInnerFunction{E,T}
    α::T
    β::T
end

function KautzInner{E}(α::T, β::T) where {E<:ContinuousTE,T<:Real}
    α < zero(α) && throw(DomainError(α, "α must be positive."))
    β < zero(β) && throw(DomainError(β, "β must be positive."))
    return KautzInner{E,T}(α, β)
end

KautzInner(α, β; te::TimeEvolution = ContinuousTE()) = KautzInner{typeof(te)}(α, β)

(K::KautzInner)(z::Number) = (z^2 - 2 * K.α * z + K.β^2) / (z^2 + 2 * K.α * z + K.β^2)


"""
    CompositeInner

Type for representing compositions of Laguerre and Kautz inner functions.
"""
struct CompositeInner{E,T,L,K} <: UnivariateInnerFunction{E,T}
    laguerre_inners::L
    kautz_inners::K
end

CompositeInner(
    Ls::A,
    Ks::B,
) where {E,T,A<:AbstractVector{<:LaguerreInner{E,T}},B<:AbstractVector{<:KautzInner{E,T}}} =
    CompositeInner{E,T,A,B}(Ls, Ks)

const CompositeLaguerreInner{E,T,A} = CompositeInner{E,T,A,<:Nothing}
const CompositeKautzInner{E,T,B} = CompositeInner{E,T,<:Nothing,B}

CompositeInner(Ls::AbstractVector{<:LaguerreInner{E,T}}) where {E,T} =
    CompositeInner{E,T,typeof(Ls),Nothing}(Ls, nothing)
# add explicit nothing constructors

CompositeInner(Ks::AbstractVector{KautzInner{E,T}}) where {E,T} =
    CompositeInner{E,T,Nothing,typeof(Ks)}(nothing, Ks)

CompositeInner(L::LaguerreInner{E,T}, K::KautzInner{E,T}) where {E,T} =
    CompositeInner([L], [K])

CompositeInner(L::LaguerreInner{T}) where {T} = CompositeInner([L])
CompositeInner(K::KautzInner{T}) where {T} = CompositeInner([K])

CompositeInner(C::CompositeInner) = C

function (C::CompositeLaguerreInner{<:ContinuousTE})(z::Complex)
    out = one(typeof(z))
    @simd ivdep for L in C.laguerre_inners
        out = out * L(z)
    end
    return out
end

function (C::CompositeKautzInner{<:ContinuousTE})(z::Complex)
    out = one(typeof(z))
    @simd ivdep for K in C.kautz_inners
        out = out * K(z)
    end
    return out
end

function (C::CompositeInner)(z::Complex)
    out = one(typeof(z))
    @simd ivdep for L in C.laguerre_inners
        out = out * L(z)
    end
    @simd ivdep for K in C.kautz_inners
        out = out * K(z)
    end
    return out
end

"""
    butter_inner(η::Number, ν::Integer)

Computes a composite inner function, of time evolution ContinuousTE, with poles matching those of
the Butterworth polynomial of order ν + 1.
"""
function butter_inner(η::Real, ν::Integer)
    n = ν + 1
    if iseven(n)
        r = n/2:-1:1
        Ls = nothing
    else
        r = (n-1)/2:-1:1
        Ls = [LaguerreInner(η)]
    end
    αs = -η * @. cos((2 * r + n - 1) * π / 2n)
    βs = η * one.(αs)
    Ks = KautzInner.(αs, βs)
    return CompositeInner{ContinuousTE,typeof(η),typeof(Ls),typeof(Ks)}(Ls, Ks)
end

ninputs(::UnivariateInnerFunction) = 1
noutputs(::UnivariateInnerFunction) = 1

nstates(L::LaguerreInner) = 1
nstates(K::KautzInner) = 2
nstates(C::CompositeLaguerreInner) = sum(nstates.(C.laguerre_inners))
nstates(C::CompositeKautzInner) = sum(nstates.(C.kautz_inners))
nstates(C::CompositeInner) =
    sum(nstates.(C.laguerre_inners)) + sum(nstates.(C.kautz_inners))

"""
    ssparams(H::AbstractInnerFunction)

Computes the state space parameters A, B, C, D for the
balanced realization of the inner function H.
"""
function ssparams(L::LaguerreInner{<:ContinuousTE,T}) where {T}
    A = -hcat(L.λ)
    B = hcat(sqrt(T(2) * L.λ))
    C = -B
    D = hcat(one(T))
    return A, B, C, D
end

function ssparams(K::KautzInner{<:ContinuousTE,T}) where {T}
    A = [-T(2)*K.α -K.β; K.β zero(T)]
    B = hcat([sqrt(T(4) * K.α); zero(T)])
    C = -B'
    D = hcat(one(T))
    return A, B, C, D
end

ssparams(H::CompositeInner) = ssparams(ssrealize(H))


"""
    ssrealize(H::AbstractInnerFunction)

Computes the balanced state-space realization of the inner function H.
"""
ssrealize(H::AbstractInnerFunction) = StateSpace{time_evolution(H)}(ssparams(H)...)
ssrealize(C::CompositeLaguerreInner) = reduce(∘, ssrealize.(C.laguerre_inners)) # can be made more efficient by pre-allocation and for loop?
ssrealize(C::CompositeKautzInner) = reduce(∘, ssrealize.(C.kautz_inners)) # can be made more efficient by pre-allocation and for loop?
ssrealize(C::CompositeInner) =
    reduce(∘, ssrealize.(C.kautz_inners)) ∘ reduce(∘, ssrealize.(C.laguerre_inners)) # can be made more efficient by pre-allocation and for loop?


∘(H2::UnivariateInnerFunction{E,T}, H1::UnivariateInnerFunction{E,T}) where {E,T} =
    ∘(CompositeInner(H2), CompositeInner(H1)) # fallback
*(H2::UnivariateInnerFunction{E,T}, H1::UnivariateInnerFunction{E,T}) where {E,T} = H2 ∘ H1

∘(H2::CompositeInner{E,T}, H1::CompositeInner{E,T}) where {E,T} = CompositeInner(
    vcat(H2.laguerre_inners, H1.laguerre_inners),
    vcat(H2.kautz_inners, H1.kautz_inners),
)

∘(H2::CompositeInner{E,T}, H1::CompositeLaguerreInner{E,T}) where {E,T} =
    CompositeInner(vcat(H2.laguerre_inners, H1.laguerre_inners), H2.kautz_inners)
∘(H2::CompositeLaguerreInner{E,T}, H1::CompositeInner{E,T}) where {E,T} = H1 ∘ H2

∘(H2::CompositeInner{E,T}, H1::CompositeKautzInner{E,T}) where {E,T} =
    CompositeInner(H2.laguerre_inners, vcat(H2.kautz_inners, H1.kautz_inners))
∘(H2::CompositeKautzInner{E,T}, H1::CompositeInner{E,T}) where {E,T} = H1 ∘ H2


∘(H2::CompositeLaguerreInner{E,T}, H1::CompositeLaguerreInner{E,T}) where {E,T} =
    CompositeInner(vcat(H2.laguerre_inners, H1.laguerre_inners))
∘(H2::CompositeLaguerreInner{E,T}, H1::CompositeKautzInner{E,T}) where {E,T} =
    CompositeInner(H2.laguerre_inners, H1.kautz_inners)
∘(H2::CompositeKautzInner{E,T}, H1::CompositeLaguerreInner{E,T}) where {E,T} = H1 ∘ H2
∘(H2::CompositeKautzInner{E,T}, H1::CompositeKautzInner{E,T}) where {E,T} =
    CompositeInner(vcat(H2.kautz_inners, H1.kautz_inners))

"""
    poles(H::AbstractInnerFunction)

Computes the poles of the inner function H.
"""
poles(L::LaguerreInner) = [complex(-L.λ)]

function poles(K::KautzInner{T}) where {T}
    if K.α > K.β
        φ = one(complex(T)) * sqrt((K.α + K.β) * (K.α - K.β))
    elseif K.α < K.β
        φ = im * sqrt((K.β + K.α) * (K.β - K.α))
    elseif K.α == K.β
        φ = one(complex(T)) * zero(K.α)
    end
    return [-K.α - φ, -K.α + φ]
end

poles(C::CompositeLaguerreInner) = reduce(vcat, poles.(C.laguerre_inners))
poles(C::CompositeKautzInner) = reduce(vcat, poles.(C.kautz_inners))
function poles(C::CompositeInner)
    pl = reduce(vcat, poles.(C.laguerre_inners))
    pk = reduce(vcat, poles.(C.kautz_inners))
    return vcat(pl, pk)
end

"""
    shift_basis(H::AbstractInnerFunction, n::Int, ω::Number)

Computes the shift basis functions up to order n evaluated at im*ω.
"""
function shift_basis!(es::AbstractVector, H::AbstractInnerFunction, n::Integer, z::Complex)

    (; A, B, C, D) = ssrealize(H)
    B = view(B, :, 1)

    ne = nstates(H)
    hz = H(z)
    ldiv!(view(es, 1:ne), lu!(z * I - A), B)
    @simd ivdep for i = 1:n
        @inbounds @views es[i*ne+1:(i+1)*ne] .= hz * es[(i-1)*ne+1:i*ne]
    end
    return es
end

function shift_basis(H::AbstractInnerFunction, n::Integer, z::Complex)
    es = zeros(typeof(z), (n + 1) * nstates(H))
    shift_basis!(es, H, n, z)
    return es
end
