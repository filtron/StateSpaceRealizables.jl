"""
    StateSpace{E,T,TA,TB,TC,TD}

Type for representing state-space models, i.e:

```math
dx = A x + B du,
y = C x + D u,
```

where d is a difference operator if E<:DiscreteTE,
or a differential if E<:ContinuousTE.
"""
struct StateSpace{E,T,TA,TB,TC,TD} <: AbstractStateSpace{E,T}
    A::TA
    B::TB
    C::TC
    D::TD
end

"""
StateSpace{E,T}(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix) where {E<:TimeEvolution, T<:Number}

Constructs a state-space model of time evolution E and numeric type T from the matrices A, B, C, D.
"""
function StateSpace{E,T}(
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
    D::AbstractMatrix,
) where {E<:TimeEvolution,T<:Number}
    A, B, C, D = convert.(AbstractMatrix{T}, (A, B, C, D))
    return StateSpace{E,T,typeof(A),typeof(B),typeof(C),typeof(D)}(A, B, C, D)
end

"""
StateSpace{E}(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix) where {E<:TimeEvolution, T<:Number}

Constructs a state-space model of time evolution E from the matrices A, B, C, D.
The inputs are converted to a common numeric type.
"""
function StateSpace{E}(
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
    D::AbstractMatrix,
) where {E<:TimeEvolution}
    T = promote_type(eltype(A), eltype(B), eltype(C), eltype(D))
    return StateSpace{E,T}(A, B, C, D)
end

ssparameter_type(
    ::StateSpace{ET,T,<:AbstractMatrix,<:AbstractMatrix,<:AbstractMatrix,<:AbstractMatrix},
) where {ET,T} = MatrixParameter()



ssparams(ss::StateSpace) = ss.A, ss.B, ss.C, ss.D

ssrealize(s::StateSpace) = s
ssrealize(A, B, C, D; te::TimeEvolution = ContinuousTE()) =
    StateSpace{typeof(te)}(A, B, C, D)

ssnegate(s::StateSpace{E}, ::MatrixParameter) where {E} =
    StateSpace{E}(s.A, s.B, -s.C, -s.D)


function _sqr_magnitude_response!(
        mag2::AbstractVector{<:AbstractMatrix},
        s::StateSpace{<:ContinuousTE},
        ::MatrixParameter,
        zs,
    )
        A, B, C, D = ssparams(s)
        E = complex.(similar(B))
        H = complex.(similar(D))
        for (i, z) in enumerate(zs)
            ldiv!(E, lu!(z * I - A), B)
            mul!(H, C, E)
            H .= D + H
            mag2[i] .= H' * H
        end
end


function _impulse_response!(
    hs::AbstractVector{<:AbstractMatrix},
    s::StateSpace{<:ContinuousTE},
    ::MatrixParameter,
    ts,
)
    A, B, C, D = ssparams(s)
    E = similar(B)
    for (i, t) in enumerate(ts)
        mul!(E, LinearAlgebra.exp!(A * t), B)
        mul!(hs[i], C, E)
        hs[i] .= hs[i] + D
    end
end