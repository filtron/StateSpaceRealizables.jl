"""
    ProperStateSpace{E,T,TA,TB,TC}

Type for representing proper state-space models, i.e:

```math
dx = A x + B du,
y = C x,
```

where d is a difference operator if E<:DiscreteTE,
or a differential if E<:ContinuousTE.
"""
struct ProperStateSpace{E,T,TA,TB,TC} <: AbstractStateSpace{E,T}
    A::TA
    B::TB
    C::TC
end

function ProperStateSpace{E,T}(
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
) where {E<:TimeEvolution,T<:Number}
    A, B, C = convert.(AbstractMatrix{T}, (A, B, C))
    return ProperStateSpace{E,T,typeof(A),typeof(B),typeof(C)}(A, B, C)
end

function ProperStateSpace{E}(
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
) where {E<:TimeEvolution}
    T = promote_type(eltype(A), eltype(B), eltype(C))
    return ProperStateSpace{E,T}(A, B, C)
end


ssparameter_type(
    ::ProperStateSpace{ET,T,<:AbstractMatrix,<:AbstractMatrix,<:AbstractMatrix},
) where {ET,T} = MatrixParameter()

isproper(::ProperStateSpace) = true


"""
    ssparams(s::StateSpaceRealizable)

Computes the state-space parameters A, B, C, D of the
state-space model s.
For proper state-space models, the return value of D is nothing.
"""
ssparams(ss::ProperStateSpace) = ss.A, ss.B, ss.C, nothing


"""
    ssrealize(s::StateSpaceRealizable)

Computes a state-space realization of the state-space model s.
Returns the input if it is already of a state-space realized type,
e.g. StateSpace / ProperStateSpace.
"""
ssrealize(s::ProperStateSpace) = s
ssrealize(A, B, C; te::TimeEvolution = ContinuousTE()) =
    ProperStateSpace{typeof(te)}(A, B, C)


function _sqr_magnitude_response!(
        mag2::AbstractVector{<:AbstractMatrix},
        s::ProperStateSpace{<:ContinuousTE},
        ::MatrixParameter,
        zs,
    )
        A, B, C, D = ssparams(s)
        E = complex.(similar(B))
        H = zeros(eltype(zs), noutputs(s), ninputs(s))
        for (i, z) in enumerate(zs)
            ldiv!(E, lu!(z * I - A), B)
            mul!(H, C, E)
            H .= H
            mag2[i] .= H' * H
        end
end

function _impulse_response!(
    hs::AbstractVector{<:AbstractMatrix},
    s::ProperStateSpace{<:ContinuousTE},
    ::MatrixParameter,
    ts,
)
    A, B, C, D = ssparams(s)
    E = similar(B)
    for (i, t) in enumerate(ts)
        mul!(E, LinearAlgebra.exp!(A * t), B)
        mul!(hs[i], C, E)
    end
end