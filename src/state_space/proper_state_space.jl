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

const MatrixProperStateSpace{ET,T,TA,TB,TC} = ProperStateSpace{ET,T,TA,TB,TC} where {ET,T,TA<:AbstractMatrix,TB<:AbstractMatrix,TC<:AbstractMatrix}

ninputs(s::MatrixProperStateSpace) = size(s.B, 2)
nstates(s::MatrixProperStateSpace) = size(s.A, 1)
noutputs(s::MatrixProperStateSpace) = size(s.C, 1)



isproper(::ProperStateSpace) = IsProper()


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


poles(s::MatrixProperStateSpace) =  complex.(eigvals(s.A))


-(s::MatrixProperStateSpace{E}) where {E} =
    ProperStateSpace{E}(s.A, s.B, -s.C)



function sqr_magnitude_response!(
        mag2::AbstractVector{<:AbstractMatrix},
        s::MatrixProperStateSpace{<:ContinuousTE},
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

function impulse_response!(
    hs::AbstractVector{<:AbstractMatrix},
    s::MatrixProperStateSpace{<:ContinuousTE},
    ts,
)
    A, B, C, D = ssparams(s)
    E = similar(B)
    for (i, t) in enumerate(ts)
        mul!(E, LinearAlgebra.exp!(A * t), B)
        mul!(hs[i], C, E)
    end
end