"""
    StateSpaceRealizable{E<:TimeEvolution}

Abstract type for representing systems that may be realized as a state-space model.
"""
abstract type StateSpaceRealizable{E<:TimeEvolution} end


time_evolution(::StateSpaceRealizable{E}) where {E} = E

# also add Realization as a type parameter.

"""
    AbstractStateSpace{E,T<:Number}

Abstract type for represinting systems in state-space form.
"""
abstract type AbstractStateSpace{E,T<:Number} <: StateSpaceRealizable{E} end

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

abstract type StateSpaceParameterType end
struct MatrixParameter end

ssparameter_type(
    ::StateSpace{ET,T,<:AbstractMatrix,<:AbstractMatrix,<:AbstractMatrix,<:AbstractMatrix},
) where {ET,T} = MatrixParameter()
ssparameter_type(
    ::ProperStateSpace{ET,T,<:AbstractMatrix,<:AbstractMatrix,<:AbstractMatrix},
) where {ET,T} = MatrixParameter()

ninputs(s::AbstractStateSpace) = ninputs(s, ssparameter_type(s))
nstates(s::AbstractStateSpace) = nstates(s, ssparameter_type(s))
noutputs(s::AbstractStateSpace) = noutputs(s, ssparameter_type(s))

ninputs(s::AbstractStateSpace, ::MatrixParameter) = size(s.B, 2)
nstates(s::AbstractStateSpace, ::MatrixParameter) = size(s.B, 1)
noutputs(s::AbstractStateSpace, ::MatrixParameter) = size(s.C, 1)

isproper(s::AbstractStateSpace) = iszero(s.D)
isproper(::ProperStateSpace) = true

poles(s::AbstractStateSpace) = poles(s, ssparameter_type(s))
poles(s::AbstractStateSpace, ::MatrixParameter) = complex.(eigvals(s.A))


"""
    ssparams(s::StateSpaceRealizable)

Computes the state-space parameters A, B, C, D of the
state-space model s.
For proper state-space models, the return value of D is nothing.
"""
ssparams(ss::StateSpace) = ss.A, ss.B, ss.C, ss.D
ssparams(ss::ProperStateSpace) = ss.A, ss.B, ss.C, nothing

"""
    ssrealize(s::StateSpaceRealizable)

Computes a state-space realization of the state-space model s.
Returns the input if it is already of a state-space realized type,
e.g. StateSpace / ProperStateSpace.
"""
ssrealize(s::StateSpace) = s
ssrealize(s::ProperStateSpace) = s

ssrealize(A, B, C, D; te::TimeEvolution = ContinuousTE()) =
    StateSpace{typeof(te)}(A, B, C, D)
ssrealize(A, B, C; te::TimeEvolution = ContinuousTE()) =
    ProperStateSpace{typeof(te)}(A, B, C)



+(s1::AbstractStateSpace{E,T}, s2::AbstractStateSpace{E,T}) where {E,T} =
    ssadd(s1, ssparameter_type(s1), s2, ssparameter_type(s2))

function ssadd(
    s1::StateSpace{E,T},
    ::MatrixParameter,
    s2::StateSpace{E,T},
    ::MatrixParameter,
) where {E,T}
    A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
    B = vcat(s2.B, s1.B)
    C = hcat(s2.C, s1.C)
    D = s2.D + s1.D
    return StateSpace{E}(A, B, C, D)
end

function ssadd(
    s1::StateSpace{E,T},
    ::MatrixParameter,
    s2::ProperStateSpace{E,T},
    ::MatrixParameter,
) where {E,T}
    A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
    B = vcat(s2.B, s1.B)
    C = hcat(s2.C, s1.C)
    D = s1.D
    return StateSpace{E}(A, B, C, D)
end

function ssadd(
    s1::ProperStateSpace{E,T},
    ::MatrixParameter,
    s2::StateSpace{E,T},
    ::MatrixParameter,
) where {E,T}
    A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
    B = vcat(s2.B, s1.B)
    C = hcat(s2.C, s1.C)
    D = s2.D
    return StateSpace{E}(A, B, C, D)
end

function ssadd(
    s1::ProperStateSpace{E,T},
    ::MatrixParameter,
    s2::ProperStateSpace{E,T},
    ::MatrixParameter,
) where {E,T}
    A = [s2.A zeros(T, nstates(s2), nstates(s1)); zeros(nstates(s1), nstates(s2)) s1.A]
    B = vcat(s2.B, s1.B)
    C = hcat(s2.C, s1.C)
    return ProperStateSpace{E}(A, B, C, D)
end


-(s::AbstractStateSpace) = ssnegate(s, ssparameter_type(s))
-(s1::AbstractStateSpace{E,T}, s2::AbstractStateSpace{E,T}) where {E,T} = +(s1, -s2)

ssnegate(s::StateSpace{E}, ::MatrixParameter) where {E} =
    StateSpace{E}(s.A, s.B, -s.C, -s.D)
ssnegate(s::ProperStateSpace{E}, ::MatrixParameter) where {E} =
    ProperStateSpace{E}(s.A, s.B, -s.C)


∘(s2::AbstractStateSpace{E,T}, s1::AbstractStateSpace{E,T}) where {E,T} =
    sscompose(s2, ssparameter_type(s2), s1, ssparameter_type(s1)) # prefer this one
*(s2::AbstractStateSpace{E,T}, s1::AbstractStateSpace{E,T}) where {E,T} = ∘(s2, s1)  # keep this one to not break things

function sscompose(
    s2::StateSpace{E,T},
    ::MatrixParameter,
    s1::StateSpace{E,T},
    ::MatrixParameter,
) where {E,T}
    A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
    B = vcat(s1.B, s2.B * s1.D)
    C = hcat(s2.D * s1.C, s2.C)
    D = s2.D * s1.D
    return StateSpace{E}(A, B, C, D)
end

function sscompose(
    s2::StateSpace{E,T},
    ::MatrixParameter,
    s1::ProperStateSpace{E,T},
    ::MatrixParameter,
) where {E,T}
    A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
    B = vcat(s1.B, zeros(T, nstates(s2), ninputs(s1)))
    C = hcat(s2.D * s1.C, s2.C)
    return ProperStateSpace{E}(A, B, C)
end

function sscompose(
    s2::ProperStateSpace{E,T},
    ::MatrixParameter,
    s1::StateSpace{E,T},
    ::MatrixParameter,
) where {E,T}
    A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
    B = vcat(s1.B, s2.B * s1.D)
    C = hcat(zero(s1.C), s2.C)
    return ProperStateSpace{E}(A, B, C)
end

function sscompose(
    s2::ProperStateSpace{E,T},
    ::MatrixParameter,
    s1::ProperStateSpace{E,T},
    ::MatrixParameter,
) where {E,T}
    A = [s1.A zeros(T, nstates(s1), nstates(s2)); s2.B*s1.C s2.A]
    B = vcat(s1.B, zeros(T, nstates(s2), ninputs(s2)))
    C = hcat(zeros(T, noutputs(s1), nstates(s1)), s2.C)
    return ProperStateSpace{E}(A, B, C)
end




sqr_magnitude_response!(mag2, s::AbstractStateSpace, zs) =
    _sqr_magnitude_response!(mag2, s, ssparameter_type(s), zs)

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

sqr_magnitude_response(s::AbstractStateSpace, zs) =
    _sqr_magnitude_response(s::AbstractStateSpace, ssparameter_type(s), zs)

function _sqr_magnitude_response(s::AbstractStateSpace, ::MatrixParameter, zs)
    n = length(zs)
    m = ninputs(s)
    mag2 = [zeros(eltype(zs), m, m) for i = 1:n]
    sqr_magnitude_response!(mag2, s, zs)
    return mag2
end


impulse_response!(hs, s::AbstractStateSpace, ts) =
    _impulse_response!(hs, s, ssparameter_type(s), ts)

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

impulse_response(s::AbstractStateSpace, ts) = _impulse_response(s, ssparameter_type(s), ts)

function _impulse_response(s::AbstractStateSpace, ::MatrixParameter, ts)
    n = length(ts)
    m = noutputs(s)
    hs = [zeros(eltype(ts), noutputs(s), ninputs(s)) for i = 1:n]
    impulse_response!(hs, s, ts)
    return hs
end