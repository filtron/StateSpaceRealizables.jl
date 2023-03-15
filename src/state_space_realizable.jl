"""
    StateSpaceRealizable{E<:TimeEvolution}

Abstract type for representing systems that may be realized as a state-space model.
"""
abstract type StateSpaceRealizable{E<:TimeEvolution} end


time_evolution(::StateSpaceRealizable{E}) where {E} = E


"""
    AbstractStateSpace{E,T<:Number}

Abstract type for represinting systems in state-space form.
"""
abstract type AbstractStateSpace{E,T<:Number} <: StateSpaceRealizable{E} end


abstract type StateSpaceParameterType end
struct MatrixParameter end

ninputs(s::AbstractStateSpace) = ninputs(s, ssparameter_type(s))
nstates(s::AbstractStateSpace) = nstates(s, ssparameter_type(s))
noutputs(s::AbstractStateSpace) = noutputs(s, ssparameter_type(s))

ninputs(s::AbstractStateSpace, ::MatrixParameter) = size(s.B, 2)
nstates(s::AbstractStateSpace, ::MatrixParameter) = size(s.B, 1)
noutputs(s::AbstractStateSpace, ::MatrixParameter) = size(s.C, 1)

isproper(s::AbstractStateSpace) = iszero(s.D)

poles(s::AbstractStateSpace) = poles(s, ssparameter_type(s))
poles(s::AbstractStateSpace, ::MatrixParameter) = complex.(eigvals(s.A))

sqr_magnitude_response!(mag2, s::AbstractStateSpace, zs) =
    _sqr_magnitude_response!(mag2, s, ssparameter_type(s), zs)

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

impulse_response(s::AbstractStateSpace, ts) = _impulse_response(s, ssparameter_type(s), ts)

function _impulse_response(s::AbstractStateSpace, ::MatrixParameter, ts)
    n = length(ts)
    m = noutputs(s)
    hs = [zeros(eltype(ts), noutputs(s), ninputs(s)) for i = 1:n]
    impulse_response!(hs, s, ts)
    return hs
end