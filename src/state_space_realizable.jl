"""
    StateSpaceRealizable{E<:TimeEvolution}

Abstract type for representing systems that may be realized as a state-space model.
"""
abstract type StateSpaceRealizable{E<:TimeEvolution} end

time_evolution(::StateSpaceRealizable{E}) where {E} = E

# trait
abstract type Properness end
struct IsProper end
struct IsNotProper end
function isproper end


abstract type LTIStateSpaceRealizable{E} <: StateSpaceRealizable{E} end


"""
    AbstractStateSpace{E,T<:Number}

Abstract type for represinting systems in state-space form.
"""
abstract type AbstractStateSpace{E,T<:Number} <: LTIStateSpaceRealizable{E} end

#=
abstract type StateSpaceParameterType end
struct MatrixParameter <: StateSpaceParameterType end
=#

isproper(s::AbstractStateSpace) = iszero(s.D)


function _sqr_magnitude_response(s::AbstractStateSpace, zs)
    n = length(zs)
    m = ninputs(s)
    mag2 = [zeros(eltype(zs), m, m) for i = 1:n]
    sqr_magnitude_response!(mag2, s, zs)
    return mag2
end

function _impulse_response(s::AbstractStateSpace, ts)
    n = length(ts)
    hs = [zeros(eltype(ts), noutputs(s), ninputs(s)) for i = 1:n]
    impulse_response!(hs, s, ts)
    return hs
end