"""
    TimeEvolution

Abstract type denoting the time evolution of a dynamical system.
See also ContinuousTime, DiscreteTime.
"""
abstract type TimeEvolution end

"""
    ContinuousTE <: TimeEvolution

Type denoting that a system evolves in continuous-time.
"""
struct ContinuousTE <: TimeEvolution end


# ADD ?:
# abstract type AbstractDiscreteTime <: TImeEvolution end
# struct SampledTE <: AbstractTimeEvolution end

"""
    DiscreteTE <: TimeEvolution

Type denoting that a system evolves in discrete-time.
"""
struct DiscreteTE <: TimeEvolution end
