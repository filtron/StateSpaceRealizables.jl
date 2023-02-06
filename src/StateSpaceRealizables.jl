module StateSpaceRealizables

using LinearAlgebra, ArrayInterface

import Base: ∘, *, +, -, ^, getproperty, eltype, length, size, ==, isequal, isapprox

include("time_evolution.jl")
export TimeEvolution, ContinuousTE, DiscreteTE

include("realization_methods.jl")
export RealizationMethod, AnyRealization, Balanced, InternallyBalanced

include("state_space_realizable.jl")
export StateSpaceRealizable,
    time_evolution,
    AbstractStateSpace,
    StateSpace,
    ProperStateSpace,
    ninputs,
    nstates,
    noutputs,
    isproper,
    poles,
    ssparams,
    ssrealize,
    sqr_magnitude_response,
    sqr_magnitude_response!,
    impulse_response,
    impulse_response!

include("inner_functions.jl")
export AbstractInnerFunction,
    UnivariateInnerFunction,
    LaguerreInner,
    KautzInner,
    CompositeInner,
    CompositeLaguerreInner,
    CompositeKautzInner,
    butter_inner,
    shift_basis!,
    shift_basis!

end
