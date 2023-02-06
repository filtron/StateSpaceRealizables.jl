

"""
    time_evolution(::StateSpaceRealizable)

Computes the type of time evolution of the input, i.e.
ContinuousTE or DiscreteTE.
"""
function time_evolution end


# these assume that all possible realization method yield the same dimensionality results.

"""
    ninputs(s::StateSpaceRealizable)

Computes the number of inputs of the state-space realization of s.
"""
function ninputs end

"""
    noutputs(s::StateSpaceRealizable)

Computes the number of outputs of the state-space realization of s.
"""
function noutputs end


"""
    nstates(s::StateSpaceRealizable)

Computes the number of states of the state-space realization of s.
"""
function nstates end


"""
    isproper(s::StateSpaceRealizable)

Computes a Bool indicating whether s may be realized as a proper state-space model.
"""
function isproper end
