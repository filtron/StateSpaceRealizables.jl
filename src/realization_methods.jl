abstract type RealizationMethod end

"""
    AnyRealization <: RealizationMethod

Type for describing an uknown state-space realization method 
"""
struct AnyRealization <: RealizationMethod end

"""
    Balanced <: RealizationMethod

Type for describing balanced realizations, i.e.
reachability and controllability Grammians are equal.  
"""
struct Balanced <: RealizationMethod end

"""
    InternallyBalanced <: RealizationMethod

Type for denoting internally balanced realizations, i.e. 
reachability Grammian is identity. 
"""
struct InternallyBalanced <: RealizationMethod end
