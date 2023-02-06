using StateSpaceRealizables, LinearAlgebra
using Test

include("test_state_space.jl")
include("test_inner_functions.jl")

@testset "StateSpaceRealizables.jl" begin
    test_state_space()
    test_inner_functions()
end
