module SymbolicDistributions

using Distributions, Symbolics, Polynomials, Expectations, SpecialFunctions

import Distributions: moment

include("moments.jl")
include("linearexpectations.jl")
include("common.jl")

export 
    moment
    linearexpectation

end
