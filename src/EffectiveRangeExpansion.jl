module EffectiveRangeExpansion

using QuadGK
using Parameters

import Base:show

export breakuprotatecut
include("kinematics.jl")

export circleintegral
export cauchyintegral, cauchyintegral′, cauchyintegral′′
export CircularIntegralMethod
export CircularIntegral, CircularSum
include("integrals.jl")

export EffectiveRangeExpansionMethod, ComplexBranchPointExpansion
export effectiverangeexpansion
include("expansion.jl")

export ERP
include("parameters.jl")

export tophysicsunits
include("conversion.jl")


end