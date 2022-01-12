## util.jl
module util

export partTypeNum

function partTypeNum(partType)

	if isa(partType, Integer)
		return partType
	end

	if lowercase(partType) in ["gas", "cells"]
		return 0
	end

	if lowercase(partType) in ["dm", "darkmatter"]
		return 1
	end

	if lowercase(partType) in ["tracer", "tracers", "tracermc", "trmc"]
		return 3
	end

	if lowercase(partType) in ["star", "stars", "stellar"]
		return 4
	end

	if lowercase(partType) in ["wind"]
		return 4
	end

	if lowercase(partType) in ["bh", "bhs", "blackhole", "blackholes"]
		return 5
	end

    error("Unknown Particle Type")

end

end
