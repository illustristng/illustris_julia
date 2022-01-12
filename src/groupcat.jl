## groupcat.jl

module groupcat

using HDF5

export gcPath
export offsetPath
export loadObjects
export loadSubhalos
export loadHalos
export loadHeader
export load
export loadSingle


function gcPath(basePath, snapNum, chunkNum=0)

	gcPath = basePath * "/groups_" * lpad(snapNum, 3, "0") * "/"
	filePath1 = gcPath * "groups_" * lpad(snapNum, 3, "0") * "." * lpad(chunkNum, 0, "0") * ".hdf5"
	filePath2 = gcPath * "fof_subhalo_tab_" * lpad(snapNum, 3, "0") * "." * lpad(chunkNum, 0, "0") * ".hdf5"

	if isfile(filePath1)
		return filePath1
	end
	return filePath2

end

function offsetPath(basePath, snapNum)

	offsetPath = basePath * "/../postprocessing/offsets/offsets_" * lpad(snapNum, 3, "0") * ".hdf5"

end

function loadObjects(basePath, snapNum, gName, nName, fields)

	result = Dict()

	if typeof(fields) == String
		fields = [fields]
	end

	f = h5open(gcPath(basePath, snapNum), "r")

	header = h5readattr(gcPath(basePath, snapNum), "Header")
	result["count"] = header["N" * nName * "_Total"]

	if (result["count"] == 0)
		println("warning: zero groups, empty return (snap=" * string(snapNum) * ").")
		return result
	end

	if isnothing(fields)
		fields = keys(f[gName])
	end

	for field in fields
		if (field in keys(f[gName])) == false
			error("Group Catalog does not have requested field [" * field * "]!")
		end

		shape = [i for i in size(f[gName][field])]
		shape[end] = result["count"]

		result[field] = zeros(eltype(f[gName][field]), Tuple(shape))

	end

	close(f)

	wOffset = 1

	for i in 1:header["NumFiles"]

		f = h5open(gcPath(basePath, snapNum, i-1), "r")

		header_l = h5readattr(gcPath(basePath, snapNum, i-1), "Header")

		if header_l["N" * nName * "_ThisFile"] == 0
			continue
		end

		for field in fields

			if (field in keys(f[gName])) == false
				error("Group Catalog does not have requested field [" * field * "]!")
			end

			shape = [i for i in size(f[gName][field])]

			if length(shape) == 1
				result[field][wOffset:wOffset+shape[end]-1] = f[gName][field][1:shape[end]]
			else
				result[field][:, wOffset:wOffset+shape[end]-1] = f[gName][field][:, 1:shape[end]]
			end

		end
        
        shape = [i for i in size(f[gName][fields[1]])]

		wOffset = wOffset + shape[end]

		close(f)

	end

	if length(fields) == 1
		return result[fields[1]]
	end

	return result

end

function loadSubhalos(basePath, snapNum, fields = nothing)

	return loadObjects(basePath, snapNum, "Subhalo", "subgroups", fields)

end

function loadHalos(basePath, snapNum, fields = nothing)

	return loadObjects(basePath, snapNum, "Group", "groups", fields)

end

function loadHeader(basePath, snapNum)

	header = h5readattr(gcPath(basePath, snapNum), "Header")

end

function load(basePath, snapNum)

	r = Dict()

	r["subhalos"] = loadSubhalos(basePath, snapNum)
	r["halos"] = loadHalos(basePath, snapNum)
	r["header"] = loadHeader(basePath, snapNum)

	return r

end

function loadSingle(basePath, snapNum; haloID=-1, subhaloID=-1)

	if (haloID < 0 && subhaloID < 0) || (haloID >= 0 && subhaloID >= 0)
		error("Must specify either haloID or subhaloID (and not both).")
	end

	if subhaloID >= 0
		gName = "Subhalo"
		searchID = subhaloID
	else
		gName = "Group"
		searchID = haloID
	end

	if occursin("fof_subhalo", gcPath(basePath, snapNum))
		f = h5open(offsetPath(basePath, snapNum), "r")
		offsets = f["FileOffsets/" * gName][:]
		close(f)
	else
		header = h5readattr(gcPath(basePath, snapNum), "Header")
		offsets = header["FileOffsets_" * gName]
	end

	@. offsets = searchID - offsets
	fileNum = maximum( findall(x->x>=0, offsets) )
	groupOffset = offsets[fileNum] + 1

	result = Dict()

	f = h5open(gcPath(basePath, snapNum, fileNum - 1), "r")
	for haloProp in keys(f[gName])
		if ndims(f[gName][haloProp]) == 1
			result[haloProp] = f[gName][haloProp][groupOffset]
		else
			result[haloProp] = f[gName][haloProp][:,groupOffset]
		end
	end
	close(f)

	return result

end

end
