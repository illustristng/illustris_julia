## snapshot.jl

module snapshot

using HDF5

include("groupcat.jl")
include("util.jl")

export snapPath
export getNumPart
export loadSubset
export getSnapOffsets
export loadSubhalo
export loadHalo

function snapPath(basePath, snapNum, chunkNum=0)

    snapPath = basePath * "/snapdir_" * lpad(snapNum, 3, "0") * "/"
    filePath = snapPath * "snap_" * lpad(snapNum, 3, "0") * "." * lpad(chunkNum, 0, "0") * ".hdf5"

    return filePath

end

function getNumPart(header)

    nTypes = 6

    nPart = zeros(Float64, nTypes)

    @. nPart = header["NumPart_Total"] | (header["NumPart_Total_HighWord"] << 32)

    return nPart

end

function loadSubset(basePath, snapNum, partType, fields = nothing, subset = nothing)

    result = Dict()

    ptNum = util.partTypeNum(partType)
    gName = "PartType" * string(ptNum)

    if typeof(fields) == String
        fields = [fields]
    end

    f = h5open(snapPath(basePath, snapNum), "r")

    header = h5readattr(snapPath(basePath, snapNum), "Header")
    nPart = getNumPart(header)

    if isnothing(subset) == false

        offsetsThisType = subset["offsetType"][ptNum+1] .- subset["snapOffsets"][:, ptNum+1]

        fileNum = maximum( findall(x->x>=0, offsetsThisType) ) - 1
        fileOff = offsetsThisType[fileNum+1]
        numToRead = subset["lenType"][ptNum+1]

    else

        fileNum = 0
        fileOff = 0
        numToRead = nPart[ptNum+1]

    end

    result["count"] = numToRead

    if (numToRead == 0)
		println("warning: no particles of requested type, empty return.")
		return result
	end
    
    i = 1
    while (gName in keys(f)) == false
        f = h5open(snapPath(basePath, snapNum, i), "r")
        i = i + 1
    end
    
    if isnothing(fields)
		fields = keys(f[gName])
	end
    
    for field in fields
        
        if( field in keys(f[gName]) ) == false
            error("Particle Type [" * string(ptNum) * "] does not have field [" * field * "]")
        end
        
        shape = [i for i in size(f[gName][field])]
        shape[end] = numToRead
        
        result[field] = zeros(eltype(f[gName][field]), Tuple(shape))
        
    end
    
    close(f)
    
    wOffset = 1
    origNumToRead = numToRead
    
    while (numToRead > 0)
        f = h5open(snapPath(basePath, snapNum, fileNum), "r")
        
        if (gName in keys(f)) == false
            close(f)
            fileNum = fileNum + 1
            fileOff = 0
            continue
        end
        
        header_l = h5readattr(snapPath(basePath, snapNum, fileNum), "Header")
        
        numTypeLocal = header_l["NumPart_ThisFile"][ptNum+1]
        
        if isnothing(subset) == true
            
            numToReadLocal = numTypeLocal
        else
            
            if(numToRead > numTypeLocal)
                numToReadLocal = numTypeLocal
            else
                numToReadLocal = numToRead
            end
        end
        
        if (fileOff + numToReadLocal > numTypeLocal)
            numToReadLocal = numTypeLocal - fileOff
        end
        
        for field in fields
            
            shape = [i for i in size(f[gName][field])]
            
            if length(shape) == 1
                result[field][wOffset:wOffset+numToReadLocal-1] = f[gName][field][fileOff+1:fileOff+numToReadLocal]
            else
                result[field][:,wOffset:wOffset+numToReadLocal-1] = f[gName][field][:,fileOff+1:fileOff+numToReadLocal]
            end
            
        end
                
        wOffset = wOffset + numToReadLocal
        numToRead = numToRead - numToReadLocal
        fileNum = fileNum + 1
        fileOff = 0
        
        close(f)
        
    end
    
    if origNumToRead != (wOffset-1)
        error("Read [" * string(wOffset-1) * "] particles, but was expecting [" * string(origNumToRead) * "]")
        
    end
    
    if length(fields) == 1
		return result[fields[1]]
	end
    
    return result
        
end

function getSnapOffsets(basePath, snapNum, i, t)
    
    r = Dict()
    
    if occursin("fof_subhalo", groupcat.gcPath(basePath, snapNum))
        f = h5open(groupcat.offsetPath(basePath, snapNum), "r")
        groupFileOffsets = f["FileOffsets/" * t][:]
        r["snapOffsets"] = transpose(f["FileOffsets/SnapByType"][:,:])
        close(f)
    else
        header_l = h5readattr(groupcat.gcPath(basePath, snapNum), "Header")
        groupFileOffsets = header_l["FileOffsets_" * t]
        r["snapOffsets"] = header_l["FileOffsets_Snap"]
    end    
    
    groupFileOffsets = Int64(i) .- groupFileOffsets
    fileNum = maximum( findall(x->x>=0, groupFileOffsets) )
    groupOffset = groupFileOffsets[fileNum]
    
    f = h5open(groupcat.gcPath(basePath, snapNum, fileNum-1), "r")
    r["lenType"] = f[t][t * "LenType"][:, groupOffset+1]
    close(f)
    
    if occursin("fof_subhalo", groupcat.gcPath(basePath, snapNum))
        f = h5open(groupcat.offsetPath(basePath, snapNum), "r")
        r["offsetType"] = f[t * "/SnapByType"][:, i+1]
        close(f)
    else
        f = h5open(groupcat.gcPath(basePath, snapNum, fileNum), "r")
        r["offsetType"] = f["Offsets"][t * "_SnapByType"][:, groupOffset+1]
    end
    
    return r 
    
end 


function loadSubhalo(basePath, snapNum, i, partType, fields = nothing)
    subset = getSnapOffsets(basePath, snapNum, i, "Subhalo")
    result = Dict()

    ptNum = util.partTypeNum(partType)
    gName = "PartType" * string(ptNum)

    if typeof(fields) == String
        fields = [fields]
    end

    f = h5open(snapPath(basePath, snapNum), "r")

    header = h5readattr(snapPath(basePath, snapNum), "Header")
    nPart = getNumPart(header)

    if isnothing(subset) == false

        offsetsThisType = subset["offsetType"][ptNum+1] .- subset["snapOffsets"][:, ptNum+1]

        fileNum = maximum( findall(x->x>=0, offsetsThisType) ) - 1
        fileOff = offsetsThisType[fileNum+1]
        numToRead = subset["lenType"][ptNum+1]

    else

        fileNum = 0
        fileOff = 0
        numToRead = nPart[ptNum+1]

    end

    result["count"] = numToRead

    if (numToRead == 0)
		println("warning: no particles of requested type, empty return.")
		return result
	end
    
    i = 1
    while (gName in keys(f)) == false
        f = h5open(snapPath(basePath, snapNum, i), "r")
        i = i + 1
    end
    
    if isnothing(fields)
		fields = keys(f[gName])
	end
    
    for field in fields
        
        if( field in keys(f[gName]) ) == false
            error("Particle Type [" * string(ptNum) * "] does not have field [" * field * "]")
        end
        
        shape = [i for i in size(f[gName][field])]
        shape[end] = numToRead
        
        result[field] = zeros(eltype(f[gName][field]), Tuple(shape))
        
    end
    
    close(f)
    
    wOffset = 1
    origNumToRead = numToRead
    
    while (numToRead > 0)
        f = h5open(snapPath(basePath, snapNum, fileNum), "r")
        
        if (gName in keys(f)) == false
            close(f)
            fileNum = fileNum + 1
            fileOff = 0
            continue
        end
        
        header_l = h5readattr(snapPath(basePath, snapNum, fileNum), "Header")
        
        numTypeLocal = header_l["NumPart_ThisFile"][ptNum+1]
        
        if isnothing(subset) == true
            
            numToReadLocal = numTypeLocal
        else
            
            if(numToRead > numTypeLocal)
                numToReadLocal = numTypeLocal
            else
                numToReadLocal = numToRead
            end
        end
        
        if (fileOff + numToReadLocal > numTypeLocal)
            numToReadLocal = numTypeLocal - fileOff
        end
        
        for field in fields
            
            shape = [i for i in size(f[gName][field])]
            
            if length(shape) == 1
                result[field][wOffset:wOffset+numToReadLocal-1] = f[gName][field][fileOff+1:fileOff+numToReadLocal]
            else
                result[field][:,wOffset:wOffset+numToReadLocal-1] = f[gName][field][:,fileOff+1:fileOff+numToReadLocal]
            end
            
        end
                
        wOffset = wOffset + numToReadLocal
        numToRead = numToRead - numToReadLocal
        fileNum = fileNum + 1
        fileOff = 0
        
        close(f)
        
    end
    
    if origNumToRead != (wOffset-1)
        error("Read [" * string(wOffset-1) * "] particles, but was expecting [" * string(origNumToRead) * "]")
        
    end
    
    if length(fields) == 1
		return result[fields[1]]
	end
    
    return result
end

function loadHalo(basePath, snapNum, i, partType, fields = nothing)
    subset = getSnapOffsets(basePath, snapNum, i, "Group")
    result = Dict()

    ptNum = util.partTypeNum(partType)
    gName = "PartType" * string(ptNum)

    if typeof(fields) == String
        fields = [fields]
    end

    f = h5open(snapPath(basePath, snapNum), "r")

    header = h5readattr(snapPath(basePath, snapNum), "Header")
    nPart = getNumPart(header)

    if isnothing(subset) == false

        offsetsThisType = subset["offsetType"][ptNum+1] .- subset["snapOffsets"][:, ptNum+1]

        fileNum = maximum( findall(x->x>=0, offsetsThisType) ) - 1
        fileOff = offsetsThisType[fileNum+1]
        numToRead = subset["lenType"][ptNum+1]

    else

        fileNum = 0
        fileOff = 0
        numToRead = nPart[ptNum+1]

    end

    result["count"] = numToRead

    if (numToRead == 0)
		println("warning: no particles of requested type, empty return.")
		return result
	end
    
    i = 1
    while (gName in keys(f)) == false
        f = h5open(snapPath(basePath, snapNum, i), "r")
        i = i + 1
    end
    
    if isnothing(fields)
		fields = keys(f[gName])
	end
    
    for field in fields
        
        if( field in keys(f[gName]) ) == false
            error("Particle Type [" * string(ptNum) * "] does not have field [" * field * "]")
        end
        
        shape = [i for i in size(f[gName][field])]
        shape[end] = numToRead
        
        result[field] = zeros(eltype(f[gName][field]), Tuple(shape))
        
    end
    
    close(f)
    
    wOffset = 1
    origNumToRead = numToRead
    
    while (numToRead > 0)
        f = h5open(snapPath(basePath, snapNum, fileNum), "r")
        
        if (gName in keys(f)) == false
            close(f)
            fileNum = fileNum + 1
            fileOff = 0
            continue
        end
        
        header_l = h5readattr(snapPath(basePath, snapNum, fileNum), "Header")
        
        numTypeLocal = header_l["NumPart_ThisFile"][ptNum+1]
        
        if isnothing(subset) == true
            
            numToReadLocal = numTypeLocal
        else
            
            if(numToRead > numTypeLocal)
                numToReadLocal = numTypeLocal
            else
                numToReadLocal = numToRead
            end
        end
        
        if (fileOff + numToReadLocal > numTypeLocal)
            numToReadLocal = numTypeLocal - fileOff
        end
        
        for field in fields
            
            shape = [i for i in size(f[gName][field])]
            
            if length(shape) == 1
                result[field][wOffset:wOffset+numToReadLocal-1] = f[gName][field][fileOff+1:fileOff+numToReadLocal]
            else
                result[field][:,wOffset:wOffset+numToReadLocal-1] = f[gName][field][:,fileOff+1:fileOff+numToReadLocal]
            end
            
        end
                
        wOffset = wOffset + numToReadLocal
        numToRead = numToRead - numToReadLocal
        fileNum = fileNum + 1
        fileOff = 0
        
        close(f)
        
    end
    
    if origNumToRead != (wOffset-1)
        error("Read [" * string(wOffset-1) * "] particles, but was expecting [" * string(origNumToRead) * "]")
        
    end
    
    if length(fields) == 1
		return result[fields[1]]
	end
    
    return result
end

end
