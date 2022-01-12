## sublink.jl

module sublink

using HDF5
using Glob

include("groupcat.jl")
include("util.jl")

export treePath
export treeOffsets
export subLinkOffsets
export loadTree
export maxPastMass
export numMergers

function treePath(basePath, treeName, chunkNum=0)
    
    fileName = "tree_extended." * string(chunkNum) * ".hdf5"
    treePath = "/trees/" * string(treeName) * "/"
    filePath1 = basePath * treePath
    filePath2 = basePath * "/postprocessing" * treePath
    filePath3 = basePath * "/../postprocessing" * treePath
    
    if size( glob(fileName, filePath1) )[end] > 0
        return filePath1 * fileName
    end
    
    if size( glob(fileName, filePath2) )[end] > 0
        return filePath2 * fileName
    end
    
    if size( glob(fileName, filePath3) )[end] > 0
        return filePath3 * fileName
    end
    
    error("Could not construct treePath from the supplied input.")
    
end

function treeOffsets(basePath, snapNum, i, treeName)
    
    if (occursin("fof_subhalo", groupcat.gcPath(basePath, snapNum))) || treeName == "SubLink_gal"
        
        f = h5open(groupcat.offsetPath(basePath, snapNum), "r")
        groupFileOffsets = f["FileOffsets/Subhalo"][:]
        close(f)
        
        offsetFile = groupcat.offsetPath(basePath, snapNum)
        prefix = "Subhalo/" * treeName * "/"
        
        groupOffset = i
        
    else
        
        header = h5readattr(groupcat.gcPath(basePath, snapNum), "Header")
        groupFileOffsets = header["FileOffsets_Subhalo"]
        
        groupFileOffsets = Int64(i) .- groupFileOffsets
        fileNum = maximum( findall(x->x>=0, groupFileOffsets) )
        groupOffset = groupFileOffsets[fileNum]
        
        offsetFile = groupcat.gcPath(basePath, snapNum, fileNum-1)
        prefix = "Offsets/Subhalo_Sublink"
        
    end
    
    f = h5open(offsetFile, "r")
    RowNum = f[prefix * "RowNum"][groupOffset+1]
    LastProgID = f[prefix * "LastProgenitorID"][groupOffset+1]
    SubhaloID = f[prefix * "SubhaloID"][groupOffset+1]
    close(f)
    
    return RowNum, LastProgID, SubhaloID
    
    
end    

offsetCache = Dict()

function subLinkOffsets(basePath, treeName, cache=true)
    
    if cache == true
        cache = offsetCache
    end
    
    if isa(cache, Dict) == true
        path = basePath * "/" * treeName
        
        if haskey(cache, path) == true
            return cache[path]
        end
        
    end    
    
    search_path = treePath(basePath, treeName, "*")
    trun_search_path = chop(search_path, tail = 20)
    fileName = "tree_extended.*.hdf5"
    numTreeFiles = size( glob(fileName, trun_search_path) )[end]
    if numTreeFiles == 0
        error("No tree files found!")
    end
    offsets = zeros(Int64, numTreeFiles)
    
    for i in 1:numTreeFiles-1
        f = h5open(treePath(basePath, treeName, i-1), "r")
        offsets[i+1] = offsets[i] + size(f["SubhaloID"][:])[end]
        close(f)
    end
    
    if isa(cache, Dict) == true
        cache[path] = offsets
    end
    
    return offsets

end

function loadTree(basePath, snapNum, i, fields=nothing, onlyMPB=false, treeName="SubLink", cache=true)
    
    RowNum, LastProgID, SubhaloID = treeOffsets(basePath, snapNum, i, treeName)
    
    if RowNum == -1
        print("Warning, empty return")
        return nothing
    end
    
    rowStart = RowNum
    rowEnd = RowNum + (LastProgID - SubhaloID)
    nRows = rowEnd - rowStart + 1
    
    if typeof(fields) == String
        fields = [fields]
    end
    
    offsets = subLinkOffsets(basePath, treeName, cache)
    
    rowOffsets = rowStart .- offsets
    
    fileNum = maximum( findall(x->x>=0, rowOffsets) ) - 1
    fileOff = rowOffsets[fileNum+1]
    
    if onlyMPB == true
        f = h5open(treePath(basePath, treeName, fileNum), "r")
        MainLeafProgenitorID = f["MainLeafProgenitorID"][fileOff+1]
        close(f)
        
        rowEnd = RowNum + (MainLeafProgenitorID - SubhaloID)
        nRows = rowEnd - rowStart + 1
        
    end
    
    result = Dict()
    
    f = h5open(treePath(basePath, treeName, fileNum), "r")
    
    if isnothing(fields)
		fields = keys(f)
	end
    
    if fileOff + nRows > size(f["SubfindID"][:])[end]
        error("Should not occur. Each tree is contained within a single file.")
    end
    
    for field in fields
        
        if( field in keys(f) ) == false
            error("SubLink tree does not have field [" * field * "]")
        end
        
        shape = [i for i in size(f[field])]
        
        if length(shape) == 1
            result[field] = f[field][fileOff+1:fileOff+nRows]
        else
            result[field] = f[field][:,fileOff+1:fileOff+nRows]
        end
        
    end
    
    if length(fields) == 1
		return result[fields[1]]
	end
    
    return result
    
end

function maxPastMass(tree, index, partType="stars")
    
    ptNum = util.partTypeNum(partType)
    
    branchSize = tree["MainLeafProgenitorID"][index+1] - tree["SubhaloID"][index+1] + 1
    masses = tree["SubhaloMassType"][ptNum+1, index+1:index+branchSize]
    return maximum(masses)
        
end

function numMergers(tree, minMassRatio=1e-10, massPartType="stars", index=0)
    
    reqFields = ["SubhaloID", "NextProgenitorID", "MainLeafProgenitorID", "FirstProgenitorID", "SubhaloMassType"]
    
    for field in reqFields
        
        if (field in keys(tree)) == false
            
            error("Error: Input tree does not contain the field: " * field)
            
        end
        
    end
    
    numMergers = 0
    invMassRatio = 1.0 / minMassRatio
    
    rootID = tree["SubhaloID"][index+1]
    fpID = tree["FirstProgenitorID"][index+1]
    
    while fpID != -1
        
        fpIndex = index + (fpID - rootID)
        fpMass  = maxPastMass(tree, fpIndex, massPartType)
        
        npID = tree["NextProgenitorID"][fpIndex+1]
        
        while npID != -1
            
            npIndex = index + (npID - rootID)
            npMass  = maxPastMass(tree, npIndex, massPartType)
            
            if (fpMass > 0.0) & (npMass > 0.0)
                
                ratio = npMass / fpMass
                
                if (ratio >= minMassRatio) & (ratio <=invMassRatio)
                    
                    numMergers = numMergers + 1
                    
                end
                
            end
            
            npID = tree["NextProgenitorID"][npIndex+1]
            
        end
        
        fpID = tree["FirstProgenitorID"][fpIndex+1]
        
    end
    
    return numMergers
            
end

end
