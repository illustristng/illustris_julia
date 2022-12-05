## lhalotree.jl

module lhalotree

using HDF5

include("groupcat.jl")

export treePath
export treeOffsets
export singleNodeFlat
export recProgenitorFlat
export loadTree

function treePath(basePath, chunkNum=0)
    
    filePaths = [ basePath * "/trees/treedata/" * "trees_sf1_135." * string(chunkNum) * ".hdf5", 
                  basePath * "/../postprocessing/trees/LHaloTree/trees_sf1_099." * string(chunkNum) * ".hdf5", 
                  basePath * "/../postprocessing/trees/LHaloTree/trees_sf1_080." * string(chunkNum) * ".hdf5" 
                ]
    
    for filePath in filePaths
        if isfile(filePath)
            return filePath
        end
    end

    throw(ArgumentError("LHaloTree file not found!"))
    
end

function treeOffsets(basePath, snapNum, i)
    
    if occursin("fof_subhalo", groupcat.gcPath(basePath, snapNum))
        
        f = h5open(groupcat.offsetPath(basePath, snapNum), "r")
        groupFileOffsets = f["FileOffsets/Subhalo"][:]
        close(f)
        
        offsetFile = groupcat.offsetPath(basePath, snapNum)
        prefix = "Subhalo/LHaloTree/"
        
        groupOffset = i
       
    else
        
        header = h5readattr(groupcat.gcPath(basePath, snapNum), "Header")
        groupFileOffsets = header["FileOffsets_Subhalo"]
        
        groupFileOffsets = Int64(i) .- groupFileOffsets
        fileNum = maximum( findall(x->x>=0, groupFileOffsets) )
        groupOffset = groupFileOffsets[fileNum]
        
        offsetFile = groupcat.gcPath(basePath, snapNum, fileNum-1)
        prefix = "Offsets/Subhalo_LHaloTree"
        
    end
    
    f = h5open(offsetFile, "r")
    TreeFile = f[prefix*"File"][groupOffset+1]
    TreeIndex = f[prefix*"Index"][groupOffset+1]
    TreeNum = f[prefix*"Num"][groupOffset+1]
    close(f)
    
    return TreeFile, TreeIndex, TreeNum
    
end

function singleNodeFlat(conn, index, data_in, data_out, count, onlyMPB)
    
    data_out[count+1] = data_in[index+1]
    
    count = count + 1
    count = recProgenitorFlat(conn, index, data_in, data_out, count, onlyMPB)
    
    return count
    
end

function recProgenitorFlat(conn, start_index, data_in, data_out, count, onlyMPB)
    
    firstProg = conn["FirstProgenitor"][start_index+1]
    
    if firstProg < 0
        return count
    end
    
    count = singleNodeFlat(conn, firstProg, data_in, data_out, count, onlyMPB)
    
    if onlyMPB == false
        nextProg = conn["NextProgenitor"][firstProg+1]
        
        while nextProg >= 0
            
            count = singleNodeFlat(conn, nextProg, data_in, data_out, count, onlyMPB)
            
            nextProg = conn["NextProgenitor"][nextProg+1]
            
        end
        
    end
    
    firstProg = conn["FirstProgenitor"][firstProg+1]
    
    return count
    
end

function loadTree(basePath, snapNum, i, fields=nothing, onlyMPB=false)
    
    TreeFile, TreeIndex, TreeNum = treeOffsets(basePath, snapNum, i)
    
    if TreeNum == -1
        print("Warning, empty return")
        return nothing
    end
    
    gName = "Tree" * string(TreeNum)
    nRows = nothing
    
    if typeof(fields) == String
        fields = [fields]
    end
    
    fTree = h5open(treePath(basePath, TreeFile), "r")
    
    if isnothing(fields)
		fields = keys(fTree[gName])
	end
    
    for field in fields
        if field in keys(fTree[gName]) == false
            error("Error: Requested field [" * field * "] not in tree.")
        end
    end
    
    connFields = ["FirstProgenitor", "NextProgenitor"]
    conn = Dict()
    
    for field in connFields
        conn[field] = fTree[gName][field][:]
    end
    
    shape = [i for i in size(conn["FirstProgenitor"])]
    dummy = zeros(Int32, Tuple(shape))
    nRows = singleNodeFlat(conn, TreeIndex, dummy, dummy, 0, onlyMPB)
    
    result = Dict()
    
    for field in fields
        if nRows < 1000
            full_data = fTree[gName][field]
        else
            full_data = fTree[gName][field][:]
        end
        
        dtype = eltype(fTree[gName][field])
        shape = [i for i in size(fTree[gName][field])]
        shape[end] = nRows
        
        data = zeros(dtype, Tuple(shape))
        
        count = singleNodeFlat(conn, TreeIndex, full_data, data, 0, onlyMPB)
        
        result[field] = data
        
    end
    
    close(fTree)
    
    if length(fields) == 1
		return result[fields[1]]
	end
    
    return result
    
end

end
