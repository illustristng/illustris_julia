## cartesian.jl

module cartesian

using HDF5

export cartPath
export getNumPixel
export loadSubset



function cartPath(basePath, cartNum, chunkNum=0)
    filePath_list = [
        basePath * "/cartesian_" * lpad(cartNum, 3, "0") * "/cartesian_" * lpad(cartNum, 3, "0") * "." *  lpad(chunkNum, 0, "0")  * ".hdf5"
    ]

    for filePath in filePath_list
        if isfile(filePath)
            return filePath
        end
    end

    error("No cartesian file found!")
end

function getNumPixel(header)
    return header["NumPixels"]
end

function loadSubset(basePath, cartNum, fields=nothing, bbox=nothing)
    
    result = Dict()

    # make sure fields is not a single element
    if isa(fields, AbstractString)
        fields = [fields]
    end

    # load header from first chunk
    
    f = h5open(cartPath(basePath, cartNum), "r")
    header = h5readattr(cartPath(basePath, cartNum), "Header")
    nPix = getNumPixel(header)
    close(f)

    # decide global read size, starting file chunk, and starting file chunk offset
    if bbox !== nothing
        load_all = false
        start_i, start_j, start_k = bbox[1]
        end_i, end_j, end_k = bbox[2]
        @assert start_i >= 0
        @assert start_j >= 0
        @assert start_k >= 0
        @assert end_i < nPix
        @assert end_j < nPix
        @assert end_k < nPix
    else
        load_all = true
        bbox = [[0, 0, 0], [nPix - 1, nPix - 1, nPix - 1]]
    end
    
    numToRead = (bbox[2][1] - bbox[1][1] + 1) * (bbox[2][2] - bbox[1][2] + 1) * (bbox[2][3] - bbox[1][3] + 1)

    if numToRead == 0
        return result
    end

    f = h5open(cartPath(basePath, cartNum, 0), "r")
    # if fields not specified, load everything; otherwise check entry
    if fields === nothing
        fields = keys(f)
        deleteat!(fields, findall(x->x=="Header",fields))
    else
        for field in fields
            # verify existence
            if !haskey(f, field)
                throw(Exception("Cartesian output does not have field [$field]"))
            end
        end
    end

    for field in fields
        # replace local length with global
        shape = [i for i in size(f[field])]
        shape[end] = numToRead

        # allocate within return dict
        dtype = eltype(f[field])
        result[field] = zeros(dtype, Tuple(shape))
    end
    close(f)

    # loop over chunks
    wOffset = 1
    fileOffset = 0
    origNumToRead = numToRead
    fileNum = 0

    while numToRead > 0
        f = h5open(cartPath(basePath, cartNum, fileNum), "r")
        # set local read length for this file chunk, truncate to be within the local size
        numPixelsLocal = size(f[fields[1]])[1]

        if load_all
            pixToReadLocal = trues(numPixelsLocal)
            numToReadLocal = numPixelsLocal
        else
            local_pixels_index = collect(fileOffset + 1:fileOffset + numPixelsLocal)
            local_pixels_i = local_pixels_index .÷ (nPix^2)
            local_pixels_j = (local_pixels_index .- local_pixels_i .* nPix^2) .÷ nPix
            local_pixels_k = local_pixels_index .- local_pixels_i .* nPix^2 .- local_pixels_j .* nPix

            pixToReadLocal = (local_pixels_i .≥ bbox[1][1]) .& (local_pixels_i .≤ bbox[2][1]) .&
                             (local_pixels_j .≥ bbox[1][2]) .& (local_pixels_j .≤ bbox[2][2]) .&
                             (local_pixels_k .≥ bbox[1][3]) .& (local_pixels_k .≤ bbox[2][3])
            numToReadLocal = count(pixToReadLocal)
        end

        # loop over each requested field for this particle type
        for field in fields
            shape = [i for i in size(f[field])]
            dset = read(f, field)
            
            if length(shape) == 1
                result[field][wOffset:wOffset + numToReadLocal - 1] = dset[pixToReadLocal]
            else
                result[field][:,wOffset:wOffset + numToReadLocal - 1] = dset[:,pixToReadLocal]
            end
        end

        wOffset += numToReadLocal
        numToRead -= numToReadLocal

        fileOffset += numPixelsLocal
        fileNum += 1
        close(f)
    end

    # verify we read the correct number
    if origNumToRead != (wOffset-1)
        error("Read [" * string(wOffset-1) * "] particles, but was expecting [" * string(origNumToRead) * "]")
    end

    if length(fields) == 1
		return result[fields[1]]
	end

    return result
end

end
