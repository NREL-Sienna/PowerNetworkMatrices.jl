"""
Serialize the PTDF to an HDF5 file.

# Arguments
- `ptdf::PTDF`: matrix
- `filename::AbstractString`: File to create
- `compress::Bool`: Whether to enabled compression, defaults to true.
- `compression_level::Int`: Compression level to use if compression is enabled.
- `force::Bool`: Whether to overwrite the file if it already exists, defaults to false.
"""
function to_hdf5(
    ptdf::PTDF,
    filename::AbstractString;
    compress = true,
    compression_level = 3,
    force = false,
)
    if isfile(filename) && !force
        error("$filename already exists. Choose a different name or set force=true.")
    end

    lookup1_keys = sort!(collect(keys(ptdf.lookup[1])))
    lookup2_keys = sort!(collect(keys(ptdf.lookup[2])))
    lookup1_values = [ptdf.lookup[1][x] for x in lookup1_keys]
    lookup2_values = [ptdf.lookup[2][x] for x in lookup2_keys]
    ref_bus_positions = collect(ptdf.ref_bus_positions)
    num_rows = length(first(values(ptdf.subnetworks)))
    subnetworks_keys = sort!(collect(keys(ptdf.subnetworks)))
    subnetworks = Matrix{Int}(undef, num_rows, length(subnetworks_keys))

    for (col, key) in enumerate(subnetworks_keys)
        vals = collect(ptdf.subnetworks[key])
        if length(vals) != num_rows
            error("Mismatch in number of rows: num_rows=$num_rows key=$key $vals")
        end
        subnetworks[:, col] = vals
    end

    HDF5.h5open(filename, "w") do file
        _to_hdf5_dataset(file, ptdf.data, compress, compression_level)
        HDF5.attributes(file)["data_type"] = _data_type_as_string(ptdf.data)
        HDF5.attributes(file)["tol"] = ptdf.tol[]
        file["axes1"] = ptdf.axes[1]
        file["axes2"] = ptdf.axes[2]
        file["lookup1_keys"] = lookup1_keys
        file["lookup1_values"] = lookup1_values
        file["lookup2_keys"] = lookup2_keys
        file["lookup2_values"] = lookup2_values
        file["ref_bus_positions"] = ref_bus_positions
        file["subnetworks_keys"] = subnetworks_keys
        file["subnetworks"] = subnetworks
    end

    return
end

"""
Deserialize a PTDF from an HDF5 file.

# Arguments
- `::Type{PTDF}`:
- `filename::AbstractString`: File containing a serialized PTDF.
"""
function from_hdf5(::Type{PTDF}, filename::AbstractString)
    HDF5.h5open(filename, "r") do file
        data = _from_hdf5_dataset(_read_data_type(file), file)
        tol = Base.RefValue(HDF5.read(HDF5.attributes(file)["tol"]))
        axes1 = read(file["axes1"])
        axes2 = read(file["axes2"])
        axes = (axes1, axes2)
        lookup1_keys = read(file["lookup1_keys"])
        lookup1_values = read(file["lookup1_values"])
        lookup2_keys = read(file["lookup2_keys"])
        lookup2_values = read(file["lookup2_values"])
        lookup1 = Dict(k => v for (k, v) in zip(lookup1_keys, lookup1_values))
        lookup2 = Dict(k => v for (k, v) in zip(lookup2_keys, lookup2_values))
        lookup = (lookup1, lookup2)
        ref_bus_positions = Set(read(file["ref_bus_positions"]))
        subnetworks_keys = read(file["subnetworks_keys"])
        subnetworks_matrix = read(file["subnetworks"])
        subnetworks = Dict{Int, Set{Int}}()
        for (col, key) in enumerate(subnetworks_keys)
            subnetworks[key] = Set(subnetworks_matrix[:, col])
        end

        PTDF(
            data,
            axes,
            lookup,
            subnetworks,
            ref_bus_positions,
            tol,
            RadialNetworkReduction(),
        )
    end
end

_data_type_as_string(::Matrix) = "Matrix"
_data_type_as_string(::SparseArrays.SparseMatrixCSC) = "SparseMatrixCSC"

function _read_data_type(file)
    data_type = HDF5.read(HDF5.attributes(file)["data_type"])
    if data_type == "Matrix"
        return Matrix
    elseif data_type == "SparseMatrixCSC"
        return SparseArrays.SparseMatrixCSC
    else
        error("Unsupported type: $(data_type)")
    end
end

function _to_hdf5_dataset(file, data::Matrix, compress, compression_level)
    if compress
        file["data", shuffle = (), deflate = compression_level] = data
    else
        file["data"] = data
    end
end

function _to_hdf5_dataset(
    file,
    data::SparseArrays.SparseMatrixCSC,
    compress,
    compression_level,
)
    for field in (:colptr, :rowval, :nzval)
        key = string(field)
        val = getproperty(data, field)
        if compress
            file[key, shuffle = (), deflate = compression_level] = val
        else
            file[key] = val
        end
    end
    HDF5.attributes(file)["m"] = data.m
    HDF5.attributes(file)["n"] = data.n
end

function _from_hdf5_dataset(::Type{<:Matrix}, file)
    return read(file["data"])
end

function _from_hdf5_dataset(::Type{SparseArrays.SparseMatrixCSC}, file)
    return SparseArrays.SparseMatrixCSC(
        HDF5.read(HDF5.attributes(file)["m"]),
        HDF5.read(HDF5.attributes(file)["n"]),
        read(file["colptr"]),
        read(file["rowval"]),
        read(file["nzval"]),
    )
end
