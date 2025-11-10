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
    subnetworks_keys = sort!(collect(keys(ptdf.subnetwork_axes)))
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
        file["subnetworks_keys"] = subnetworks_keys
        for (ix, axes) in enumerate(values(ptdf.subnetwork_axes))
            file["subnetwork_bus_axis_$ix"] = axes[1]
            file["subnetwork_arc_axis_$ix"] = axes[2]
        end
        # Serialize network reduction data
        _serialize_network_reduction_data(
            file,
            ptdf.network_reduction_data,
            compress,
            compression_level,
        )
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
        axes = (axes1, Tuple.(axes2))
        lookup1_keys = read(file["lookup1_keys"])
        lookup1_values = read(file["lookup1_values"])
        lookup2_keys = read(file["lookup2_keys"])
        lookup2_values = read(file["lookup2_values"])
        lookup1 = Dict(k => v for (k, v) in zip(lookup1_keys, lookup1_values))
        lookup2 = Dict(Tuple(k) => v for (k, v) in zip(lookup2_keys, lookup2_values))
        lookup = (lookup1, lookup2)
        subnetworks_keys = read(file["subnetworks_keys"])
        subnetwork_axes = Dict{Int, Tuple{Vector{Int}, Vector{Tuple{Int, Int}}}}()
        for (ix, key) in enumerate(subnetworks_keys)
            bus_axis = read(file["subnetwork_bus_axis_$ix"])
            arc_axis = read(file["subnetwork_arc_axis_$ix"])
            subnetwork_axes[ix] = (bus_axis, Tuple.(arc_axis))
        end
        # Deserialize network reduction data
        network_reduction_data = _deserialize_network_reduction_data(file)
        PTDF(
            data,
            axes,
            lookup,
            subnetwork_axes,
            tol,
            network_reduction_data,
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

# Helper functions for NetworkReductionData serialization
function _serialize_network_reduction_data(
    file::HDF5.File,
    nrd::NetworkReductionData,
    compress::Bool,
    compression_level::Int,
)
    # Check if network reduction data is empty
    if isempty(nrd)
        # Don't create the group if there's no data
        return
    end

    nrd_group = HDF5.create_group(file, "network_reduction_data")

    # Serialize simple sets and maps
    _serialize_set(nrd_group, "irreducible_buses", nrd.irreducible_buses)
    _serialize_set(nrd_group, "removed_buses", nrd.removed_buses)
    _serialize_set_of_tuples(nrd_group, "removed_arcs", nrd.removed_arcs)

    # Serialize bus maps
    _serialize_bus_reduction_map(nrd_group, nrd.bus_reduction_map)
    _serialize_dict_int_int(
        nrd_group,
        "reverse_bus_search_map",
        nrd.reverse_bus_search_map,
    )

    # Serialize added maps
    _serialize_dict_tuple_complex(nrd_group, "added_branch_map", nrd.added_branch_map)
    _serialize_dict_int_complex(
        nrd_group,
        "added_admittance_map",
        nrd.added_admittance_map,
    )

    # Serialize branch maps (these contain PSY objects, so we store metadata)
    _serialize_direct_branch_map(nrd_group, nrd.direct_branch_map)
    _serialize_parallel_branch_map(nrd_group, nrd.parallel_branch_map)
    _serialize_series_branch_map(nrd_group, nrd.series_branch_map, compress, compression_level)
    _serialize_transformer3W_map(nrd_group, nrd.transformer3W_map)

    # Serialize direct_branch_name_map
    _serialize_dict_string_tuple(
        nrd_group,
        "direct_branch_name_map",
        nrd.direct_branch_name_map,
    )

    # Serialize reduction container
    _serialize_reduction_container(nrd_group, nrd.reductions)

    return
end

function _deserialize_network_reduction_data(file::HDF5.File)
    # Check if network reduction data group exists
    if !haskey(file, "network_reduction_data")
        return NetworkReductionData()
    end

    nrd_group = file["network_reduction_data"]

    # Deserialize simple sets and maps
    irreducible_buses = _deserialize_set(nrd_group, "irreducible_buses")
    removed_buses = _deserialize_set(nrd_group, "removed_buses")
    removed_arcs = _deserialize_set_of_tuples(nrd_group, "removed_arcs")

    # Deserialize bus maps
    bus_reduction_map = _deserialize_bus_reduction_map(nrd_group)
    reverse_bus_search_map = _deserialize_dict_int_int(nrd_group, "reverse_bus_search_map")

    # Deserialize added maps
    added_branch_map = _deserialize_dict_tuple_complex(nrd_group, "added_branch_map")
    added_admittance_map = _deserialize_dict_int_complex(nrd_group, "added_admittance_map")

    # Deserialize branch maps
    direct_branch_map = _deserialize_direct_branch_map(nrd_group)
    parallel_branch_map = _deserialize_parallel_branch_map(nrd_group)
    series_branch_map = _deserialize_series_branch_map(nrd_group)
    transformer3W_map = _deserialize_transformer3W_map(nrd_group)

    # Deserialize direct_branch_name_map
    direct_branch_name_map = _deserialize_dict_string_tuple(nrd_group, "direct_branch_name_map")

    # Deserialize reduction container
    reductions = _deserialize_reduction_container(nrd_group)

    # Create NetworkReductionData with all fields
    # Note: reverse maps and derived fields are intentionally left empty
    # as they would contain PSY objects which cannot be serialized
    return NetworkReductionData(
        irreducible_buses = irreducible_buses,
        bus_reduction_map = bus_reduction_map,
        reverse_bus_search_map = reverse_bus_search_map,
        direct_branch_map = direct_branch_map,
        reverse_direct_branch_map = Dict{PSY.ACTransmission, Tuple{Int, Int}}(),
        parallel_branch_map = parallel_branch_map,
        reverse_parallel_branch_map = Dict{PSY.ACTransmission, Tuple{Int, Int}}(),
        series_branch_map = series_branch_map,
        reverse_series_branch_map = Dict{PSY.ACTransmission, Tuple{Int, Int}}(),
        transformer3W_map = transformer3W_map,
        reverse_transformer3W_map = Dict{ThreeWindingTransformerWinding, Tuple{Int, Int}}(),
        removed_buses = removed_buses,
        removed_arcs = removed_arcs,
        added_admittance_map = added_admittance_map,
        added_branch_map = added_branch_map,
        all_branch_maps_by_type = Dict{String, Any}(),
        reductions = reductions,
        name_to_arc_map = Dict{
            Type,
            DataStructures.SortedDict{String, Tuple{Tuple{Int, Int}, String}},
        }(),
        filters_applied = Dict{Type, Function}(),
        direct_branch_name_map = direct_branch_name_map,
    )
end

# Serialization helpers for basic types
function _serialize_set(group::HDF5.Group, name::String, s::Set{Int})
    if isempty(s)
        HDF5.attributes(group)[name * "_empty"] = true
    else
        group[name] = collect(s)
    end
end

function _deserialize_set(group::HDF5.Group, name::String)
    if get(HDF5.attributes(group), name * "_empty", false)
        return Set{Int}()
    end
    if !haskey(group, name)
        return Set{Int}()
    end
    return Set{Int}(read(group[name]))
end

function _serialize_set_of_tuples(
    group::HDF5.Group,
    name::String,
    s::Set{Tuple{Int, Int}},
)
    if isempty(s)
        HDF5.attributes(group)[name * "_empty"] = true
    else
        tuples_array = reduce(hcat, [[t[1], t[2]] for t in s])
        group[name] = tuples_array
    end
end

function _deserialize_set_of_tuples(group::HDF5.Group, name::String)
    if get(HDF5.attributes(group), name * "_empty", false)
        return Set{Tuple{Int, Int}}()
    end
    if !haskey(group, name)
        return Set{Tuple{Int, Int}}()
    end
    tuples_array = read(group[name])
    return Set{Tuple{Int, Int}}([(tuples_array[1, i], tuples_array[2, i]) for i in 1:size(tuples_array, 2)])
end

function _serialize_bus_reduction_map(
    group::HDF5.Group,
    map::Dict{Int, Set{Int}},
)
    if isempty(map)
        HDF5.attributes(group)["bus_reduction_map_empty"] = true
        return
    end

    # Store as JSON for complex nested structure
    json_str = JSON3.write(Dict(string(k) => collect(v) for (k, v) in map))
    HDF5.attributes(group)["bus_reduction_map"] = json_str
end

function _deserialize_bus_reduction_map(group::HDF5.Group)
    if get(HDF5.attributes(group), "bus_reduction_map_empty", false)
        return Dict{Int, Set{Int}}()
    end
    if !haskey(HDF5.attributes(group), "bus_reduction_map")
        return Dict{Int, Set{Int}}()
    end

    json_str = read(HDF5.attributes(group)["bus_reduction_map"])
    json_dict = JSON3.read(json_str)
    return Dict{Int, Set{Int}}(
        parse(Int, k) => Set{Int}(collect(v)) for (k, v) in json_dict
    )
end

function _serialize_dict_int_int(
    group::HDF5.Group,
    name::String,
    dict::Dict{Int, Int},
)
    if isempty(dict)
        HDF5.attributes(group)[name * "_empty"] = true
        return
    end

    keys_arr = collect(keys(dict))
    vals_arr = [dict[k] for k in keys_arr]
    group[name * "_keys"] = keys_arr
    group[name * "_values"] = vals_arr
end

function _deserialize_dict_int_int(group::HDF5.Group, name::String)
    if get(HDF5.attributes(group), name * "_empty", false)
        return Dict{Int, Int}()
    end
    if !haskey(group, name * "_keys")
        return Dict{Int, Int}()
    end

    keys_arr = read(group[name * "_keys"])
    vals_arr = read(group[name * "_values"])
    return Dict{Int, Int}(k => v for (k, v) in zip(keys_arr, vals_arr))
end

function _serialize_dict_tuple_complex(
    group::HDF5.Group,
    name::String,
    dict::Dict{Tuple{Int, Int}, Complex{Float32}},
)
    if isempty(dict)
        HDF5.attributes(group)[name * "_empty"] = true
        return
    end

    keys_arr = reduce(hcat, [[k[1], k[2]] for k in keys(dict)])
    vals_real = [real(dict[Tuple(keys_arr[:, i])]) for i in 1:size(keys_arr, 2)]
    vals_imag = [imag(dict[Tuple(keys_arr[:, i])]) for i in 1:size(keys_arr, 2)]

    group[name * "_keys"] = keys_arr
    group[name * "_values_real"] = vals_real
    group[name * "_values_imag"] = vals_imag
end

function _deserialize_dict_tuple_complex(group::HDF5.Group, name::String)
    if get(HDF5.attributes(group), name * "_empty", false)
        return Dict{Tuple{Int, Int}, Complex{Float32}}()
    end
    if !haskey(group, name * "_keys")
        return Dict{Tuple{Int, Int}, Complex{Float32}}()
    end

    keys_arr = read(group[name * "_keys"])
    vals_real = read(group[name * "_values_real"])
    vals_imag = read(group[name * "_values_imag"])

    return Dict{Tuple{Int, Int}, Complex{Float32}}(
        (keys_arr[1, i], keys_arr[2, i]) => Complex{Float32}(vals_real[i], vals_imag[i])
        for i in 1:size(keys_arr, 2)
    )
end

function _serialize_dict_int_complex(
    group::HDF5.Group,
    name::String,
    dict::Dict{Int, Complex{Float32}},
)
    if isempty(dict)
        HDF5.attributes(group)[name * "_empty"] = true
        return
    end

    keys_arr = collect(keys(dict))
    vals_real = [real(dict[k]) for k in keys_arr]
    vals_imag = [imag(dict[k]) for k in keys_arr]

    group[name * "_keys"] = keys_arr
    group[name * "_values_real"] = vals_real
    group[name * "_values_imag"] = vals_imag
end

function _deserialize_dict_int_complex(group::HDF5.Group, name::String)
    if get(HDF5.attributes(group), name * "_empty", false)
        return Dict{Int, Complex{Float32}}()
    end
    if !haskey(group, name * "_keys")
        return Dict{Int, Complex{Float32}}()
    end

    keys_arr = read(group[name * "_keys"])
    vals_real = read(group[name * "_values_real"])
    vals_imag = read(group[name * "_values_imag"])

    return Dict{Int, Complex{Float32}}(
        k => Complex{Float32}(vals_real[i], vals_imag[i])
        for (i, k) in enumerate(keys_arr)
    )
end

function _serialize_dict_string_tuple(
    group::HDF5.Group,
    name::String,
    dict::Dict{String, Tuple{Int, Int}},
)
    if isempty(dict)
        HDF5.attributes(group)[name * "_empty"] = true
        return
    end

    # Use JSON for string keys
    json_str = JSON3.write(Dict(k => [v[1], v[2]] for (k, v) in dict))
    HDF5.attributes(group)[name] = json_str
end

function _deserialize_dict_string_tuple(group::HDF5.Group, name::String)
    if get(HDF5.attributes(group), name * "_empty", false)
        return Dict{String, Tuple{Int, Int}}()
    end
    if !haskey(HDF5.attributes(group), name)
        return Dict{String, Tuple{Int, Int}}()
    end

    json_str = read(HDF5.attributes(group)[name])
    json_dict = JSON3.read(json_str)
    return Dict{String, Tuple{Int, Int}}(
        String(k) => (Int(v[1]), Int(v[2])) for (k, v) in json_dict
    )
end

# Serialization for branch maps containing PSY objects
# We serialize metadata about the PSY objects, not the objects themselves

# Since we cannot serialize actual PSY objects without the full PowerSystems library,
# we store the arc mappings and metadata only. This preserves the network topology
# and allows the maps to be useful for downstream applications.

function _serialize_direct_branch_map(
    group::HDF5.Group,
    map::Dict{Tuple{Int, Int}, PSY.ACTransmission},
)
    if isempty(map)
        HDF5.attributes(group)["direct_branch_map_empty"] = true
        return
    end

    # Store arc tuples and branch metadata
    arcs = collect(keys(map))
    branch_names = [PSY.get_name(map[arc]) for arc in arcs]
    branch_types = [string(typeof(map[arc])) for arc in arcs]

    arcs_array = reduce(hcat, [[arc[1], arc[2]] for arc in arcs])
    group["direct_branch_map_arcs"] = arcs_array

    # Store metadata as JSON
    metadata = Dict("names" => branch_names, "types" => branch_types)
    HDF5.attributes(group)["direct_branch_map_metadata"] = JSON3.write(metadata)
end

function _deserialize_direct_branch_map(group::HDF5.Group)
    if get(HDF5.attributes(group), "direct_branch_map_empty", false)
        return Dict{Tuple{Int, Int}, PSY.ACTransmission}()
    end
    if !haskey(group, "direct_branch_map_arcs")
        return Dict{Tuple{Int, Int}, PSY.ACTransmission}()
    end

    # For now, return empty dict as we cannot reconstruct PSY objects without the System
    # The arc information is preserved in other serialized fields
    # Users who need the full PSY objects should keep the original System
    return Dict{Tuple{Int, Int}, PSY.ACTransmission}()
end

function _serialize_parallel_branch_map(
    group::HDF5.Group,
    map::Dict{Tuple{Int, Int}, BranchesParallel},
)
    if isempty(map)
        HDF5.attributes(group)["parallel_branch_map_empty"] = true
        return
    end

    # Store arc tuples and metadata for each BranchesParallel
    arcs = collect(keys(map))
    arcs_array = reduce(hcat, [[arc[1], arc[2]] for arc in arcs])
    group["parallel_branch_map_arcs"] = arcs_array

    # For each parallel branch set, store metadata and equivalent_ybus
    metadata_list = []
    for (i, arc) in enumerate(arcs)
        bp = map[arc]
        branch_names = [PSY.get_name(br) for br in bp.branches]
        branch_types = [string(typeof(br)) for br in bp.branches]

        meta = Dict("arc_index" => i, "names" => branch_names, "types" => branch_types)

        # Serialize equivalent_ybus if it exists
        if !isnothing(bp.equivalent_ybus)
            ybus_flat = vec(bp.equivalent_ybus)
            meta["ybus_real"] = real.(ybus_flat)
            meta["ybus_imag"] = imag.(ybus_flat)
        end

        push!(metadata_list, meta)
    end

    HDF5.attributes(group)["parallel_branch_map_metadata"] =
        JSON3.write(metadata_list)
end

function _deserialize_parallel_branch_map(group::HDF5.Group)
    if get(HDF5.attributes(group), "parallel_branch_map_empty", false)
        return Dict{Tuple{Int, Int}, BranchesParallel}()
    end
    if !haskey(group, "parallel_branch_map_arcs")
        return Dict{Tuple{Int, Int}, BranchesParallel}()
    end

    # Return empty dict as we cannot reconstruct PSY objects
    return Dict{Tuple{Int, Int}, BranchesParallel}()
end

function _serialize_series_branch_map(
    group::HDF5.Group,
    map::Dict{Tuple{Int, Int}, BranchesSeries},
    compress::Bool,
    compression_level::Int,
)
    if isempty(map)
        HDF5.attributes(group)["series_branch_map_empty"] = true
        return
    end

    # Store arc tuples
    arcs = collect(keys(map))
    arcs_array = reduce(hcat, [[arc[1], arc[2]] for arc in arcs])
    group["series_branch_map_arcs"] = arcs_array

    # For each series branch set, store metadata
    metadata_list = []
    for (i, arc) in enumerate(arcs)
        bs = map[arc]

        # Collect branch information
        branch_info = []
        for (branch_type, branch_list) in bs.branches
            for br in branch_list
                push!(branch_info, Dict("name" => PSY.get_name(br), "type" => string(branch_type)))
            end
        end

        meta = Dict(
            "arc_index" => i,
            "branches" => branch_info,
            "needs_insertion_order" => bs.needs_insertion_order,
            "insertion_order" => bs.insertion_order,
            "segment_orientations" => [string(s) for s in bs.segment_orientations],
        )

        # Serialize equivalent_ybus if it exists
        if !isnothing(bs.equivalent_ybus)
            ybus_flat = vec(bs.equivalent_ybus)
            meta["ybus_real"] = real.(ybus_flat)
            meta["ybus_imag"] = imag.(ybus_flat)
        end

        push!(metadata_list, meta)
    end

    HDF5.attributes(group)["series_branch_map_metadata"] = JSON3.write(metadata_list)
end

function _deserialize_series_branch_map(group::HDF5.Group)
    if get(HDF5.attributes(group), "series_branch_map_empty", false)
        return Dict{Tuple{Int, Int}, BranchesSeries}()
    end
    if !haskey(group, "series_branch_map_arcs")
        return Dict{Tuple{Int, Int}, BranchesSeries}()
    end

    # Return empty dict as we cannot reconstruct PSY objects
    return Dict{Tuple{Int, Int}, BranchesSeries}()
end

function _serialize_transformer3W_map(
    group::HDF5.Group,
    map::Dict{Tuple{Int, Int}, ThreeWindingTransformerWinding},
)
    if isempty(map)
        HDF5.attributes(group)["transformer3W_map_empty"] = true
        return
    end

    # Store arc tuples and transformer metadata
    arcs = collect(keys(map))
    arcs_array = reduce(hcat, [[arc[1], arc[2]] for arc in arcs])
    group["transformer3W_map_arcs"] = arcs_array

    # Store metadata
    metadata_list = []
    for (i, arc) in enumerate(arcs)
        t3w = map[arc]
        meta = Dict(
            "arc_index" => i,
            "transformer_name" => PSY.get_name(t3w.transformer),
            "transformer_type" => string(typeof(t3w.transformer)),
            "winding_number" => t3w.winding_number,
        )
        push!(metadata_list, meta)
    end

    HDF5.attributes(group)["transformer3W_map_metadata"] = JSON3.write(metadata_list)
end

function _deserialize_transformer3W_map(group::HDF5.Group)
    if get(HDF5.attributes(group), "transformer3W_map_empty", false)
        return Dict{Tuple{Int, Int}, ThreeWindingTransformerWinding}()
    end
    if !haskey(group, "transformer3W_map_arcs")
        return Dict{Tuple{Int, Int}, ThreeWindingTransformerWinding}()
    end

    # Return empty dict as we cannot reconstruct PSY objects
    return Dict{Tuple{Int, Int}, ThreeWindingTransformerWinding}()
end

function _serialize_reduction_container(group::HDF5.Group, rc::ReductionContainer)
    # Serialize each reduction type
    reduction_data = Dict{String, Any}()

    if !isnothing(rc.radial_reduction)
        reduction_data["radial_reduction"] =
            Dict("irreducible_buses" => rc.radial_reduction.irreducible_buses)
    end

    if !isnothing(rc.degree_two_reduction)
        reduction_data["degree_two_reduction"] = Dict(
            "irreducible_buses" => rc.degree_two_reduction.irreducible_buses,
            "reduce_reactive_power_injectors" =>
                rc.degree_two_reduction.reduce_reactive_power_injectors,
        )
    end

    if !isnothing(rc.ward_reduction)
        reduction_data["ward_reduction"] =
            Dict("study_buses" => rc.ward_reduction.study_buses)
    end

    HDF5.attributes(group)["reduction_container"] = JSON3.write(reduction_data)
end

function _deserialize_reduction_container(group::HDF5.Group)
    if !haskey(HDF5.attributes(group), "reduction_container")
        return ReductionContainer()
    end

    json_str = read(HDF5.attributes(group)["reduction_container"])
    reduction_data = JSON3.read(json_str)

    rc = ReductionContainer()

    if haskey(reduction_data, "radial_reduction")
        rd = reduction_data["radial_reduction"]
        rc.radial_reduction = RadialReduction(
            irreducible_buses = collect(Int, rd["irreducible_buses"]),
        )
    end

    if haskey(reduction_data, "degree_two_reduction")
        d2d = reduction_data["degree_two_reduction"]
        rc.degree_two_reduction = DegreeTwoReduction(
            irreducible_buses = collect(Int, d2d["irreducible_buses"]),
            reduce_reactive_power_injectors = Bool(d2d["reduce_reactive_power_injectors"]),
        )
    end

    if haskey(reduction_data, "ward_reduction")
        wd = reduction_data["ward_reduction"]
        rc.ward_reduction = WardReduction(collect(Int, wd["study_buses"]))
    end

    return rc
end
