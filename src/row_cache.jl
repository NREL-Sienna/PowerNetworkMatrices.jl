"""
Structure used for saving the rows of the Virtual PTDF and LODF matrix.

# Arguments
- `temp_cache::Dict{Int, Union{Vector{Float64}, SparseArrays.SparseVector{Float64}}}`:
        Dictionay saving the row of the PTDF/LODF matrix
- `persistent_cache_keys::Set{Int}`:
        Set listing the rows saved in `temp_cache`
- `max_cache_size::Int`
        Defines the maximum allowed cache size (rows*row_size)
- `max_num_keys::Int`
        Defines the maximum number of keys saved (rows of the matrix)
"""
struct RowCache{T <: Union{Vector{Float64}, SparseArrays.SparseVector{Float64}}}
    temp_cache::Dict{Int, T}
    persistent_cache_keys::Set{Int}
    max_cache_size::Int
    max_num_keys::Int
end

"""
Structure used for saving the rows of the Virtual PTDF and LODF matrix.

# Arguments
- `temp_cache::Dict{Int, Union{Vector{Float64}, SparseArrays.SparseVector{Float64}}}`:
        Dictionay saving the rows of the PTDF/LODF matrix
- `persistent_cache_keys::Set{Int}`:
        Set listing the rows saved in `temp_cache`
- `max_cache_size::Int`
        Defines the maximum allowed cache size (rows*row_size)
- `max_num_keys::Int`
        Defines the maximum number of keys saved (rows of the matrix)
"""
function RowCache(max_cache_size::Int, persistent_rows::Set{Int}, row_size)
    persistent_data_size = (length(persistent_rows) + 1) * row_size
    if persistent_data_size > max_cache_size
        error(
            "The required cache size for the persisted row is larger than the max cache size. Persistent data size = $(persistent_data_size), max cache size = $(max_cache_size)",
        )
    else
        @debug "required cache for persisted values = $((length(persistent_rows) + 1)*row_size). Max cache specification = $(max_cache_size)"
    end
    max_num_keys = max(length(persistent_rows) + 1, floor(Int, max_cache_size / row_size))
    return RowCache(
        sizehint!(
            Dict{Int, Union{Vector{Float64}, SparseArrays.SparseVector{Float64}}}(),
            max_num_keys,
        ),
        persistent_rows,
        max_cache_size,
        max_num_keys,
    )
end

"""
Check if cache is empty.
"""
function Base.isempty(cache::RowCache)
    return isempty(cache.temp_cache)
end

"""
Erases the cache.
"""
function Base.empty!(cache::RowCache)
    isempty(cache.temp_cache) && return
    if !isempty(cache.persistent_cache_keys)
        @warn("Calling empty! will delete entries for the persistent rows")
    end
    empty!(cache.temp_cache)
    return
end

"""
Checks if `key` is present as a key of the dictionary in `cache`

# Arguments
- `cache::RowCache`:
        cache where data is stored.
- `key::Int`:
        row number (corresponds to the enumerated branch index).
"""
function Base.haskey(cache::RowCache, key::Int)
    return haskey(cache.temp_cache, key)
end

"""
Allocates vector as row of the matrix saved in cache.

# Arguments
- `cache::RowCache`:
        cache where the row vector is going to be saved
- `val::Union{Vector{Float64}, SparseArrays.SparseVector{Float64}}`:
        vector to be saved
- `key::Int`:
        row number (corresponding to the enumerated branch index) related to the input row vector
"""
function Base.setindex!(
    cache::RowCache,
    val::Union{Vector{Float64}, SparseArrays.SparseVector{Float64}},
    key::Int,
)
    check_cache_size!(cache)
    cache.temp_cache[key] = val
    return
end

"""
Gets the row of the stored matrix in cache.

# Arguments
- `cache::RowCache`:
        cache where the row vector is going to be saved
- `key::Int`:
        row number (corresponding to the enumerated branch index) related to the row vector.
"""
function Base.getindex(
    cache::RowCache,
    key::Int,
)::Union{Vector{Float64}, SparseArrays.SparseVector{Float64}}
    return cache.temp_cache[key]
end

"""
Shows the number of rows stored in cache
"""
function Base.length(cache::RowCache)
    return length(cache.temp_cache)
end

"""
Deletes the row of of key `key` from the stored matrix in cache.
"""
function purge_one!(cache::RowCache)
    for k in keys(cache.temp_cache)
        if k ∉ cache.persistent_cache_keys
            delete!(cache.temp_cache, k)
            break
        end
    end
    return
end

"""
Check saved rows in cache and deletes those ones that are not present in "persistent_cache_keys".
"""
function check_cache_size!(cache::RowCache)
    if length(cache.temp_cache) > cache.max_num_keys
        purge_one!(cache)
    end
    return
end
