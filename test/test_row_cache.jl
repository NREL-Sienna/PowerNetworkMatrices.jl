@testset "Test RowCache functions" begin

    # dummy rows
    key = 1
    row = ones((24,))

    # define cache struct
    cache = PNM.RowCache(
        5 * length(row) * sizeof(Float64),
        Set([1, 2, 3]),
        length(row) * sizeof(Float64),
    )

    # test: Base.isempty
    @test isempty(cache) == true

    # test: Base.haskey and Base.getindex
    cache[key] = row
    @test cache[key] == row
    @test haskey(cache, 1) == true

    # test: Base.length (number of rows stored)
    @test length(cache) == 1

    # test: purge_one!
    PNM.purge_one!(cache)
    @test length(cache) == 1
    @test haskey(cache, 5) == false

    # test: Base.setindex, check_cache_size!
    for i in 1:5
        cache[i] = row
    end
    cache[6] = row
    @test length(cache) == 5

    # test: Base.empty!
    empty!(cache)
    @test length(cache) == 0
end