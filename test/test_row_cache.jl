@testset "RowCache functions" begin

    # dummy rows
    key1 = 1
    row1 = ones((24,))

    key5 = 5
    row5 = ones((24,))

    # define cache struct
    cache = PNM.RowCache(
        5 * length(row) * sizeof(Float64),
        Set([1, 2, 3]),
        length(row) * sizeof(Float64),
    )

    # test: Base.isempty
    @test isempty(cache) == true

    # test: Base.haskey and Base.getindex
    cache[key1] = row1
    @test cache[key1] == row1
    @test haskey(cache, 1) == true

    # test: Base.length (number of rows stored)
    @test length(cache) == 1

    # test: purge_one!
    PNM.purge_one!(cache)
    @test length(cache) == 1
    @test haskey(cache, 5) == false

    # test: Base.setindex, check_cache_size!
    for i in 1:5
        cache[i] = row1
    end
    cache[6] = row1
    @test length(cache) == 5

    # test: Base.empty!
    empty!(cache)
    @test length(cache) == 0
end
