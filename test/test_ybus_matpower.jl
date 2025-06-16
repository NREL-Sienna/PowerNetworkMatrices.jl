using DelimitedFiles

@testset "Ybus - ACTIVSg10k" begin
    sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg10k_sys")
    matpower_cols =
        readdlm(joinpath(pwd(), "test", "test_data", "ybus_10k_cols.csv"), Int64)
    matpower_rows =
        readdlm(joinpath(pwd(), "test", "test_data", "ybus_10k_rows.csv"), Int64)
    matpower_vals_re = readdlm(joinpath(pwd(), "test", "test_data", "ybus_10k_vals_re.csv"))
    matpower_vals_im = readdlm(joinpath(pwd(), "test", "test_data", "ybus_10k_vals_im.csv"))
    ybus_pnm = Ybus(sys)
    @test nnz(ybus_pnm.data) ==
          length(filter(!iszero, matpower_vals_re .+ im .* matpower_vals_im))
    for (col, row, val_re, val_im) in
        zip(matpower_cols, matpower_rows, matpower_vals_re, matpower_vals_im)
        @test isapprox(ybus_pnm.data[row, col], Complex(val_re, val_im))
    end
end

@testset "Ybus - CATS" begin
    sys = System(joinpath(pwd(), "test", "test_data", "CaliforniaTestSystem.m"))
    matpower_cols =
        readdlm(joinpath(pwd(), "test", "test_data", "ybus_cats_cols.csv"), Int64)
    matpower_rows =
        readdlm(joinpath(pwd(), "test", "test_data", "ybus_cats_rows.csv"), Int64)
    matpower_vals_re =
        readdlm(joinpath(pwd(), "test", "test_data", "ybus_cats_vals_re.csv"))
    matpower_vals_im =
        readdlm(joinpath(pwd(), "test", "test_data", "ybus_cats_vals_im.csv"))
    ybus_pnm = Ybus(sys)
    @test nnz(ybus_pnm.data) ==
          length(filter(!iszero, matpower_vals_re .+ im .* matpower_vals_im))
    for (col, row, val_re, val_im) in
        zip(matpower_cols, matpower_rows, matpower_vals_re, matpower_vals_im)
        @test isapprox(ybus_pnm.data[row, col], Complex(val_re, val_im))
    end
end

    end 
end 