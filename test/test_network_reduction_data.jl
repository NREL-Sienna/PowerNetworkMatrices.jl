@testset "Test parallel branch naming" begin
    l1 = Line(nothing)
    set_name!(l1, "A33-1")
    l2 = Line(nothing)
    set_name!(l2, "A33-2")
    bp = PNM.BranchesParallel([l1, l2])
    # Use common substring at start of name if possible:
    @test PNM.get_name(bp) == "A33-double_circuit"
    set_name!(l1, "B1")
    set_name!(l2, "C2")
    bp = PNM.BranchesParallel([l1, l2])
    # Otherwise concatonate names:
    @test PNM.get_name(bp) == "B1_C2_double_circuit"
end


