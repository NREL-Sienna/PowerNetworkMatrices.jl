import LinearAlgebra
import SparseArrays
import Random

const _AA_TEST_SEED = 0xABCDEF

function _random_spd(n::Int; density::Float64 = 0.05, scale::Real = 4)
    M = SparseArrays.sprandn(n, n, density)
    return SparseArrays.sparse(M + M' + scale * n * LinearAlgebra.I)
end

@testset "AccelerateWrapper: smoke" begin
    PNM._has_apple_accelerate_backend() || return
    Random.seed!(_AA_TEST_SEED)
    n = 200
    ABA = _random_spd(n)

    cache = PNM.AccelerateWrapper.aa_factorize(ABA)
    @test PNM.AccelerateWrapper.is_factored(cache)

    b = randn(n)
    x = copy(b)
    PNM.AccelerateWrapper.solve!(cache, x)
    @test isapprox(ABA * x, b; atol = 1e-9)

    B = SparseArrays.sprandn(n, 120, 0.03)
    for j in 1:5:120
        B[:, j] .= 0
    end
    SparseArrays.dropzeros!(B)
    out = zeros(n, 120)
    PNM.AccelerateWrapper.solve_sparse!(cache, B, out)
    @test isapprox(out, ABA \ Matrix(B); atol = 1e-9)
end

@testset "AccelerateWrapper: _create_factorization dispatch" begin
    PNM._has_apple_accelerate_backend() || return
    Random.seed!(_AA_TEST_SEED)
    n = 200
    ABA = _random_spd(n)

    K_klu = PNM._create_factorization(PNM.KLUSolver(), ABA)
    K_aa = PNM._create_factorization(PNM.AppleAccelerateSolver(), ABA)
    @test K_klu isa PNM.KLULinSolveCache{Float64}
    @test K_aa isa PNM.AAFactorCache

    b_klu = randn(n)
    b_aa = copy(b_klu)
    PNM._solve_factorization(K_klu, b_klu)
    PNM._solve_factorization(K_aa, b_aa)
    @test isapprox(b_klu, b_aa; atol = 1e-9)
end

@testset "AccelerateWrapper: with_solver resolves to concrete method" begin
    PNM._has_apple_accelerate_backend() || return
    Random.seed!(_AA_TEST_SEED)
    n = 200
    ABA = _random_spd(n)
    K_aa = PNM._create_factorization(PNM.AppleAccelerateSolver(), ABA)

    work = [zeros(n)]
    tmp = [zeros(n)]
    lk = ReentrantLock()
    b = randn(n)
    copyto!(work[1], b)

    result = PNM.with_solver(K_aa, work, tmp, lk) do K, wba, _td
        PNM._solve_factorization(K, wba)
        return copy(wba)
    end
    @test isapprox(ABA * result, b; atol = 1e-9)

    m = which(
        PNM.with_solver,
        (
            typeof(identity),
            PNM.AAFactorCache,
            Vector{Vector{Float64}},
            Vector{Vector{Float64}},
            ReentrantLock,
        ),
    )
    @test Base.unwrap_unionall(m.sig).parameters[3] === PNM.AAFactorCache
end

@testset "AccelerateWrapper: KLU vs AA per-column solve parity" begin
    PNM._has_apple_accelerate_backend() || return
    Random.seed!(_AA_TEST_SEED)
    n = 200
    nbranches = 40
    ABA = _random_spd(n)
    BA = SparseArrays.sprandn(n, nbranches, 0.03)

    K_klu = PNM._create_factorization(PNM.KLUSolver(), ABA)
    K_aa = PNM._create_factorization(PNM.AppleAccelerateSolver(), ABA)

    work_klu = zeros(n)
    work_aa = zeros(n)
    for col in 1:nbranches
        work_klu .= Vector(BA[:, col])
        copyto!(work_aa, work_klu)
        PNM._solve_factorization(K_klu, work_klu)
        PNM._solve_factorization(K_aa, work_aa)
        @test isapprox(work_klu, work_aa; atol = 1e-9)
    end
end
