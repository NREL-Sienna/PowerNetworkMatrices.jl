const SUPPORTED_LINEAR_SOLVERS = ["KLU", "MKLPardiso", "Dense"]
const GODERYA_MAX_PERFORMANCE_NODE = 15_000

const KiB = 1024
const MiB = KiB * KiB
const GiB = MiB * KiB
const MAX_CACHE_SIZE_MiB = 100
const ROW_PERSISTENT_CACHE_WARN = 1 * GiB
