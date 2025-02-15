const SUPPORTED_LINEAR_SOLVERS = ["KLU", "MKLPardiso", "Dense"]

const KiB = 1024
const MiB = KiB * KiB
const GiB = MiB * KiB
const MAX_CACHE_SIZE_MiB = 100
const ROW_PERSISTENT_CACHE_WARN = 1 * GiB

DEFAULT_LODF_CHUNK_SIZE = 18_000

IS.@scoped_enum(NetworkReductionTypes, RADIAL = 1, WARD = 2, RADIAL_AND_WARD = 3)
