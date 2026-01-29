# NREL-Sienna Programming Practices

This document outlines general programming practices, conventions, and guidelines for all NREL-Sienna packages.

## Julia Compatibility

Sienna packages typically support Julia **^1.9** or later.

## Performance Requirements

**Priority:** Critical
**Reference:** [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)

### Anti-Patterns to Avoid

#### Type Instability
Functions must return consistent concrete types.
- ❌ Bad: `f(x) = x > 0 ? 1 : 1.0`
- ✅ Good: `f(x) = x > 0 ? 1.0 : 1.0`
- Check: Use `@code_warntype`

#### Abstract Field Types
Struct fields must have concrete types or be parameterized.
- ❌ Bad: `struct Foo; data::AbstractVector; end`
- ✅ Good: `struct Foo{T<:AbstractVector}; data::T; end`

#### Untyped Containers
- ❌ Bad: `Vector{Any}()`, `Vector{Real}()`
- ✅ Good: `Vector{Float64}()`, `Vector{Int}()`

#### Non-const Globals
- ❌ Bad: `THRESHOLD = 0.5`
- ✅ Good: `const THRESHOLD = 0.5`

#### Unnecessary Allocations
Patterns to follow:
- Use views instead of copies (`@view`, `@views`)
- Pre-allocate arrays instead of `push!` in loops
- Use in-place operations (functions ending with `!`)

#### Captured Variables
Avoid closures that capture variables causing boxing. Solution: pass variables as function arguments instead.

#### Splatting Penalty
Avoid splatting (`...`) in performance-critical code.

#### Abstract Return Types
Avoid returning Union types or abstract types.

### Best Practices

- Use `@inbounds` when bounds are verified
- Use broadcasting (dot syntax) for element-wise operations
- Avoid try-catch in hot paths
- Use function barriers to isolate type instability
- Use sparse matrix operations throughout (`SparseMatrixCSC`)
- Leverage specialized factorization methods for sparse linear solves

**Note:** Apply these guidelines with judgment. Not every function is performance-critical. Focus optimization efforts on hot paths and frequently called code.

## Code Conventions

**Style Guide:** [NREL-Sienna Style Guide](https://nrel-sienna.github.io/InfrastructureSystems.jl/stable/style/)

### Formatter
- **Tool**: JuliaFormatter
- **Command**: `julia -e 'include("scripts/formatter/formatter_code.jl")'`

### Key Rules
- **Constructors**: use `function Foo()` not `Foo() = ...`
- **Asserts**: prefer `InfrastructureSystems.@assert_op` over `@assert`
- **Globals**: UPPER_CASE for constants
- **Exports**: all exports in main module file
- **Comments**: complete sentences, describe why not how
- **Sparse matrices**: use `SparseMatrixCSC` throughout, avoid dense when possible

## Documentation Practices

**Framework:** [Diataxis](https://diataxis.fr/)
**Sienna Guide:** [Documentation Best Practices](https://nrel-sienna.github.io/InfrastructureSystems.jl/stable/docs_best_practices/explanation/)

### Docstring Requirements
- **Scope**: all elements of public interface
- **Include**: function signatures and arguments list
- **Automation**: `DocStringExtensions.TYPEDSIGNATURES`
- **See also**: add links for functions with same name (multiple dispatch)

### API Docs
- **Public**: `docs/src/api/public.md` using `@autodocs` with `Public=true, Private=false`
- **Internals**: `docs/src/api/internals.md`

## Common Tasks

```bash
# Develop locally
julia --project=test -e 'using Pkg; Pkg.develop(path=".")'

# Run tests
julia --project=test test/runtests.jl

# Build documentation
julia --project=docs docs/make.jl

# Format code
julia -e 'include("scripts/formatter/formatter_code.jl")'

# Check formatting
git diff --exit-code

# Instantiate test environment
julia --project=test -e 'using Pkg; Pkg.instantiate()'
```

## Contribution Workflow

- **Branch naming**: `feature/description` or `fix/description` (branches in main repo)
- **Main branch**: `main`

### PR Process
1. Create a feature branch in the main repo
2. Make changes following the style guide
3. Run formatter before committing
4. Ensure tests pass
5. Submit pull request

## General Troubleshooting

### Type Instability
- **Symptom**: Poor performance, many allocations
- **Diagnosis**: Use `@code_warntype` on suspect function
- **Solution**: See Performance Requirements → Anti-Patterns to Avoid

### Formatter Fails
- **Symptom**: Formatter command returns error
- **Solution**: `julia -e 'include("scripts/formatter/formatter_code.jl")'`

### Test Failures
- **Symptom**: Tests fail unexpectedly
- **Solution**: `julia --project=test -e 'using Pkg; Pkg.instantiate()'`

## AI Agent Code Generation Guidelines

### Priorities
- Performance matters - use concrete types in hot paths
- Use sparse matrices (`SparseMatrixCSC`) by default when appropriate
- Apply anti-patterns list with judgment (not exhaustively everywhere)
- Run formatter on all changes
- Add docstrings to public interface elements
- Consider type stability in performance-critical functions

### When Modifying Code
- Read existing code patterns before making changes
- Maintain consistency with existing style
- Prefer failing fast with clear errors over silent failures
- Test with various input scenarios
