# How to Apply Network Reduction

This guide shows you how to apply network reduction techniques to simplify your power system model.

## Prerequisites

- `PowerNetworkMatrices.jl` installed
- A power system model with reduction candidates

## Types of Network Reduction

`PowerNetworkMatrices.jl` supports two main reduction techniques:

1. **Radial Branch Reduction** - Removes radial branches
2. **Degree Two Reduction** - Reduces buses with degree two

## Applying Radial Reduction

Radial branches can be reduced to simplify the network while preserving essential characteristics using [`RadialReduction`](@ref):

```julia
using PowerNetworkMatrices
const PNM = PowerNetworkMatrices

# Apply radial reduction
reduced_sys = PNM.RadialReduction(sys)
```

### When to Use Radial Reduction

- System has significant radial branches
- You want to reduce computational complexity
- Radial connections don't affect your analysis

## Applying Degree Two Reduction

Buses with exactly two connections can be eliminated using [`DegreeTwoReduction`](@ref):

```julia
# Apply degree two reduction
reduced_sys = PNM.DegreeTwoReduction(sys)
```

### When to Use Degree Two Reduction

- Many pass-through buses exist
- You're focusing on key interconnection points
- Computational efficiency is important

## Combining Reductions

You can apply multiple reductions sequentially:

```julia
# First apply radial reduction
temp_sys = PNM.RadialReduction(sys)

# Then apply degree two reduction
fully_reduced_sys = PNM.DegreeTwoReduction(temp_sys)
```

## Computing Matrices on Reduced Systems

After reduction, compute network matrices normally using [`PTDF`](@ref):

```julia
# Reduce the system
reduced_sys = PNM.RadialReduction(sys)

# Compute PTDF on reduced system
ptdf_reduced = PNM.PTDF(reduced_sys)
```

## Validating Reductions

Always verify that reductions preserve important system properties:

```julia
# Check system size before reduction
original_bus_count = length(get_components(Bus, sys))

# Apply reduction
reduced_sys = PNM.RadialReduction(sys)

# Check system size after reduction
reduced_bus_count = length(get_components(Bus, reduced_sys))

@info "Reduced from $original_bus_count to $reduced_bus_count buses"
```

## Limitations and Considerations

### Radial Reduction

- Only affects radial (single-connection) branches
- Preserves power flow at non-radial buses
- May not be suitable for all analysis types

### Degree Two Reduction

- Only affects buses with exactly two connections
- Assumes linear flow characteristics
- Check if critical measurement points are preserved

## Related Topics

- [Radial Reduction Tutorial](@ref) - Detailed walkthrough
- [Degree Two Reduction Tutorial](@ref) - Detailed walkthrough
- [Network Reduction Theory](@ref) - Understanding the mathematics
