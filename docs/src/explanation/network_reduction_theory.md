# Network Reduction Theory

This document explains the theory and mathematics behind network reduction techniques in power systems.

## Why Network Reduction?

Power system networks can be very large, with thousands of buses and branches. Network reduction techniques simplify these networks while preserving essential characteristics for analysis.

### Benefits of Reduction:

 1. **Computational Efficiency**: Smaller matrices are faster to compute and invert
 2. **Focus**: Reduces complexity to focus on areas of interest
 3. **Scalability**: Makes large-scale studies tractable
 4. **Clarity**: Simplifies visualization and understanding

### Preservation Goals:

Network reduction aims to preserve:

  - Power flow relationships at retained buses
  - Impedance relationships between key points
  - Stability characteristics (when applicable)
  - Essential topology features

## Radial Branch Reduction

### What is a Radial Branch?

A radial branch is one that connects to a bus that has only one connection to the rest of the network. Think of it as a "dead-end" in the network.

```
Main Network --- Bus A --- Bus B (radial)
                         \
                          Bus C (radial)
```

In this example, buses B and C are radial - they each connect to only one other bus.

### Mathematical Basis

Radial buses can be reduced because power flow to/from them is uniquely determined by the connection bus. The reduction involves:

 1. **Identifying radial buses**: Buses with degree = 1
 2. **Transferring loads**: Moving load/generation to the connection point
 3. **Removing the branch**: Eliminating the radial connection

### Why Radial Reduction Works

For a radial bus $k$ connected to bus $j$:

$$P_k = \frac{V_j V_k}{X_{jk}} \sin(\theta_j - \theta_k)$$

Under DC approximation with $V_j \approx V_k \approx 1$:

$$P_k \approx \frac{\theta_j - \theta_k}{X_{jk}}$$

Since $P_k$ is determined by the load at bus $k$, the angle $\theta_k$ is uniquely determined:

$$\theta_k = \theta_j - P_k X_{jk}$$

The radial bus angle is completely determined by its parent bus, so it can be eliminated by transferring the load.

### Implications of Radial Reduction

**Preserved:**

  - Power flow at non-radial buses
  - Voltage angles at non-radial buses
  - System topology at the core network

**Changed:**

  - Number of buses and branches (reduced)
  - Detailed behavior at radial locations
  - Local voltage profiles

**When to Use:**

  - Radial connections are not of interest for analysis
  - Focusing on transmission backbone
  - Computational efficiency is important

**When to Avoid:**

  - Studying distribution feeders (often radial)
  - Detailed voltage analysis needed
  - Radial buses have critical measurements

## Degree Two Reduction (Kron Reduction)

### What is a Degree Two Bus?

A degree two bus connects exactly two other buses, acting as a "pass-through" point.

```
Bus A --- Bus B (degree 2) --- Bus C
```

Bus B has degree 2 - it connects A and C but has no other connections.

### Mathematical Basis: Kron Reduction

Kron reduction is a systematic method to eliminate buses from the admittance matrix.

Given admittance matrix:

$$\begin{bmatrix} I_r \\ I_e \end{bmatrix} = \begin{bmatrix} Y_{rr} & Y_{re} \\ Y_{er} & Y_{ee} \end{bmatrix} \begin{bmatrix} V_r \\ V_e \end{bmatrix}$$

where subscript $r$ denotes retained buses and $e$ denotes eliminated buses.

If $I_e = 0$ (no injection at eliminated buses):

$$Y_{ee} V_e = -Y_{er} V_r$$
$$V_e = -Y_{ee}^{-1} Y_{er} V_r$$

Substituting back:

$$I_r = (Y_{rr} - Y_{re} Y_{ee}^{-1} Y_{er}) V_r = Y_{reduced} V_r$$

The reduced admittance matrix is:

$$Y_{reduced} = Y_{rr} - Y_{re} Y_{ee}^{-1} Y_{er}$$

### For a Single Degree Two Bus

Consider bus $k$ connecting buses $i$ and $j$:

```
Bus i ----(y_ik)---- Bus k ----(y_kj)---- Bus j
```

The equivalent admittance directly between $i$ and $j$ after eliminating $k$ is:

$$y_{ij}^{new} = y_{ij}^{old} + \frac{y_{ik} \cdot y_{kj}}{y_{ik} + y_{kj} + y_{kk}}$$

If there's no shunt at bus $k$ ($y_{kk} = 0$):

$$y_{ij}^{new} = y_{ij}^{old} + \frac{y_{ik} \cdot y_{kj}}{y_{ik} + y_{kj}}$$

This is analogous to combining series impedances in circuit theory.

### Physical Interpretation

Eliminating a degree two bus:

  - Combines the series impedances of the two connecting branches
  - Creates an equivalent direct connection
  - Preserves the overall impedance between endpoints

### Why Degree Two Reduction Works

The key insight is that a degree two bus with no injection ($P_k = 0$) serves only to pass power through. Its voltage angle is determined by:

$$P_k = \frac{\theta_i - \theta_k}{X_{ik}} + \frac{\theta_j - \theta_k}{X_{kj}} = 0$$

Solving for $\theta_k$:

$$\theta_k = \frac{X_{kj} \theta_i + X_{ik} \theta_j}{X_{ik} + X_{kj}}$$

The angle is a weighted average of its neighbors, so eliminating it and creating a direct equivalent connection preserves the flow relationship.

### Implications of Degree Two Reduction

**Preserved:**

  - Power flows at retained buses
  - Voltage angles at retained buses
  - Overall impedance relationships
  - Equivalence for power flow studies

**Changed:**

  - Number of buses (reduced)
  - Detailed behavior at eliminated bus
  - Branch topology (new equivalent branches created)

**When to Use:**

  - Pass-through buses are not of analytical interest
  - Focusing on injection/load buses
  - Simplifying large transmission systems
  - Reducing computational burden

**When to Avoid:**

  - The degree two bus has measurements or controls
  - Detailed branch flows needed at that location
  - Bus has significant shunt elements
  - Studying protection or relay settings

## Combining Reductions

Multiple reduction techniques can be applied sequentially:

 1. **Order matters**: Applying radial reduction first may expose new degree two buses
 2. **Iterative approach**: After each reduction, new candidates may appear
 3. **Convergence**: Process continues until no more reducible buses exist

### Example Sequence:

```
Original Network (100 buses)
    ↓
Radial Reduction (→ 85 buses, new degree-2 buses exposed)
    ↓
Degree Two Reduction (→ 70 buses, new radial buses exposed)
    ↓
Radial Reduction (→ 65 buses)
    ↓
No more reducible buses
```

## Practical Considerations

### Load and Generation

When eliminating a bus:

  - **With load/generation**: Must be transferred to retained buses
  - **Zero injection**: Can be eliminated cleanly
  - **Load distribution**: May use weighted schemes based on impedances

### Shunt Elements

Buses with shunt admittances (capacitors, reactors):

  - Complicate degree two reduction
  - May need to be distributed to neighbors
  - Can affect the equivalent admittance calculation

### Numerical Stability

  - Very small impedances can cause numerical issues
  - Parallel branches should be combined first
  - Check for near-zero denominators in calculations

### Validation

Always verify reductions:

  - Compare power flows before and after
  - Check that key quantities are preserved
  - Validate on a small test system first

## Applications in PowerNetworkMatrices.jl

### Market Studies

Reducing the network focuses analysis on:

  - Key trading hubs
  - Major transmission corridors
  - Critical generation centers

### Security Analysis

Reduction helps:

  - Speed up contingency screening
  - Focus on critical elements
  - Enable larger study scope

### Planning Studies

Benefits include:

  - Simplifying future scenarios
  - Faster iteration on designs
  - Clearer visualization of key paths

## Limitations

### What Reduction Cannot Do:

  - Preserve detailed voltage profiles everywhere
  - Capture local dynamics at eliminated buses
  - Represent phenomena requiring those buses (e.g., local stability)
  - Maintain exact AC power flow at eliminated locations

### When Full Network is Required:

  - Detailed state estimation
  - Protection coordination
  - Local voltage studies
  - Distributed generation integration at eliminated buses

## Mathematical Guarantees

Under DC power flow assumptions, network reduction:

  - **Preserves**: Impedance matrix (Schur complement)
  - **Preserves**: Power flows at retained buses
  - **Preserves**: Voltage angles at retained buses
  - **Maintains**: Linearity of the system

These guarantees hold exactly for the DC approximation and approximately for AC systems near the linearization point.

## Further Reading

For practical application, see:

  - [How to Apply Network Reduction](@ref) guide
  - [Radial Reduction Tutorial](@ref)
  - [Degree Two Reduction Tutorial](@ref)
