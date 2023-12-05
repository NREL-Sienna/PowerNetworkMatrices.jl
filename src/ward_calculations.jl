"""
structure and functions for power systems analysis based on Ward's decomposition
"""

"""
Data related to the considered Ward decomposition is here stored.
Given the set of buses under study, the original Ybus can be divided into 
sub matrices belongin to study (S) and external (E) buses.
       |Y_ss  Y_se|
Ybus = |          |
       |Y_es  Y_ee|

# Arguments
- `A::Matrix{ComplexF64}`:
        the product of Y_se * inv(Y_ee) * Y_es
- `B::Matrix{ComplexF64}`:
        the transposed product of Y_se * inv(Y_ee)
- `study_buses::Vector{Int64}`:
        vector of the study buses (all other buses will be aggregated according 
        to the Ward's decomposition)
"""
struct Ward_data
    A::Matrix{ComplexF64}
    B::Matrix{ComplexF64}
    study_buses::Vector{Int64}
end

function Ward_data(
    ybus::Ybus{Tuple{Vector{Int64}, Vector{Int64}}, Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}}},
    study_buses::Vector{Int64}
)
    # get the indeces of the study and external buses
    select_x = [ybus.lookup[1][i] for i in study_buses];
    select_y = [ybus.lookup[2][i] for i in study_buses];
    unselect_x = [ybus.lookup[1][i] for i in keys(ybus.lookup[1]) if i ∉ study_buses];
    unselect_y = [ybus.lookup[2][i] for i in keys(ybus.lookup[1]) if i ∉ study_buses];
    # first reorder the Ybus to get the intermediate sub matrices
    Y_se = ybus.data[select_x, unselect_y];
    Y_es = ybus.data[unselect_x, select_y];
    Y_ee = ybus.data[unselect_x, unselect_y];
    # get the matrices
    A = Y_se*KLU.solve!(KLU.klu(Y_ee), Matrix(Y_es));
    B = KLU.solve!(KLU.klu(Y_ee), Matrix(transpose(Y_se)));
    return Ward_data(
        A,
        B,
        study_buses 
    )
end