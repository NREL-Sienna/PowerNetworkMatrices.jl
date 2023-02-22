# Build A Matrix
# Build ABA Matrix
# Build AB Trial

function get_ac_branches(sys)
    # Filter out DC Branches here
    return sort!(
        collect(PSY.get_components(PSY.ACBranch, sys));
        by = x ->
            (PSY.get_number(PSY.get_arc(x).from), PSY.get_number(PSY.get_arc(x).to)),
    )
end

function get_buses(sys)
    return sort!(collect(PSY.get_components(PSY.Bus, sys)); by = x -> PSY.get_number(x))
end

function validate_linear_solver(linear_solver::String)
    if linear_solver ∉ SUPPORTED_LINEAR_SOLVERS
        error(
            "Invalid linear solver. Supported linear solvers are: $(SUPPORTED_LINEAR_SOLVERS)",
        )
    end
    return
end
