@kwdef mutable struct ReductionContainer
    radial::Union{Nothing, RadialReduction} = nothing
    degree_two_reduction::Union{Nothing, DegreeTwoReduction} = nothing
    ward_reduction::Union{Nothing, WardReduction} = nothing
end

has_radial_reduction(rb::ReductionContainer) = !isnothing(rb.radial)
has_degree_two_reduction(rb::ReductionContainer) = !isnothing(rb.degree_two_reduction)
has_ward_reduction(rb::ReductionContainer) = !isnothing(rb.ward_reduction)

function validate_reduction_type(
    reduction::T,
    prior_reductions::ReductionContainer,
) where {T <: NetworkReduction}
    prior_reduction_types = [typeof(x) for x in prior_reductions]
    if length(prior_reductions) == 0
        return
    else
        if T ∈ prior_reduction_types
            throw(IS.DataFormatError("$T is applied twice to the same system"))
        end
        if WardReduction ∈ prior_reduction_types
            throw(
                IS.DataFormatError(
                    "$T reduction is applied after Ward reduction. Ward reduction must be applied last.",
                ),
            )
        end
        if T == RadialReduction
            if DegreeTwoReduction ∈ prior_reduction_types
                throw(
                    IS.DataFormatError(
                        "When applying both RadialReduction and DegreeTwoReduction, RadialReduction must be applied first",
                    ),
                )
            end
        end
    end
end
