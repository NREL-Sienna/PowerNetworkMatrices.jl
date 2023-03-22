var documenterSearchIndex = {"docs":
[{"location":"tutorials/intro_page/#SIIP-Examples","page":"Tutorials","title":"SIIP-Examples","text":"","category":"section"},{"location":"tutorials/intro_page/","page":"Tutorials","title":"Tutorials","text":"All the tutorials for the SIIP project are part of a separate repository SIIP-Examples.","category":"page"},{"location":"quick_start_guide/#Quick-Start-Guide","page":"Quick Start Guide","title":"Quick Start Guide","text":"","category":"section"},{"location":"api/internal/#Internal-API","page":"Internal API Reference","title":"Internal API","text":"","category":"section"},{"location":"api/internal/","page":"Internal API Reference","title":"Internal API Reference","text":"Modules = [PowerNetworkMatrices]\nPublic = false","category":"page"},{"location":"api/internal/#PowerNetworkMatrices.get_data-Tuple{PowerNetworkMatrices.PowerNetworkMatrix}","page":"Internal API Reference","title":"PowerNetworkMatrices.get_data","text":"get_data(\n    mat::PowerNetworkMatrices.PowerNetworkMatrix\n) -> Dict{Int64, Array{Float64, N} where N}\n\n\nreturns the raw array data of the PowerNetworkMatrix\n\n\n\n\n\n","category":"method"},{"location":"api/internal/#PowerNetworkMatrices.get_lookup-Tuple{PowerNetworkMatrices.PowerNetworkMatrix}","page":"Internal API Reference","title":"PowerNetworkMatrices.get_lookup","text":"get_lookup(\n    mat::PowerNetworkMatrices.PowerNetworkMatrix\n) -> Tuple{Dict, Dict}\n\n\nreturns the lookup tuple of the `PowerNetworkMatrix`. The first entry corresponds\nto the first dimension and the second entry corresponds to the second dimension. For\ninstance in Ybus the first dimension is buses and second dimension is buses too, and in\nPTDF the first dimension is branches and the second dimension is buses\n\n\n\n\n\n","category":"method"},{"location":"#PowerNetworkMatrices.jl","page":"Welcome Page","title":"PowerNetworkMatrices.jl","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"CurrentModule = PowerNetworkMatrices","category":"page"},{"location":"#Overview","page":"Welcome Page","title":"Overview","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"PowerNetworkMatrices.jl documentation and code are organized according to the needs of different users depending on their skillset and requirements. In broad terms there are three categories:","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"Modeler: Users that want to run a particular analysis or experiment and use PowerNetworkMatrices.jl to develop data sets.\nModel Developer: Users that want to develop custom components and structs in order to exploit PowerNetworkMatrices.jl features to produce custom data sets.\nCode Base Developers: Users that want to add new core functionalities or fix bugs in the core capabilities of PowerNetworkMatrices.jl.","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"PowerNetworkMatrices.jl is an active project under development, and we welcome your feedback, suggestions, and bug reports.","category":"page"},{"location":"#Installation","page":"Welcome Page","title":"Installation","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"The latest stable release of PowerNetworkMatrices can be installed using the Julia package manager with","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"] add PowerNetworkMatrices","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"For the current development version, \"checkout\" this package with","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"] add PowerNetworkMatrices#master","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"PowerNetworkMatrices has been developed as part of the Scalable Integrated Infrastructure Planning (SIIP) initiative at the U.S. Department of Energy's National Renewable Energy Laboratory (NREL).","category":"page"},{"location":"api/public/#Public-API-Reference","page":"Public API Reference","title":"Public API Reference","text":"","category":"section"},{"location":"api/public/","page":"Public API Reference","title":"Public API Reference","text":"Modules = [PowerNetworkMatrices]\nPublic = true","category":"page"},{"location":"api/public/#PowerNetworkMatrices.AdjacencyMatrix","page":"Public API Reference","title":"PowerNetworkMatrices.AdjacencyMatrix","text":"Nodal incidence matrix (Adjacency) is an N x N matrix describing a power system with N buses. It represents the directed connectivity of the buses in a power system.\n\nThe AdjacencyMatrix Struct is indexed using the Bus Numbers, no need for them to be sequential\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.AdjacencyMatrix-Tuple{Any, Vector{PowerSystems.Bus}}","page":"Public API Reference","title":"PowerNetworkMatrices.AdjacencyMatrix","text":"Builds a AdjacencyMatrix from a collection of buses and branches. The return is an N x N AdjacencyMatrix Array indexed with the bus numbers.\n\nKeyword arguments\n\ncheck_connectivity::Bool: Checks connectivity of the network using Goderya's algorithm\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerNetworkMatrices.AdjacencyMatrix-Tuple{PowerSystems.System}","page":"Public API Reference","title":"PowerNetworkMatrices.AdjacencyMatrix","text":"Builds a AdjacencyMatrix from the system. The return is an N x N AdjacencyMatrix Array indexed with the bus numbers.\n\nKeyword arguments\n\ncheck_connectivity::Bool: Checks connectivity of the network using Goderya's algorithm\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerNetworkMatrices.LODF","page":"Public API Reference","title":"PowerNetworkMatrices.LODF","text":"Line Outage Distribution Factors (LODFs) are a sensitivity measure of how a change in a line’s flow affects the flows on other lines in the system.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.LODF-2","page":"Public API Reference","title":"PowerNetworkMatrices.LODF","text":"Builds the LODF matrix from a group of branches and nodes. The return is a LOLDF array indexed with the branch name.\n\nKeyword arguments\n\ndist_slack::Vector{Float64}: Vector of weights to be used as distributed slack bus.   The distributed slack vector has to be the same length as the number of buses\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.LODF-3","page":"Public API Reference","title":"PowerNetworkMatrices.LODF","text":"Builds the LODF matrix from a system. The return is a LOLDF array indexed with the branch name.\n\nKeyword arguments\n\ndist_slack::Vector{Float64}: Vector of weights to be used as distributed slack bus.   The distributed slack vector has to be the same length as the number of buses\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.PTDF","page":"Public API Reference","title":"PowerNetworkMatrices.PTDF","text":"Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.\n\nThe PTDF struct is indexed using the Bus numbers and branch names\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.PTDF-Tuple{Any, Vector{PowerSystems.Bus}}","page":"Public API Reference","title":"PowerNetworkMatrices.PTDF","text":"Builds the PTDF matrix from a group of branches and nodes. The return is a PTDF array indexed with the bus numbers.\n\nKeyword arguments\n\ndist_slack::Vector{Float64}: Vector of weights to be used as distributed slack bus.   The distributed slack vector has to be the same length as the number of buses\nlinear_solver::String: Linear solver to be used. Options are \"Dense\", \"KLU\" and \"MKLPardiso\ntol::Float64: Tolerance to eliminate entries in the PTDF matrix (default eps())\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerNetworkMatrices.PTDF-Tuple{PowerSystems.System}","page":"Public API Reference","title":"PowerNetworkMatrices.PTDF","text":"Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.\n\nKeyword arguments\n\ndist_slack::Vector{Float64}: Vector of weights to be used as distributed slack bus.   The distributed slack vector has to be the same length as the number of buses\nlinear_solver::String: Linear solver to be used. Options are \"Dense\", \"KLU\" and \"MKLPardiso\ntol::Float64: Tolerance to eliminate entries in the PTDF matrix (default eps())\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerNetworkMatrices.VirtualPTDF","page":"Public API Reference","title":"PowerNetworkMatrices.VirtualPTDF","text":"Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.\n\nThe PTDF struct is indexed using the Bus numbers and branch names\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.VirtualPTDF-2","page":"Public API Reference","title":"PowerNetworkMatrices.VirtualPTDF","text":"Builds the PTDF matrix from a group of branches and nodes. The return is a PTDF array indexed with the bus numbers.\n\nKeyword arguments\n\ndist_slack::Vector{Float64}: Vector of weights to be used as distributed slack bus.   The distributed slack vector has to be the same length as the number of buses\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.VirtualPTDF-3","page":"Public API Reference","title":"PowerNetworkMatrices.VirtualPTDF","text":"Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.\n\nKeyword arguments\n\ndist_slack::Vector{Float64}: Vector of weights to be used as distributed slack bus.   The distributed slack vector has to be the same length as the number of buses\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.Ybus","page":"Public API Reference","title":"PowerNetworkMatrices.Ybus","text":"Nodal admittance matrix (Ybus) is an N x N matrix describing a power system with N buses. It represents the nodal admittance of the buses in a power system.\n\nThe Ybus Struct is indexed using the Bus Numbers, no need for them to be sequential\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.Ybus-2","page":"Public API Reference","title":"PowerNetworkMatrices.Ybus","text":"Builds a Ybus from a collection of buses and branches. The return is a Ybus Array indexed with the bus numbers and the branch names.\n\nKeyword arguments\n\ncheck_connectivity::Bool: Checks connectivity of the network using Goderya's algorithm\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PowerNetworkMatrices.Ybus-Tuple{PowerSystems.System}","page":"Public API Reference","title":"PowerNetworkMatrices.Ybus","text":"Builds a Ybus from the system. The return is a Ybus Array indexed with the bus numbers and the branch names.\n\nKeyword arguments\n\ncheck_connectivity::Bool: Checks connectivity of the network\n\n\n\n\n\n","category":"method"},{"location":"api/public/#PowerNetworkMatrices.validate_connectivity-Tuple{PowerSystems.System}","page":"Public API Reference","title":"PowerNetworkMatrices.validate_connectivity","text":"validate_connectivity(sys::PowerSystems.System) -> Bool\n\n\nChecks the network connectivity of the system.\n\nKeyword arguments\n\nconnectivity_method::Function = goderya_connectivity: Specifies the method used as Goderya's algorithm (goderya_connectivity) or depth first search/network traversal (dfs_connectivity)\nNote that the default Goderya method is more efficient, but is resource intensive and may not scale well on large networks.\n\n\n\n\n\n","category":"method"}]
}
