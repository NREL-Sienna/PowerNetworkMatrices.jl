precompile = @timed using PowerNetworkMatrices

open("precompile_time.txt", "a") do io
    write(io, "| $(ARGS[1]) | $(precompile.time) |\n")
end

using PowerSystems
using PowerSystemCaseBuilder
using Logging

configure_logging(; console_level = Logging.Error)

sys = build_system(MatpowerTestSystems, "matpower_ACTIVSg10k_sys")

try
    _, time_build_ptdf1, _, _ = @timed PTDF(sys)
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build PTDF First | $(time_build_ptdf1) |\n")
    end
    _, time_build_ptdf2, _, _ = @timed PTDF(sys)
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build PTDF Second | $(time_build_ptdf2) |\n")
    end
catch e
    @error exception = (e, catch_backtrace())
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])- Build PTDF | FAILED TO TEST |\n")
    end
end

try
    _, time_build_ybus1, _, _ = @timed Ybus(sys)
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build Ybus First | $(time_build_ybus1) |\n")
    end
    _, time_build_ybus2, _, _ = @timed Ybus(sys)
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build Ybus Second | $(time_build_ybus2) |\n")
    end
catch e
    @error exception = (e, catch_backtrace())
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])- Build Ybus | FAILED TO TEST |\n")
    end
end

try
    _, time_build_LODF1, _, _ = @timed LODF(sys)
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build LODF First | $(time_build_LODF1) |\n")
    end
    _, time_build_LODF2, _, _ = @timed LODF(sys)
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build LODF Second | $(time_build_LODF2) |\n")
    end
catch e
    @error exception = (e, catch_backtrace())
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])- Build LODF | FAILED TO TEST |\n")
    end
end
