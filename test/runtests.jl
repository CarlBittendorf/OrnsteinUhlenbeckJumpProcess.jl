using OrnsteinUhlenbeckJumpProcess
using Test, Dates
using Aqua

@testset "OrnsteinUhlenbeckJumpProcess.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OrnsteinUhlenbeckJumpProcess)
    end
end
