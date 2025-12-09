using OrnsteinUhlenbeckJumpProcess
using Test
using Aqua

@testset "OrnsteinUhlenbeckJumpProcess.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(OrnsteinUhlenbeckJumpProcess)
    end
    # Write your tests here.
end
