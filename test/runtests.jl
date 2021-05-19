using PBWT
using Test

@testset "PBWT.jl" begin
    H =  Array{Int8, 2}([0 1 0 1 0 1; 1 1 0 0 0 1; 1 1 1 1 1 1; 0 1 1 1 1 0; 0 0 0 0 0 0; 1 0 0 0 1 0; 1 1 0 0 0 1; 0 1 0 1 1 0])
    ppa = Array{Int32, 2}([1 1 5 5 5 5; 2 4 6 6 6 2; 3 5 1 1 2 7; 4 8 4 8 7 1; 5 2 8 2 1 6; 6 3 2 7 8 8; 7 6 3 4 4 4; 8 7 7 3 3 3])
    div = Array{Int32, 2}([1 2 3 4 5 6; 1 1 2 2 2 3; 1 1 3 3 3 1; 1 1 1 1 1 5; 1 2 1 2 5 6; 1 1 2 1 1 5; 1 1 1 4 4 4; 1 1 1 2 2 2])
    @test build_prefix_and_divergence_arrays(H) == (ppa, div)
end
