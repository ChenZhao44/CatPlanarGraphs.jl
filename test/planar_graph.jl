@testset "PlanarGraph.jl" begin
    rz = planar_rz()
    @test incident(rz, 1, :face) == [2, 4, 30]
    @test ϕ(rz, incident(rz, 1, :ϕ)[1]) == 1
end
