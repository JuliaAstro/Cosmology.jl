using Cosmology: H, cosmology, hubble_dist0, hubble_dist, hubble_time0, hubble_time
using DynamicQuantities: @u_str
using Test

@testset "Utilities" begin
    c = cosmology(h = 0.7)
    @test hubble_time(c, 0) ≈ hubble_time0(c)
    @test hubble_dist(c, 0) ≈ hubble_dist0(c)
    @test H(c, 0) ≈ 70u"km/s/Constants.Mpc"
end
