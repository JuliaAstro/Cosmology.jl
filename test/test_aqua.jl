using Aqua

@testset "Aqua tests" begin
    using Cosmology
    using Aqua

    Aqua.test_all(Cosmology)
end
