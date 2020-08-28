using Posets, Test

@testset "energia no local" begin
    @test energia_no_local([1,2,3], sortperm([1,2,3])) == 0
    @test energia_no_local([1,3,2], sortperm([1,2,3])) == 1
end
