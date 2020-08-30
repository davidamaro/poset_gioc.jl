using Posets, Test

@testset "energia no local" begin
    @test energia_no_local([1,2,3], sortperm([1,2,3])) == 0
    @test energia_no_local([1,3,2], sortperm([1,2,3])) == 1

    mat = zeros(Int,(3,3))
    mat[1,2] = 1
    mat[1,3] = 1
    @test energia_no_local([1,2,3], mat) ==  0
    @test energia_no_local([3,1,2], mat) ==  1
end

@testset "matriz" begin
    mat = zeros(Int,(3,3))
    mat[1,2] = 1
    mat[1,3] = 1

    @test crear_matriz(sortperm.([[1,3,2], [1,2,3]])) ==  mat

    mat = zeros(Int,(3,3))
    mat[1,3] = 1
    mat[1,2] = 1
    mat[3,2] = 1
    @test crear_matriz(sortperm([1,3,2])) ==  mat
end

@testset "reduccion transitiva" begin
    mat = zeros(Int,(3,3))
    mat[1,2] = 1
    mat[2,3] = 1

    @test matriz_rutas(mat) == [0 1 1; 0 0 1; 0 0 0]
    mat = zeros(Int,(3,3))
    mat[1,2] = 1
    mat[2,3] = 1
    mat = matriz_rutas(mat)

    @test reduccion_transitiva(mat) == [0 1 0; 0 0 1; 0 0 0]
end
