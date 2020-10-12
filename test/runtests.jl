using Posets, Test
import LightGraphs: DiGraph
import SparseArrays: spzeros
import StatsBase: countmap
import LinearAlgebra: norm

@testset "energia local" begin
    @test energia_local([1,2,3], lista_posets_3[2]) == 0
    @test energia_local([1,3,2], lista_posets_3[2]) == 0
    @test energia_local([2,1,3], lista_posets_3[2]) == 1
    @test energia_local([2,3,1], lista_posets_3[2]) == 0
    @test energia_local([3,1,2], lista_posets_3[2]) == 0
    @test energia_local([3,2,1], lista_posets_3[2]) == 1
end

@testset "energia no local" begin
    @test energia_no_local([1,2,3], sortperm([1,2,3])) == 0
    @test energia_no_local([1,3,2], sortperm([1,2,3])) == 1

    mat = zeros(Int,(3,3))
    mat[1,2] = 1
    mat[1,3] = 1
    @test energia_no_local([1,2,3], mat) ==  0
    @test energia_no_local([3,1,2], mat) ==  1
end

@testset "matriz asociada a una lista de rankings" begin
  @test matriz_interseccion_rankings([[1,2,3], [1,3,2]]) == [ 0 1 1; 0 0 0; 0 0 0 ]
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

@testset "densidad de un poset" begin
  @test densidad_exacta(lista_posets_3[1],2,energia_no_local) == [6,0,0,0,0,0]
  @test densidad_exacta(lista_posets_3[2],2,energia_no_local) == [3,3,0,0,0,0]
  @test densidad_exacta(lista_posets_3[3],2,energia_no_local) == [2,2,2,0,0,0]
  @test densidad_exacta(lista_posets_3[4],2,energia_no_local) == [1,2,2,1,0,0]
end

@testset "union de posets" begin
  @test matriz_union_rankings([[1,2,3], [1,3,2]], binario = false) == [0 2 2; 0 0 1; 0 1 0]
  @test matriz_union_rankings([[1,2,3], [1,3,2]]) == [0 1 1; 0 0 1; 0 1 0]
end

@testset "aristas resaltadas" begin
  ejemplo_1 = zeros(Int, (3,3))
  ejemplo_1[1,2]= 1
  ejemplo_1[1,3] = 1
  ejemplo_2 = zeros(Int, (3,3))
  ejemplo_2[1,2]= 1
  ejemplo_3 = zeros(Int, (3,3))
  ejemplo_3[2,3]= 1

  @test resaltados(ejemplo_1|>DiGraph, ejemplo_2|>DiGraph) == [2]
  @test resaltados(ejemplo_1|>DiGraph, ejemplo_3|>DiGraph) == [1]
end

@testset "probando aciclicas" begin
    ciclica = spzeros(Int64,6,6)
    ciclica[1,2] = 1
    ciclica[2,3] = 1
    ciclica[2,4] = 1
    ciclica[4,5] = 1
    ciclica[4,6] = 1
    ciclica[5,6] = 1
    ciclica[6,3] = 1
    ciclica[6,4] = 1
    @test !isacyclic(ciclica)
    aciclica = spzeros(Int64,6,6)
    aciclica[1,2] = 1
    aciclica[2,3] = 1
    aciclica[2,4] = 1
    aciclica[4,5] = 1
    aciclica[4,6] = 1
    aciclica[5,6] = 1
    aciclica[6,3] = 1
    @test isacyclic(aciclica)
    chava_fea_1 = spzeros(Int64,5,5)
    chava_fea_1[1,5] = 1
    chava_fea_1[4,1] = 1
    chava_fea_1[4,2] = 1
    chava_fea_1[5,4] = 1
    chava_fea_1[5,3] = 1
    @test !isacyclic(chava_fea_1)
    chava_guapa_1 = spzeros(Int64,5,5)
    chava_guapa_1[1,5] = 1
    chava_guapa_1[4,2] = 1
    chava_guapa_1[5,4] = 1
    chava_guapa_1[5,3] = 1
    @test isacyclic(chava_guapa_1)
    chava_guapa_2 = spzeros(Int64,5,5)
    chava_guapa_2[1,2] = 1
    chava_guapa_2[2,3] = 1
    chava_guapa_2[3,4] = 1
    chava_guapa_2[4,5] = 1
    @test isacyclic(chava_guapa_2)
    chava_fea_2 = spzeros(Int64,5,5)
    chava_fea_2[1,2] = 1
    chava_fea_2[2,3] = 1
    chava_fea_2[3,4] = 1
    chava_fea_2[4,5] = 1
    chava_fea_2[5,3] = 1
    @test !isacyclic(chava_fea_2)
    chava_rara_1 = spzeros(Int64,5,5)
    chava_rara_1[1,2] = 1
    chava_rara_1[1,3] = 1
    chava_rara_1[4,5] = 1
    @test isacyclic(chava_rara_1)
    chava_rara_2 = spzeros(Int64,5,5)
    chava_rara_2[1,2] = 1
    chava_rara_2[2,3] = 1
    chava_rara_2[3,1] = 1
    chava_rara_2[4,5] = 1
    @test !isacyclic(chava_rara_2)
end

@testset "total variation 3 nodos" begin
    function boba()
        lista_posetsrandom = Array{Int64,2}[]
        for i in 1:10000
            #push!(lista_posetsrandom, caminata_poset(4,16) |> Array |> reduccion_transitiva)
            push!(lista_posetsrandom, caminata_poset(3,9))
        end

        #listaaumentada = map(x -> matriz_rutas(Array(x)) |> reduccion_transitiva, lista_posetsrandom);
        map(x -> matriz_rutas(Array(x)) |> reduccion_transitiva, lista_posetsrandom)
    end
    listaaumentada = boba()
    valores = countmap(listaaumentada) |> values |> collect;
    @test sum(abs.((valores./sum(valores)).-(1/19)))/2 < 0.1
end

@testset "total variation 4 nodos" begin
    function boba()
        lista_posetsrandom = Array{Int64,2}[]
        for i in 1:100000
            #push!(lista_posetsrandom, caminata_poset(4,16) |> Array |> reduccion_transitiva)
            push!(lista_posetsrandom, caminata_poset(4,4^2))
        end

        #listaaumentada = map(x -> matriz_rutas(Array(x)) |> reduccion_transitiva, lista_posetsrandom);
        map(x -> matriz_rutas(Array(x)) |> reduccion_transitiva, lista_posetsrandom)
    end
    listaaumentada = boba()
    valores = countmap(listaaumentada) |> values |> collect;
    @test sum(abs.((valores./sum(valores)).-(1/219)))/2 < 0.1
end

@testset "varianza y promedio de posicion" begin
    @test posicionpromedio([[1,2,3], [1,3,2]], 1) ≈ 1.0
    @test posicionvarianza([[1,2,3], [1,3,2]], 1) ≈ 0.0
end

@testset "probando la transitividad iterativa" begin
  banano = [1. 1. 1.; .214 1. .71; .27 .91 1.]
  @test norm(iteraciontransitiva(banano) - [1. 1. 1.; .27 1. .71; .27 .91 1.]) < 10^(-5)
end

@testset "coeficiente de correlacion" begin
  @test pearson([1,2,3],[1,3,2])         ≈ 1/2
  @test pearson([1,2,3,4],[1,3,2,4])     ≈ 12/15
  @test pearson([1,2,3,4,5],[1,3,2,4,5]) ≈ 9/10
  @test pearson([1,2,3,4],[3,2,1,4])     ≈ 1/5
end

@testset "coeficiente de comparabilidad" begin
  @test matriz_interseccion_rankings([[1,2,3], [1,3,2]]) |> matriz_rutas |> numeroincompatibilidades == 1
  @test matriz_interseccion_rankings([[1,2,3,4,5], [1,2,4,5,3]]) |> matriz_rutas |> numeroincompatibilidades == 2
end

