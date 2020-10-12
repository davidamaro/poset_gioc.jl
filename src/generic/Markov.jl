import SparseArrays: spzeros, SparseMatrixCSC
import LightGraphs: is_cyclic, DiGraph, transitiveclosure, transitivereduction, adjacency_matrix, topological_sort_by_dfs
export isacyclic, caminata_poset, caminata_poset_4
export listaposetsaleatorios
export encontrarminimo
import LinearAlgebra: norm
export caminatale

function complement(x::Array{Int64,1}, y::Array{Int64,1})
    todos = vcat([x, y]...) |> unique
    filter(z -> (in(z, x) && !in(z, y)), todos)
end

function isacyclic!(mat::SparseMatrixCSC{Int64,Int64})
    len = mat.n
    nodos = 1:len |> collect

    while length(nodos) > 0
        valoresnulos = findall(!iszero, mat)
        nodos_consalida = Int[]
        for v in valoresnulos
            push!(nodos_consalida, v.I[1])
        end
        nodos_consalida|>unique!
        malo = complement(nodos, nodos_consalida)
        if length(malo) == 0
            return false
        end
        nodos = filter!(x -> x != malo[1], nodos)

        restantes = (filter(x -> (x.I[2] == malo[1]) || (x.I[1] == malo[1]), valoresnulos))
        for r in restantes
            mat[r.I...] = 0
        end
    end
    return true
end

function isacyclic(matinput::SparseMatrixCSC{Int64,Int64})
    mat = deepcopy(matinput)
    len = mat.n
    nodos = 1:len |> collect

    while length(nodos) > 0
        valoresnulos = findall(!iszero, mat)
        nodos_consalida = Int[]
        for v in valoresnulos
            push!(nodos_consalida, v.I[1])
        end
        nodos_consalida|>unique!
        malo = complement(nodos, nodos_consalida)
        if length(malo) == 0
            return false
        end
        nodos = filter!(x -> x != malo[1], nodos)

        restantes = (filter(x -> (x.I[2] == malo[1]) || (x.I[1] == malo[1]), valoresnulos))
        for r in restantes
            mat[r.I...] = 0
        end
    end
    return true
end

function aristarandom(n::Int64)
    i = rand(1:n)
    remanentes = filter!(x -> x!=i, 1:n |> collect)
    j = rand(remanentes)
    i,j
end

function cardinality(dag::Array{Int64,2})
    t = dag |> matriz_rutas |> sum #|> DiGraph |> transitiveclosure |> adjacency_matrix |> sum
    r = dag |> reduccion_transitiva |> sum#|> DiGraph |> transitivereduction |> adjacency_matrix |> sum
    2^(t - r)
end

function cardinality(dag::SparseMatrixCSC{Int64,Int64})
    t = dag |> matriz_rutas |> sum #|> DiGraph |> transitiveclosure |> adjacency_matrix |> sum
    r = dag |> reduccion_transitiva |> sum#|> DiGraph |> transitivereduction |> adjacency_matrix |> sum
    2^(t - r)
end

@doc Markdown.doc"""
    caminata_poset(numeronodos, pasos; verbose)
> Devuelve un DAG

# Examples:
```
julia> caminata_poset(3,100)
```
"""
function caminata_poset(n::Int64,pasos::Int64;verbose::Bool=false)
    original  =  spzeros(Int64, n,n)
    operacion =  spzeros(Int64, n,n)
    for i in 1:pasos
        x,y = aristarandom(n)
        if verbose
            @show x,y
        end

        if original[x,y] == 0
            operacion[x,y] = 1
            if !isacyclic(operacion)
                operacion[x,y] = 0
            end
        else
            operacion[x,y] = 0
        end


        cor,cop = cardinality.([original,operacion])

        proba   = rand()

        if verbose
            @show original, operacion
        end

        if proba < minimum([1.0, cor/cop])
            original  = operacion |> copy
        else
            operacion = original |> copy
        end
    end
    original
end

function caminata_poset_4(n::Int64,pasos::Int64;verbose::Bool=false)
    original  =  spzeros(Int64, n,n)
    #operacion =  spzeros(Int64, n,n)
    operacion = original |> deepcopy
    for i in 1:pasos
        
        x,y = rand([(1,2), (1,3), (1,4), (2,1), (2,3), (2,4), (3,1), (3,2), (3,4), (4,1), (4,2), (4,3)])
        if verbose
            @show x,y
        end

        if original[x,y] == 0
            operacion[x,y] = 1
            if is_cyclic(operacion |> DiGraph)
                operacion[x,y] = 0
            end
        else
            operacion[x,y] = 0
        end


        cor,cop = cardinality.([original,operacion])

        proba = rand()

        if verbose
            @show norm(original - operacion), cor/cop
        end

        if proba < minimum([1.0, cor/cop])
            original = deepcopy(operacion)
        else
            operacion = deepcopy(original)
        end
    end
    original
end

@doc Markdown.doc"""
    caminata_poset(numeronodos, pasos; verbose)
> Devuelve un DAG

# Examples:
```
julia> listaposetsaleatorios(3,10)
```
"""
function listaposetsaleatorios(n::Int64, steps::Int64; m::Int64 = n^2)
  lista_posetsrandom = SparseMatrixCSC{Int64,Int64}[]
  for i in 1:steps
    push!(lista_posetsrandom, caminata_poset(n,n^2))
  end
  lista_posetsrandom
end

@doc Markdown.doc"""
    encontrarminimo(mat::Array{Float64,2},n::Int64, tope::Int64)

> Utilizado un numero de pruebas `tope` encuentra el minimo usando 
> posets aleatorios generados usando `caminata_poset`.

# Examples:
```
julia> mat = [0.3, 0.3 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0];
julia> encontrarminimo(mat, 4, 100)
```
"""
function encontrarminimo(mat::Array{Float64,2},n::Int64, tope::Int64)
    minimo::Float64 = 100.
    bono::SparseMatrixCSC{Int64,Int64} = caminata_poset(n, n^2)
    for i in 1:tope
        pp = caminata_poset(n, n^2) |> reduccion_transitiva
        nn = norm(mat - pp)
        if minimo > nn#norm(union_promedio - pp)
            minimo = nn
            bono = pp |> copy
        end
    end
    minimo, bono
end


@doc Markdown.doc"""
    caminatale(mat::Array{Int64,1}, pasos::Int64)
> Return a random linear extension as an Array{Int64,1}.

# Example
```
julia> mat = [0 1 1; 0 0 0; 0 0 0]
julia> caminatale(mat, 200)

julia> posetaleatorio = caminata_poset(3,200)
julia> linext_aleatoria = posetaleatorio |> Array |> x->caminatale(x, 200)
```
"""
function caminatale(mat::Array{Int64,2}, pasos::Int64)
    le::Array{Int64,1} = mat |> DiGraph |> topological_sort_by_dfs
    n::Int64 = le |> length
    matrizrutas::Array{Int64,2} = mat |> matriz_rutas
    tmp::Array{Int64,1} = copy(le)

    for _ in 1:pasos
        p = rand([0,1])
        if p == 0
            i,j = aristarandom(n)

            t = tmp[i]
            tmp[i] = tmp[j]
            tmp[j] = t

            if islinearextension(matrizrutas, tmp)
                le = copy(tmp)
            end
        else
            continue
        end
    end
    le
end

function islinearextension(mat::Array{Int64,2} #= matriz de rutas =#, lista::Array{Int64,1})
    n::Int64 = lista |> length

    for i in 1:n-1, j in i+1:n
        if mat[lista[j], lista[i]] == 1
            return false
        end
    end

    return true
end
