import SparseArrays: spzeros, SparseMatrixCSC
import LightGraphs: is_cyclic, DiGraph, transitiveclosure, transitivereduction, adjacency_matrix
export isacyclic, caminata_poset, caminata_poset_4
import LinearAlgebra: norm

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
    t = dag |> DiGraph |> transitiveclosure |> adjacency_matrix |> sum
    r = dag |> DiGraph |> transitivereduction |> adjacency_matrix |> sum
    2^(t - r)
end

function cardinality(dag::SparseMatrixCSC{Int64,Int64})
    t = dag |> DiGraph |> transitiveclosure |> adjacency_matrix |> sum
    r = dag |> DiGraph |> transitivereduction |> adjacency_matrix |> sum
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

        proba = rand()

        if verbose
            @show original, operacion
        end

        if proba < minimum([1.0, cor/cop])
            original = operacion
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
            #@show original, operacion, cor/cop
        end

        if proba < minimum([1.0, cor/cop])
            original = deepcopy(operacion)
        else
            operacion = deepcopy(original)
        end
    end
    original
end
