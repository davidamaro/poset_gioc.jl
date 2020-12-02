export monte_carlo, wang_landau, densidad_exacta, positive
export Simulacion, crear_matriz, new_mc
export matriz_rutas, reduccion_transitiva
export ranking_natural, matriz_union_rankings
export rankings_random, determinar_minimo_6, determinar_minimo_5
export matriz_interseccion_rankings
export condorcet
export determinar_minimos
export derecha_abajo
export posicionpromedio, posicionvarianza
export pearson
export sensibilidad
export generarmatriz, iteraciontransitiva, convertidor
export filter
export numeroincompatibilidades, gradocoincidencia1
export m², m³, pareja_matrizadyacencia, mn
export fuzzy
export equivalencias, membresia

import Statistics: mean, median
import Combinatorics: permutations
import Base.+
import LinearAlgebra: norm
import SparseArrays: spzeros, SparseMatrixCSC

function +(a::Array{Int64,2},b::Tuple{Int64,Int64})
    x,y = b
    if a[y,x] == 0
        a[x,y] += 1
    else
        a[x,y] = 0
        a[y,x] = 0
    end
    a
end

@doc Markdown.doc"""
    energia_por_pasos_p(ranking_test, rankings_ordenados; mt = false)
# Ejemplos:
```
julia> sim_a = Simulacion(5000
                         ,50
                         ,50000
                         ,1
                         ,1
                         ,0.2
                         ,0.5
                         ,energia_no_local
                         ,1
                         ,1)
```
"""
mutable struct Simulacion
  tmax::Int
  tdiscard::Int
  ntprint::Int
  dtprint::Int
  jprint::Int
  beta::Float64
  problarge::Float64

  metodo_energia::Function
  minimo_energia::Int
  maximo_energia::Int
  fac::Float64
  fmin::Float64
end

################################################################################
#                                                                              #
#                                Monte Carlo                                   #
#                                                                              #
################################################################################
function monte_carlo(lista_ranks_ordenados,sim::Simulacion)
  t = 0.
  #tmax = 5000.
  tmax = sim.tmax
  numero_nodos = lista_ranks_ordenados |> first |> length
  var = collect(1:numero_nodos)
  #problarge = 0.5
  problarge = sim.problarge
  tabtotal = zeros(Int,numero_nodos^2)
  tdiscard = sim.tdiscard
  ntprint = sim.ntprint
  dtprint = sim.dtprint
  jprint = sim.jprint
  tprint = [tdiscard + l*dtprint for l in 1:ntprint]

  metodo_energia = sim.metodo_energia

  while (t<=tmax)
    ind = rand(1:numero_nodos-1)
    prob = rand()

    if prob <= problarge
      ind1 = rand(ind+1:numero_nodos)
    else
      ind1 = ind+1
    end
    # energia es la suma de las diferencias
    # entre las energias del ranking volteado
    # y el viejo
     v = deepcopy(var)
     tmp = v[ind]
     v[ind] = v[ind1]
     v[ind1] = tmp
     energia = metodo_energia(v, lista_ranks_ordenados)
               - metodo_energia(var, lista_ranks_ordenados)
    #

#    v = var[ind:ind1]
#    temp = [tab1[f21(v[i], v[i+1], numero_nodos)] for i in 1:length(v)-1] 
#    energia = temp |> sum 
    delte = sim.beta*energia
    prob = rand()

    if delte <= 0. || prob <= exp(-delte)
#      var[ind] = v[end]
#      var[ind1] = v[1]
        var = v
    end

    t = t+1.0/numero_nodos
    if t-tdiscard >= tprint[jprint]
      jprint += 1
      for k in 1:numero_nodos^2
        l1, l2 = ind_1d_a_2d(k,numero_nodos)
        if sortperm(var)[l1] < sortperm(var)[l2]
          tabtotal[k] += 1
        end
      end
    end
  end
  tabtotal
end
 
################################################################################
#                                                                              #
#                                Wang Landau                                   #
#                                                                              #
################################################################################
function wang_landau(lista_ranks_ordenados, sim::Simulacion)
  numero_nodos = lista_ranks_ordenados |> first |> length
  orden_normal = collect(1:numero_nodos)
  var = orden_normal[lista_ranks_ordenados[1]]
  len = lista_ranks_ordenados |> length

  metodo_energia = sim.metodo_energia
  

  ener    = metodo_energia(var, lista_ranks_ordenados)
  ntot    = sim.maximo_energia#10000 #energia maxima
  dtprint = sim.dtprint
  hist    = zeros(Int,ntot)
  hist1   = zeros(Int,ntot)
  offset  = sim.minimo_energia
  weight  = zeros(ntot)

  hist[ener + offset]   = 1
  weight[ener + offset] = 1

  jprint = sim.jprint
  t      = 0
  fac    = sim.fac 
  f      = 1.
  fmin   = sim.fmin 

  while f > fmin
    #indice_aleatorio = rand(1:len)
    t    = t + 1/(len*numero_nodos)
    ind  = rand(1:numero_nodos-1)
    ind1 = ind + 1

    helper       = deepcopy(var)
    tmp          = helper[ind]
    helper[ind]  = helper[ind1]
    helper[ind1] = tmp
    ener1        = ener + ( metodo_energia(helper, lista_ranks_ordenados) - metodo_energia(var   , lista_ranks_ordenados) )


    if ener <= 0 || ener > ntot || ener1 <= 0 || ener1 > ntot
      @show ener <= 0 , ener > ntot , ener1 <= 0 , ener1 > ntot, ener1, ntot
      error("se saldra de los bordes. mala energia")
    end

    delte = weight[ener1 + offset] - weight[ener + offset]

    if delte <= 0
      var  = helper
      ener = ener1
    elseif rand() <= exp(-delte)
      var  = helper
      ener = ener1
    end

    if ener <= 0 || ener > ntot
      error("energia fuera de borde")
    end

    weight[ener+offset] += f
    hist[ener+offset]   += 1

    if t >= jprint*dtprint
      jprint += 1
      hist1 = filter(x -> x > 0, hist)
      if minimum(hist1) > fac*mean(hist1)
        f = f/2.
        hist = zeros(Int,ntot)
        println("cambia f:"*string(f))
      end
    end
  end
  weight1 = filter(positive, weight)

  weight1, hist1
end

################################################################################
#                                                                              #
#                                Densidad exacta                               #
#                                                                              #
################################################################################
@doc Markdown.doc"""
    densidad_exacta(ranking_test, max, metodo)
`max` corresponde a la energia maxima (para un solo ranking).
'metodo' es una funcion.

# Ejemplos:
```
dens = densidad_exacta(sortperm.([[1,2,3,4,5], [5,1,3,4,2]]), 17, energia_por_pasos_p)
dens = densidad_exacta(matriz_poset, 17, energia_por_pasos_p)
```
"""
function densidad_exacta(lista_ranks_ordenados::Array{Array{T,1},1},max::T,metodo::Function) where T<:Integer
  num_objetos = lista_ranks_ordenados |> first |> length
  len = length(lista_ranks_ordenados)
  densidad = zeros(T,len*max);
  for i in permutations(1:num_objetos)
      ener = metodo(collect(i), lista_ranks_ordenados)
      densidad[ener+1] += 1
  end
  densidad
end

function densidad_exacta(matriz_poset::Array{T,2},max::T,metodo::Function) where T<:Integer
  num_objetos,num_objetos  = matriz_poset |> size
  densidad = zeros(T,num_objetos*max);
  matriz_poset_ruta = matriz_poset |> matriz_rutas
  for i in permutations(1:num_objetos)
      ener = metodo(collect(i), matriz_poset_ruta)
      densidad[ener+1] += 1
  end
  densidad
end
################################################################################
#                                                                              #
#                                Misc                                          #
#                                                                              #
################################################################################
positive(x::T) where T<:Real = x > 0
################################################################################
#                                                                              #
#                                Reduccion transitiva                          #
#                                                                              #
################################################################################
#source https://github.com/lucasrabiec/TransitiveAlgorithms/blob/master/src/TransitiveAlgorithms/Reduction.cs
function matriz_rutas(mat::SparseMatrixCSC{Int64,Int64})
    path_matrix = copy(mat)
    n,n = size(mat)
    for i in 1:n, j in 1:n
        if i == j
            continue
        end
        if path_matrix[j,i] > 0
            for k in 1:n
                if path_matrix[j,k] == 0
                    path_matrix[j,k] = path_matrix[i,k]
                end
            end
        end
    end
    path_matrix
end

function matriz_rutas(mat::Array{T,2}) where T <: Integer
    path_matrix = deepcopy(mat)
    n,n = size(mat)
    for i in 1:n, j in 1:n
        if i == j
            continue
        end
        if path_matrix[j,i] > 0
            for k in 1:n
                if path_matrix[j,k] == 0
                    path_matrix[j,k] = path_matrix[i,k]
                end
            end
        end
    end
    path_matrix
end

function reduccion_transitiva(mat::SparseMatrixCSC{Int64,Int64})
    mat_red = mat |> matriz_rutas
    n,n = size(mat)

    for i in 1:n, j in 1:n
        if mat_red[j,i] > 0
            for k in 1:n
                if mat_red[i,k] > 0
                    mat_red[j,k] = 0
                end
            end
        end
    end
    mat_red
end

function reduccion_transitiva(mat::Array{T,2}) where T <: Integer
    mat_red = mat |> matriz_rutas
    n,n = size(mat)

    for i in 1:n, j in 1:n
        if mat_red[j,i] > 0
            for k in 1:n
                if mat_red[i,k] > 0
                    mat_red[j,k] = 0
                end
            end
        end
    end
    mat_red
end
################################################################################
#                                                                              #
#                                Nuevo MonteCarlo                              #
#                                                                              #
################################################################################
function crear_matriz(orden::Array{T,1}) where T <: Integer
  n = orden |> length
  mat = zeros(Int, (n,n))
  #natural = 1:n |> collect

  for i in 1:n, j in i+1:n
    #tmp = natural[orden]
    #mat[tmp[i], tmp[j]] += 1
    mat[orden[i], orden[j]] += 1
  end

  mat
end

function matriz_union_rankings(lista_ranks::Array{Array{T, 1}, 1}; binario = true, promediado = false) where T <: Real
  len = length(lista_ranks)
  #mat = zeros(Float64, (len, len))
  mat = sum(crear_matriz.(lista_ranks)).*(1.0)
  if binario && !promediado
    for (ind,val) in enumerate(mat)
      if val > 0
        mat[ind] = 1
      else
        mat[ind] = 0
      end
    end
  elseif promediado && !binario
    for (ind,val) in enumerate(mat)
      if val > 0
        mat[ind] /= len
      else
        mat[ind] = 0
      end
    end
  elseif !promediado && !binario
    return mat
  else
    throw(ArgumentError("keywords invalidas"))
  end
  mat
end

function matriz_interseccion_rankings(lista_ranks::Array{Array{T, 1}, 1}; porcentaje = 1.0, binario = true) where T <: Integer
  ranks = unique(lista_ranks)
  len = ranks |> length
  mat = sum(crear_matriz.(ranks))
  if binario
    for (ind,val) in enumerate(mat)
      if val < len*porcentaje
        mat[ind] = 0
      else
        mat[ind] = 1
      end
    end
  else
    for (ind,val) in enumerate(mat)
      if val < len*porcentaje
        mat[ind] = 0
      end
    end
  end
  mat
end

#function crear_matriz(lista_ranks_ordenados::Array{Array{T, 1}, 1}) where T <: Integer
#  n = lista_ranks_ordenados[1] |> length
#  mat = zeros(Int, (n,n))
#  natural = 1:n |> collect
#  for orden in lista_ranks_ordenados
#    tmp = natural[orden]
#    for i in 1:n-1, j in i+1:n
#      mat[tmp[i], tmp[j]] = 1
#      if mat[tmp[i], tmp[j]] > 0 && mat[tmp[j],tmp[i]] > 0
#        mat[tmp[i],tmp[j]] = 0
#        mat[tmp[j],tmp[i]] = 0
#      end
#    end
#  end
#  mat
#end


function new_mc(lista_rankings, sim::Simulacion)
  t = 0.
  tmax = sim.tmax
  problarge = sim.problarge
  tdiscard = sim.tdiscard
  ntprint = sim.ntprint
  dtprint = sim.dtprint
  jprint = sim.jprint
  tprint = [tdiscard + l*dtprint for l in 1:ntprint]

  numero_nodos = lista_rankings[1] |> length
  
  distancias = crear_matriz(lista_rankings) |> reduccion_transitiva
  original = deepcopy(distancias)

  rank_paso = deepcopy(distancias)
  while (t<=tmax)
    ind = rand(1:numero_nodos-1)
    ind1 = rand(filter(x -> x != ind, 1:numero_nodos|>collect))
    @show ind, ind1

    #rank_paso[ind, ind1] += 1 # esta adiciones debe ser mas controlada
    +(rank_paso, (ind,ind1))
    n_paso       = norma_matrices(rank_paso  |> reduccion_transitiva)  + norma_matrices(reduccion_transitiva(rank_paso)  - original)*numero_nodos
    n_distancias = norma_matrices(distancias |> reduccion_transitiva) + norma_matrices(reduccion_transitiva(distancias)  - original)*numero_nodos
    energia = n_paso - n_distancias

    delte = sim.beta*energia
    prob = rand()

    if delte <= 0. || prob <= exp(-delte)
      distancias = rank_paso
    end

    t = t+1.0/numero_nodos
    if t-tdiscard >= tprint[jprint]
      jprint += 1
    end
  end
  distancias
end

function ranking_natural(n::Int)
  @assert n >= 2
  1:n |> collect
end

#############################################################
#                                                            
#                                                            
#              Rankings                                      
#                                                            
#                                                            
#############################################################
function permutacion_corta(ranking::Vector{T}) where T <: Integer
    n = ranking |> length
    i = rand(1:n-1)
    xx = deepcopy(ranking)
    tmp = xx[i]
    xx[i] = xx[i+1]
    xx[i+1] = tmp
    
    xx
end
function permutacion_larga(ranking::Vector{T}, pasos::T) where T <: Integer
    n = ranking |> length
    xx = deepcopy(ranking)
    for ii in 1:pasos
        i = rand(1:n-1)
        tmp = xx[i]
        xx[i] = xx[i+1]
        xx[i+1] = tmp
    end

    xx
end

function rankings_random(n::T, estructura::Vector{T}; ini = ranking_natural(n), pasos = 1000) where T <: Integer
    lista = [ini]
    while length(estructura) > 0
        s = pop!(estructura)
        for j in 1:s-1
            push!(lista, permutacion_corta(lista[end]))
        end
        if length(estructura) > 0
            push!(lista, permutacion_larga(lista[end], pasos))
        end
    end
    lista
end

function determinar_minimo_5(ranking_prueba)
    min = norma_matrices(ranking_prueba)
    ind_min = 0
    for (ind, poset) in enumerate(lista_posets_e_5)
        tmp = norma_matrices(ranking_prueba - poset)
        if tmp < min
            min = tmp
            ind_min = ind
        end
    end
    (min, ind_min)
end

function determinar_minimo_6(ranking_prueba; todos = false)
    min = norm(ranking_prueba)
    ind_min = 0
    for (ind, poset) in enumerate(lista_posets_e_6)
        tmp = norm(ranking_prueba - poset)
        if tmp < min
            min = tmp
            ind_min = ind
        end
    end
    (min, ind_min)
end

function determinar_minimos(ranking_prueba, lista)
    min = norm(ranking_prueba)
    ind_min = 0
    for (ind, poset) in enumerate(lista)
        tmp = norm(ranking_prueba - poset)
        if tmp < min
            min = tmp
            ind_min = ind
        end
    end
    lista_posetsminimos = Int[]
    for (i,v) in enumerate(lista)
        nn = norm(v - ranking_prueba)
        if nn ≈ min
            push!(lista_posetsminimos, i)
        end
    end
    return lista_posetsminimos
end

function condorcet(matriz::Array{T,2}) where T <: Real
  matriz_condorcet = zero(matriz)
  for (i,v) in enumerate(matriz)
    if v < .5
        matriz_condorcet[i] = 0
    else
        matriz_condorcet[i] = 1
    end
  end
  matriz_condorcet
end

function derecha_abajo(test,lista)
    alto    = length(test)
    derecha = 0
    indice = 1
    for (ind,poset) in enumerate(lista)
        if norm(poset - test) < alto && norm(poset) > derecha
            derecha = norm(poset)
            alto =  norm(poset - test)
            indice = ind
        end
    end
    indice
end

"""
```
posicionpromedio(listarankings::Array{Array{Int64,1},1},m::Int64)
```
listarankings es una lista enteros, correspondiente a los rankings.
m es un entero, que corresponde al nodo del que se interesa calcular
la posicion promedio.
"""
function posicionpromedio(listarankings::Array{Array{Int64,1},1},m::Int64)
    n::Int64 = length(listarankings[1])
    posicionesm::Array{Int64,1} = zeros(Int,n)
    for el in listarankings
        posicionesm[sortperm(el)[m]] += 1
    end
    sum(collect(1:n) .* (posicionesm/length(listarankings)))
end

"""
```
posicionvarianza(listarankings::Array{Array{Int64,1},1},m::Int64)
```
listarankings es una lista enteros, correspondiente a los rankings.
m es un entero, que corresponde al nodo del que se interesa calcular
la varianza de su posicion.
"""
function posicionvarianza(listarankings::Array{Array{Int64,1},1},m::Int64)
    n::Int64 = length(listarankings[1])
    posicionesm::Array{Int64,1} = zeros(Int,n)
    for el in listarankings
        posicionesm[sortperm(el)[m]] += 1
    end
    sqrt(sum(collect(1:n).^2 .* (posicionesm/length(listarankings))) - (sum(collect(1:n) .* (posicionesm/length(listarankings))))^2)
end

@doc Markdown.doc"""
    pearson(r1::Array{Int64,1}, r2::Array{Int64,1})
> Calcula el coeficiente de Pearson entre dos rankings.
> $d(i) = R_1(i) - R_2(i)$ en donde $R_j(i)$  es la posicion del nodo i en el ranking j.
> El coeficiente es calculado como $6 \Sigma/(n (n^2 -1))$, con $\Sigma = \sum d^2$.

# Examples:
```
julia> r1 = [1,2,3]
julia> r2 = [2,1,3]
julia> pearson(r1,r2)
```
"""
function pearson(r1::Array{Int64,1}, r2::Array{Int64,1})
    len::Int64 = r1 |> length
    1 - 6*sum([(sortperm(r1)[i] - sortperm(r2)[i])^2 for i in 1:len])/(len*(len^2 - 1))
end

@doc Markdown.doc"""
    sensibilidad(rankingscompleto::Array{Array{Int64,1},1}, rankingsincompleto::Array{Array{Int64,1},1})
> Calcula el coeficiente de Pearson entre dos rankings.
> $d(i) = R_1(i) - R_2(i)$ en donde $R_j(i)$  es la posicion del nodo i en el ranking j.
> El coeficiente es calculado como $6 \Sigma/(n (n^2 -1))$, con $\Sigma = \sum d^2$.

# Examples:
```
julia> r1 = [[1,2,3], [1,3,2], [3,1,2]]
julia> r2 = [[1,2,3], [1,3,2]]
julia> sensibilidad(r1,r2)
```
"""
function sensibilidad(rankingscompleto::Array{Array{Int64,1},1}, rankingsincompleto::Array{Array{Int64,1},1})
    mat1::Array{Int64,2} = matriz_interseccion_rankings(rankingscompleto)
    mat2::Array{Int64,2} = matriz_interseccion_rankings(rankingsincompleto)
    len::Int64 = length(mat1)
    suma::Int64 = 0

    for i in 1:len
        if (mat1[i] == 1) ⊻ (mat2[i] == 1)
            suma += 1
        end
    end

    suma
end

function convertidor(lista::Array{Int64,1})
    uno = lista |> copy
    n::Int64 = uno |> length
    @simd for i in eachindex(uno)
        @inbounds uno[i] = n + 1 - uno[i]
    end
    uno
end
function convertidor(lista::Array{Array{Int64,1},1})
    convertidor.(lista)
end
function convertidor!(uno::Array{Int64,1})
    n::Int64 = uno |> length
    @simd for i in eachindex(uno)
        @inbounds uno[i] = n + 1 - uno[i]
    end
    uno
end

function convertidor!(lista::Array{Array{Int64,1},1})
    convertidor!.(lista)
end

function generarmatriz(listarankings::Array{Array{Int64,1},1})
    n = listarankings[1] |> length
    #listarankings = sortperm.(listarankings)
    ordenados = convertidor(listarankings)
    mat = zeros(Float64, n, n)
    for i in 1:n, j in 1:n
        @inbounds mat[i,j] = porranking(i,j, ordenados)
    end
    mat
end

function porranking(a::Int64,b::Int64,ordenados::Array{Array{Int64,1},1})
    #ordenados = sortperm.(listarankings)
    sum([min(x[a], x[b]) for x in ordenados])/sum([x[a] for x in ordenados])
end

function iteraciontransitiva(input)
    mat = input |> similar
    n,_ = size(input)
    listplaceholder = zeros(Float64, n)
    for x in 1:n, y in 1:n
        listplaceholder = zeros(Float64, n)
        for w in 1:n
            listplaceholder[w] = min(input[x,w], input[w,y])
        end
        mat[x,y] = maximum(listplaceholder)
    end
    mat
end

function Base.filter(predicado, mat::Array{Int64,2})
    for i in eachindex(mat)
        if !predicado(mat[i])
            mat[i] = 0.0
        end
    end
    mat
end
function Base.filter(predicado, mat::Array{Float64,2})
    for i in eachindex(mat)
        if !predicado(mat[i])
            mat[i] = 0.0
        end
    end
    mat
end

function numeroincompatibilidades(mat::Array{Int64,2})
    incom = 0
    n,_=size(mat)
    for i in 1:n-1, j in i+1:n
        if mat[i,j]+mat[j,i] != 1
            incom +=1
        end
    end
    incom
end
function gradocoincidencia1(listaranks)
    @assert length(listaranks) == 2
    n = listaranks[1] |> length
    1 - (matriz_interseccion_rankings(listaranks) |> numeroincompatibilidades)/binomial(n,2)
end

## metodos m₂ y m₃

function m²(listarankings; alter::Bool = false)
    n = listarankings[1] |> length
    mat = zeros(Float64, n,listarankings |> length)
    output = zeros(Float64,n,2)
    if !alter
      for (i,l) in enumerate( listarankings )
          mat[:,i] = [1/x for x in sortperm(l)]
      end
    else
      for (i,l) in enumerate( listarankings )
          mat[:,i] = [n-x for x in sortperm(l)]
      end
    end
    for i in 1:n
        x,y = extrema(mat[i,:]) 
        output[i,1] = y
        output[i,2] = x
    end
    output
end

function m³(listarankings; alter::Bool = false)
    n = listarankings[1] |> length
    mat = zeros(Float64, n,listarankings |> length)
    output = zeros(Float64,n,3)
    if !alter
      for (i,l) in enumerate( listarankings )
          mat[:,i] = [1/x for x in sortperm(l)]
      end
    else
      for (i,l) in enumerate( listarankings )
          mat[:,i] = [n - x for x in sortperm(l)]
      end
    end
    for i in 1:n
        x,y = extrema(mat[i,:]) 
        m = median(mat[i,:])
        output[i,1] = y
        output[i,2] = m
        output[i,3] = x
    end
    output
end

function mn(listarankings; alter::Bool = false)
    n = listarankings[1] |> length
    mat = zeros(Float64, n,listarankings |> length)
    output = zeros(Float64,n,listarankings |> length)
    if !alter
      for (i,l) in enumerate( listarankings )
          mat[:,i] = [1/x for x in sortperm(l)]
      end
    else
      for (i,l) in enumerate( listarankings )
          mat[:,i] = [n - x for x in sortperm(l)]
      end
    end
    for i in 1:n
#        x,y = extrema(mat[i,:]) 
#        m = median(mat[i,:])
#        output[i,1] = y
#        output[i,2] = m
#        output[i,3] = x
        output[i,:] .= mat[i,:]
    end
    output
end

function pareja_matrizadyacencia(matrizparejas)
    n,_ = size(matrizparejas)
    output = zeros(Int64, n,n)
    for i in 1:n, j in 1:n
        if i == j
            continue
        end
        if all(matrizparejas[i,:] .>= matrizparejas[j,:]) && !all(matrizparejas[i,:] .== matrizparejas[j,:])
            output[i,j] = 1
        end
    end
    output
end

function fuzzy(p)
    n,m = size(p)
    matfuzzy = zeros(Float64, n,n)
    for i in 1:n, j in 1:n
        hh = sum([p[j,k] for k in 1:m])
        if hh ≈ 0.0
          matfuzzy[i,j] = 1.0
        else
          matfuzzy[i,j] = sum([min(p[i,k], p[j,k]) for k in 1:m])/hh
        end
    end
    matfuzzy
end

@doc Markdown.doc"""
equivalencia(matriz de propiedades calculada con mn)
calcula las clases de equivalencia usando la matriz
de propiedades. Usa sort.
"""
function equivalencias(propiedades)
    numnodos, numpropo = propiedades |> size
    listasimilaridades::Array{Array{Int,1},1} = [Int[] for _ in 1:numnodos]
    for i in 1:numnodos
        for j in 1:numnodos
            if norm( sort(propiedades[i,:] )- sort(propiedades[j,:] )) < 10^(-5)
                push!(listasimilaridades[i], j)
            end
        end
    end
    unique(sort.(listasimilaridades))
end

@doc Markdown.doc"""
membresia(extension lineal, poset en forma de matriz de adyacencia)
calcula el grado de membresia de una extension lineal con respecto
a un fuzzy poset.
"""
function membresia(el, mposet)
    n, _ = size(mposet)
    @assert n == length(el)
    or = el |> sortperm
    tope = 0.
    for u in 1:n, v in 1:n
        if u == v
            continue
        end
    #for u in 1:n, v in u+1:n
        este =min(mposet[u,v], 1 - (or[u] < or[v]) )
        #@show este, mposet[u,v]
        if tope < este
            tope = este
        end
        if tope ≈ 1.
            return 0
        end
    end

    1-tope
    #tope
end
