export monte_carlo, wang_landau, densidad_exacta, positive
export Simulacion, crear_matriz, new_mc
export matriz_rutas, reduccion_transitiva
import Statistics: mean
import Combinatorics: permutations

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
```
"""
function densidad_exacta(lista_ranks_ordenados::Array{Array{T,1},1},max::T,metodo::Function) where T<:Integer
  limite = lista_ranks_ordenados |> first |> length
  len = length(lista_ranks_ordenados)
  densidad = zeros(T,len*max);
  for i in permutations(1:limite)
      ener = metodo(collect(i), lista_ranks_ordenados)
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
  natural = 1:n |> collect

  for i in 1:n, j in i+1:n
    tmp = natural[orden]
    mat[tmp[i], tmp[j]] += 1
  end

  mat
end

function crear_matriz(lista_ranks_ordenados::Array{Array{T, 1}, 1}) where T <: Integer
  n = lista_ranks_ordenados[1] |> length
  mat = zeros(Int, (n,n))
  natural = 1:n |> collect
  for orden in lista_ranks_ordenados
    tmp = natural[orden]
    for i in 1:n-1, j in i+1:n
      mat[tmp[i], tmp[j]] = 1
      if mat[tmp[i], tmp[j]] > 0 && mat[tmp[j],tmp[i]] > 0
        mat[tmp[i],tmp[j]] = 0
        mat[tmp[j],tmp[i]] = 0
      end
    end
  end
  mat
end


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
  
  distancias = crear_matriz(lista_rankings) |> reduccion_transitiva#zeros(Int, (numero_nodos, numero_nodos))
  original = deepcopy(distancias)

  rank_paso = deepcopy(distancias)
  while (t<=tmax)
    ind = rand(1:numero_nodos-1)
    ind1 = rand(filter(x -> x != ind, 1:numero_nodos|>collect))
    prob = rand()

#    if prob <= problarge
#      ind1 = rand(ind+1:numero_nodos)
#    else
#      ind1 = ind+1
#    end
    # energia es la suma de las diferencias
    # entre las energias del ranking volteado
    # y el viejo
    rank_paso[ind, ind1] += 1
    n_paso       = norma_matrices(rank_paso |> reduccion_transitiva)  + norma_matrices(rank_paso - original)*0.5
    n_distancias = norma_matrices(distancias |> reduccion_transitiva) + norma_matrices(distancias - original)*0.5
    energia = n_paso - n_distancias
    #

#    v = var[ind:ind1]
#    temp = [tab1[f21(v[i], v[i+1], numero_nodos)] for i in 1:length(v)-1] 
#    energia = temp |> sum 
    delte = sim.beta*energia
    prob = rand()

    if delte <= 0. || prob <= exp(-delte)
#      var[ind] = v[end]
#      var[ind1] = v[1]
      distancias = rank_paso
    end

    t = t+1.0/numero_nodos
    if t-tdiscard >= tprint[jprint]
      jprint += 1
    end
  end
  distancias
end
