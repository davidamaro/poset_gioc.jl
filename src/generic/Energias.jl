import Primes: prime
import LinearAlgebra: norm

export energia_por_pasos_p, energia_adhoc, energia_local
export energia_no_local, energia_trucada, energia_por_pasos
export energia_por_pasos_p

const Poset = Array{T, 2} where T <: Integer

function energia_interna_posets(numero_nodos,numero_ranks,S)
  exponent = 0.0      
  U = Float64[]
  for T in 0.1:0.05:0.5
    Ener = 0
    maxL = 0.0                                           # Initialize
    #for i in 0:numero_nodos+1
    for i in 1:numero_ranks
      if S[i]!= 0 && (S[i] - Ener/T)>maxL
          maxL = S[i] - Ener/T
          Ener = Ener + numero_nodos
      end
    end
    sumdeno = 0
    sumnume = 0
    #Ener    = -2*numero_nodos
    Ener    = 0
    #for i in range(0, numero_nodos):
    #for i in 1:numero_nodos-1
    for i in 1:numero_ranks
      if S[i] != 0
          exponent = S[i] - Ener/T - maxL
      end
      sumnume += Ener*exp(exponent)
      sumdeno += exp(exponent)
      Ener     = Ener +  numero_nodos
    end
    push!(U,sumnume/sumdeno/numero_nodos)
    #U = sumnume/sumdeno/numero_nodos                    # internal energy
    #U(T)/numero_nodos
    #energ.plot(pos = (T, U) )   
  end
  intervalo_temperatura, U
end

function mi_hamming(rank, list_ranks)
  suma = 0
  for elem in list_ranks
    ja = zip(rank,elem)
    for (a,b) in ja
      if a!=b 
        suma += 1
      end
    end
  end
  suma
end


function energia_individual(rank, rankings, numero_nodos; mt = false)
  orden = Array{Int,1}[]
  for ranking in rankings
    push!(orden, sortperm(ranking))
  end

  len = rankings|>length
  num_nodos = rankings[1] |> length
  suma = 0

  if mt
    for l in orden, i in 1:length(rank)-1
        a,b = rank[i:i+1]
        if l[a] + 1 == l[b] 
          suma += (numero_nodos - i)
        end
    end
  else
    for l in orden, i in 1:length(rank)-1
        a,b = rank[i:i+1]
        if l[a] + 1 == l[b] 
          suma += 1
        end
    end
  end

  suma
end

function calcular_tabs_div(rankings, numero_nodos)
  orden = Array{Int,1}[]
  for rank in rankings
    push!(orden, sortperm(rank))
  end

  len = rankings|>length
  tab = zeros(Int,numero_nodos^2)
  for k in 1:numero_nodos^2
    for l in orden
      a,b = ind_1d_a_2d(k,numero_nodos)
      if l[a] + 1 == l[b] 
        tab[k] += 1
      end
    end
  end

  tab1 = zeros(Int,numero_nodos^2)
  for k in 1:numero_nodos^2
    tab1[k] = tab[k]-tab[ind_1d_inverso(k,numero_nodos)]
  end

  errlist = Int[]
  for k in 1:numero_nodos^2
    a,b = ind_1d_a_2d(k,numero_nodos)
    if tab[k] + tab[ind_1d_inverso(k,numero_nodos)] == len || a == b
      continue
    else
      push!(errlist, k)
    end
  end
  if length(errlist) > 0
    @warn "hay errores"
  end
  tab,tab1
end

function calcular_tabs_energia(rank, tab, numero_nodos)
  suma = 0
  for i in 1:length(rank)-1, j in i+1:length(rank)
    suma += tab[ind_2d_a_1d(rank[i], rank[j], numero_nodos)]
  end
  suma
end
function calcular_tabs(rankings, numero_nodos)
  orden = Array{Int,1}[]
  for rank in rankings
    push!(orden, sortperm(rank))
  end

  len = rankings|>length
  tab = zeros(Int,numero_nodos^2)
  for k in 1:numero_nodos^2
    for l in orden
      a,b = ind_1d_a_2d(k,numero_nodos)
      if l[a] < l[b]
        tab[k] += 1
      end
    end
  end

  tab1 = zeros(Int,numero_nodos^2)
  for k in 1:numero_nodos^2
    tab1[k] = tab[k]-tab[ind_1d_inverso(k,numero_nodos)]
  end

  errlist = Int[]
  for k in 1:numero_nodos^2
    a,b = ind_1d_a_2d(k,numero_nodos)
    if tab[k] + tab[ind_1d_inverso(k,numero_nodos)] == len || a == b
      continue
    else
      push!(errlist, k)
    end
  end
  if length(errlist) > 0
    error("hay errores")
  end
  tab,tab1
end


function ind_1d_a_2d(i, numero_nodos)
  Int(floor((i-1)/numero_nodos)+1), 1+mod(i-1,numero_nodos)
end

function ind_2d_a_1d(i,j,numero_nodos) 
  (i-1)*numero_nodos+ j
end

function ind_1d_inverso(k, numero_nodos)
  a,b = ind_1d_a_2d(k,numero_nodos)
  ind_2d_a_1d(b,a,numero_nodos)
end

function energia_individual_inicial(rank, rankings, numero_nodos)
  orden = Array{Int,1}[]
  for ranking in rankings
    push!(orden, sortperm(ranking))
  end

  len = rankings|>length
  suma = 0

  for i in 1:length(rank)-1
    for l in orden 
      a,b = rank[i:i+1]
      if l[a] + 1 == l[b] 
        suma += 1
      end
    end
    suma = len - suma
  end

  suma
end

function energia_adhoc(rank, rankings_ordenados; mt = false)

  len = rankings_ordenados|>length
  num_nodos = rank |> length
  suma = 0

  if mt
    for rank_ordered in rankings_ordenados, ind in 1:length(rank)-1
        a,b = rank[ind:ind+1]
        if rank_ordered[b] + 1 == rank_ordered[a] 
          suma += (num_nodos - ind)
        end
    end
  else
    for rank_ordered in rankings_ordenados, ind in 1:length(rank)-1
        a,b = rank[ind:ind+1]
        if rank_ordered[a] + 1 == rank_ordered[b] 
          suma -= 1
        end
        if rank_ordered[b] + 1 == rank_ordered[a] 
          suma += 1
        end
    end
  end

  suma
end

@doc Markdown.doc"""
`function energia_local(ranking_test, rankings_ordenados; mt = false)`
# Ejemplo
```
julia > energia_local([3,1,2,4,5], sortperm.([[1,2,3,4,5], [1,2,4,5,3]]), mt = true)
2
```
"""
function energia_local(ranking_test, rankings_ordenados; mt = false)
  #len = rankings_ordenados[1]|>length
  len = length(ranking_test)
  suma = 0
  if mt 
    for orden in rankings_ordenados
        for r in 1:len-1
            if !(orden[ranking_test[r]] < orden[ranking_test[r + 1]])
                suma += (len - r)
            end
        end
    end
  else
    for orden in rankings_ordenados
        for r in 1:len-1
            if !(orden[ranking_test[r]] < orden[ranking_test[r + 1]])
                suma += 1
            end
        end
    end
  end
  suma
end

function energia_local(ranking_test, matriz::Poset)
  #len = rankings_ordenados[1]|>length
  len = length(ranking_test)
  suma = 0

  for r in 1:len-1
      #if !(matriz[ranking_test[r]] < matriz[ranking_test[r + 1]])
     # if (matriz[ranking_test[r]] > matriz[ranking_test[r + 1]])
      if matriz[ranking_test[r + 1], ranking_test[r]] > 0
          suma += 1
      end
  end

  suma
end

@doc Markdown.doc"""
energia_no_local
# Ejemplos

## Energia respecto a un solo ranking
```
julia > energia_no_local([1,2,3], sortperm([1,3,2]), mt = true)
1
```
## Energia respecto a un conjunto ranking
```
julia > energia_no_local([1,2,3], sortperm.([[1,2,3], [1,3,2]]), mt = true)
1
```
## Energia respecto a un poset
```
julia > energia_no_local(sortperm([2,1,3]), lista_posets_3[2])
1
```
"""
function energia_no_local(ranking_test, matriz::Array{Int,2}; mt = false)
  len = length(ranking_test)
  suma = 0
  for i in 1:len-1, j in i+1:len
      if (matriz[ranking_test[j], ranking_test[i]] > 0)
          suma += 1
      end
  end
  suma
end

function energia_no_local(rank::Array{Int64,1}, orden::Array{Int64,1}; mt = false)
    total = 0
    n = rank |> length
    
    if mt
      for i in 1:n, j in i+1:n
          if orden[rank[i]] > orden[rank[j]]
              total += i 
          end
      end
    else
      for i in 1:n, j in i+1:n
          if orden[rank[i]] > orden[rank[j]]
              total += 1
          end
      end
    end
    total
end

function energia_no_local(rank::Array{Int64,1}, rankings_ord::Array{Array{Int64,1},1}; mt = false)
    total = 0
    n = rank |> length
    if mt
      for orden in rankings_ord
          for i in 1:n, j in i+1:n
              if orden[rank[i]] > orden[rank[j]]
                  total += i
              end
          end
      end
    else
      for orden in rankings_ord
          for i in 1:n, j in i+1:n
              if orden[rank[i]] > orden[rank[j]]
                  total += 1
              end
          end
      end
    end
    total
end

function energia_trucada(ranking_b::Array{Int64,1}, ranking_ref::Array{Int64,1})
    n = ranking_b |> length
    total = n
    pos = 1

    while pos <= n
      if ranking_b[pos] != ranking_ref[pos] 
        break
      end
      total -= 1
      pos += 1
      #@show pos
    end
    
    total
end

function energia_trucada(ranking_b::Array{Int64,1}, list_ranking_ref)
    n = ranking_b |> length
    total = n * length(list_ranking_ref)
    pos = 1

    for ranking_ref in list_ranking_ref
      while pos <= n
        if ranking_b[pos] != ranking_ref[pos] 
          break
        end
        total -= 1
        pos += 1
        #@show pos
      end
      pos = 1
    end
    
    total
end

@doc Markdown.doc"""
    energia_por_pasos(ranking_test, rankings_ordenados; mt = false)
# Ejemplos:
julia> energia_por_pasos([1,2,3], sortperm.([[1,3,2]]))
1
julia> energia_por_pasos([3,1,2,4,5], sortperm.([[1,2,3,4,5], [1,2,4,5,3]]))
6
"""
function energia_por_pasos(ranking_b::Array{Int64,1}, ranking_r_ord::Array{Int64,1})
  rank_b_ord = sortperm(ranking_b)
  len = length(ranking_b)
  total = 0
  for i in 1:len 
    if rank_b_ord[i] < ranking_r_ord[i]
      #@show rank_b_ord[i], ranking_r_ord[i], i
      #total += (len - i)
      total += (ranking_r_ord[i] - rank_b_ord[i])
      #total += (len - i)
    end
  end
  total
end

function energia_por_pasos(ranking_b::Array{Int64,1}, ranking_r_ord::Array{Array{Int64,1}, 1})
  sum((x -> energia_por_pasos(ranking_b, x)).(ranking_r_ord))
end

@doc Markdown.doc"""
    energia_por_pasos_p(ranking_test, rankings_ordenados; mt = false)
# Ejemplos:
julia> energia_por_pasos_p([1,2,3], sortperm.([[1,3,2]]))
1
julia> energia_por_pasos_p([3,1,2,4,5], sortperm.([[1,2,3,4,5], [1,2,4,5,3]]))
6
"""
function energia_por_pasos_p(ranking_b::Array{T,1}, ranking_r_ord::Array{T,1}) where T <: Integer
  rank_b_ord = sortperm(ranking_b)
  len = length(ranking_b)
  total = 0
  for i in 1:len 
    if rank_b_ord[i] < ranking_r_ord[i]
      #@show rank_b_ord[i], ranking_r_ord[i], i
      #total += (len - i)
      total += sum(prime(i) for i in 1:(ranking_r_ord[i] - rank_b_ord[i]))
      #total += (len - i)
    end
  end
  total
end

function energia_por_pasos_p(ranking_b::Array{T,1}, ranking_r_ord::Array{Array{T,1}, 1}) where T <: Integer
  sum((x -> energia_por_pasos(ranking_b, x)).(ranking_r_ord))
end

#function norma_matrices(mat::Array{T,2}) where T <: Real
#  sum(mat.*mat) |> sqrt
#end
