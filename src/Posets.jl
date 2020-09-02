module Posets

using DelimitedFiles
using SimplePosets
using Statistics
using Random
using Combinatorics
using Markdown
using Primes

include("Generic.jl")

import .Generic:  energia_por_pasos_p, energia_adhoc, energia_local,
                  energia_no_local, energia_trucada, energia_por_pasos,
                  energia_por_pasos_p,
                  monte_carlo, wang_landau, Simulacion, densidad_exacta,
                  positive,
                  SG, apply_filter,
                  crear_matriz, norma_matrices, new_mc,
                  matriz_rutas, reduccion_transitiva,
                  lista_posets_3, lista_posets_4, lista_posets_5,
                  lista_posets_e_4, lista_posets_e_5, lista_posets_e_6 ,
                  ranking_natural
export energia_por_pasos_p, energia_adhoc, energia_local
export energia_no_local, energia_trucada, energia_por_pasos
export energia_por_pasos_p
export monte_carlo, wang_landau, Simulacion, densidad_exacta
export positive
export SG, apply_filter
export crear_matriz, norma_matrices, new_mc
export matriz_rutas, reduccion_transitiva
export lista_posets_3, lista_posets_4, lista_posets_5
export lista_posets_e_4, lista_posets_e_5, lista_posets_e_6 
export ranking_natural

end #module
