module Posets

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
                  new_mc,
                  matriz_rutas, reduccion_transitiva,
                  lista_posets_3, lista_posets_4, lista_posets_5,
                  lista_posets_e_4, lista_posets_e_5, lista_posets_e_6 ,
                  ranking_natural, matriz_union_rankings, matriz_interseccion_rankings,
                  rankings_random, determinar_minimo_6, determinar_minimo_5,
                  crear_matriz, condorcet,
                  resaltados, determinar_tipo,
                  determinar_minimos,
                  derecha_abajo

export energia_por_pasos_p, energia_adhoc, energia_local
export energia_no_local, energia_trucada, energia_por_pasos
export energia_por_pasos_p
export monte_carlo, wang_landau, Simulacion, densidad_exacta
export positive
export SG, apply_filter
export new_mc
export matriz_rutas, reduccion_transitiva
export lista_posets_3, lista_posets_4, lista_posets_5
export lista_posets_e_4, lista_posets_e_5, lista_posets_e_6 
export ranking_natural, matriz_union_rankings, matriz_interseccion_rankings
export rankings_random, determinar_minimo_6, determinar_minimo_5
export crear_matriz
export resaltados, determinar_tipo
export condorcet
export determinar_minimos
export derecha_abajo

end #module
