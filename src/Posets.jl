module Posets

using DelimitedFiles
using SimplePosets
using Statistics
using Random
using Combinatorics
using Markdown
using Primes

#Base.include(Main, "Generic.jl")
include("Generic.jl")

import .Generic:  energia_por_pasos_p, energia_adhoc, energia_adecuada,
                  energia_no_local, energia_trucada, energia_por_pasos,
                  energia_por_pasos_p,
                  monte_carlo, wang_landau, Simulacion, densidad_exacta,
                  positive,
                  SG, apply_filter,
                  crear_matriz, norma_matrices, new_mc
export energia_por_pasos_p, energia_adhoc, energia_adecuada
export energia_no_local, energia_trucada, energia_por_pasos
export energia_por_pasos_p
export monte_carlo, wang_landau, Simulacion, densidad_exacta
export positive
export SG, apply_filter
export crear_matriz, norma_matrices, new_mc

end #module
