export monte_carlo, wang_landau, densidad_exacta, positive
export Simulacion
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

