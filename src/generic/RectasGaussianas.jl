import Distributions: MvNormal
export generapuntuaciones_gaussian
export comparativaruidosa
export matrizposet

@doc Markdown.doc"""
`function randomsphere_point(n)`
  1. `n` es la dimension de la esfera en donde se toma el punto.
  Devuelve un punto en el octante generalizado.
# Ejemplo
```
julia > randomsphere_point(2)
punto en S1.
```
"""
function randomsphere_point(n)
    vec = [randn() for _ in 1:n]
    while any(vec.<0.0)
        vec = [randn() for _ in 1:n]
    end
    #@show vec
    vec/sqrt(dot(vec,vec))
end

@doc Markdown.doc"""
`function proyectar(vector,vectorrecta)`
Devuelve el vector resultado de proyectar `vector` sobre 
`vectorrecta`.
# Ejemplo
```
julia > sinejemplo
```
"""
function proyectar(vector,vectorrecta)
    (dot(vector,vectorrecta)/dot(vectorrecta,vectorrecta))*vectorrecta
end

@doc Markdown.doc"""
`function proporcion(vector,vectorrecta)`
Devuelve el valor de la proyeccion de `vector` sobre 
`vectorrecta`.
# Ejemplo
```
julia > sinejemplo
```
"""
function proporcion(vector,vectorrecta)
    (proyectar(vector, vectorrecta)./vectorrecta)[1]
end

@doc Markdown.doc"""
`function normalizacion(lista)`
  Toma una lista de valores con un rango arbitrario
  y los devuelve __normalizados__ en un rango 
  de $0$ a $1$.
# Ejemplo
```
julia > normalizacion([1,2,3])
[0.0, 0.5, 1.0]
```
"""
function normalizacion(lista)
    x,y = extrema(lista)

    [1-(y - u)/(y-x) for u in lista]
end

@doc Markdown.doc"""
`function minicomparativa(p,q)`
  compara dos puntos, `p` y `p`. Si $p_i \leq q_i$
  entonces $pRq$, por lo que se devuelve $1$.
  En otro caso, se devuelve $0$.
# Ejemplo
```
julia > minicomparativa([1,2,3],[0,1,2])
1
julia > minicomparativa([1,2,3],[0,1,4])
0
```
"""
function minicomparativa(p,q)
    if all(p .<= q)
        return 1
    else
        return 0
    end
end

@doc Markdown.doc"""
`function matrizposet(lista)`
  calcula el _poset de puntos_, es decir, recibe una lista
  de puntos que compara coordenada a coordenada. Si todos son mayores
  o iguales, se establece un enlace.
# Ejemplo
```
julia > matrizposet()
```
"""
function matrizposet(lista, metodo)
    dim, cantidadpuntos = size(lista)
    mat = zeros(Int64, cantidadpuntos,cantidadpuntos)
    for i in 1:cantidadpuntos, j in 1:cantidadpuntos
        if i == j
            continue
        end
        mat[i,j] = metodo(lista[:,i], lista[:,j])
    end
    mat
end

@doc Markdown.doc"""
`function generapuntuaciones_gaussian(numerorankings, numeronodos, dim)`
  Calcula el _poset de puntos_ y las puntuaciones.
  1. `numerorankings` es el numero de rankings que se generara. Denota el 
    numero de rectas usadas.
  2. `numeronodos` numero de nodos a usar. Corresponde al numero de
    puntos gaussianos que se usan.
  3. `dim` es la dimension de los puntos a usar.
# Ejemplo
```
julia > generapuntuaciones_gaussian(numerorankings, numeronodos, dim)
```
"""
function generapuntuaciones_gaussian(numerorankings, numeronodos, dim;
                                     safe::Bool = false,ruido::Bool = false, matnodos = 1, matruido = 1)

    listapuntos = [randomsphere_point(dim) for _ in 1:numerorankings];

    puntosnodos = rand(MvNormal([0 for _ in 1:dim],matnodos),numeronodos);

    if !ruido
      posetdepuntos = (matrizposet(puntosnodos, minicomparativa))
    else
      metodo = (x,y) -> comparativaruidosa(x,y,matruido=matruido)
      posetdepuntos = (matrizposet(puntosnodos, metodo))
    end

    if safe
      if !(posetdepuntos |> isacyclic)
        throw(error("no es un poset"))
      end
    end

    bloquenormalizado = [[proporcion(puntosnodos[:,i], j) for i in 1:numeronodos] |> normalizacion for j in listapuntos];

    posetdepuntos, hcat(bloquenormalizado...)
end

function comparativaruidosa(a::Float64, b::Float64)
    x,y = [randn() for _ in 1:2]
    a + x <= b + y
end
function comparativaruidosa(a::Array{Float64,1}, b::Array{Float64,1}; matruido = 1)
    n = length(a)
    @assert n == length(b)

    x,y = rand(MvNormal([0 for _ in 1:n],matruido),2)
    all((a .+ x) .<= (b .+ y))
end
