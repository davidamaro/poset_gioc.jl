import LightGraphs: SimpleDiGraph, nv, has_edge, inneighbors, outneighbors, edges
export resaltados, determinar_tipo

function resaltados(g1::SimpleDiGraph{T},g2::SimpleDiGraph{T}) where T <: Integer
    aristas_g2 = zeros(T, length(edges(g2)))
    for (i,ed) in enumerate(edges(g2))
        if has_edge(g1, ed)
            aristas_g2[i] = 2
        else
            aristas_g2[i] = 1
        end
    end
    aristas_g2
end

function determinar_tipo(grafica::SimpleDiGraph{T}) where T <: Integer
    len = nv(grafica)
    tipos = zeros(Int, len)
    for i in 1:len
        x = (inneighbors(grafica,i) |> length) > 0
        y = (outneighbors(grafica,i) |> length) > 0
        if x && y
            tipos[i] = 2
        elseif x & !y
            tipos[i] = 3
        else
            tipos[i] = 1
        end
    end
    tipos
end
