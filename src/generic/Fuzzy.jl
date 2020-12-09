export Σcount, fentropia, fs
export αcut_poset, permpuntuaciones, m3, m2

function finter(fp1, fp2)
    @assert length(fp1) == length(fp2)
    [min(fp1[i], fp2[i]) for i in 1:length(fp1)]
end

function funion(fp1, fp2)
    @assert length(fp1) == length(fp2)
    [max(fp1[i], fp2[i]) for i in 1:length(fp1)]
end

finverso(fp) = [1-fp[i] for i in 1:length(fp)]

function fs(fp1, fp2)
    sum(finter(fp1, fp2))/sum(fp1)
end

Σcount(fp1) = sum(fp1)

fentropia(x) = fs(funion(x,finverso(x)), finter(x, finverso(x)))

function αcut_poset(poset, α)
    copia = copy(poset)
    a = filter(x -> x > α, copia)
    n,_ = size(a)
    b = zeros(Int, n,n)
    for i in 1:n, j in i+1:n
        if a[i,j] > 0 && a[j,i] == 0
            b[i,j] = 1
        elseif a[j,i] > 0 && a[i,j] == 0
            b[j,i] = 1
        end
    end
    b
end

function permpuntuaciones(mat, k)
    mmat = copy( mat )
    n,m = mmat |> size
    #@show n,m
    for i in 1:m
        for _ in 1:k
            p = rand(1:n-1)
            #@show p
            #mmat[:, i][[p, p+1]] .= mat[:, i][[p+1, p]]
            t = mmat[p, i]
            mmat[p,i] = mmat[p+1,i]
            mmat[p+1,i] = t
        end

    end
    mmat
end

function m3(mat; alter::Bool = false)
    n,m = mat |> size
    output = zeros(Float64,n,3)

    for i in 1:n
        x,y = extrema(mat[i,:])
        m = median(mat[i,:])
        output[i,1] = y
        output[i,2] = m
        output[i,3] = x
    end
     output
end

function m2(mat; alter::Bool = false)
    n,m = mat |> size
    output = zeros(Float64,n,2)

    for i in 1:n
        x,y = extrema(mat[i,:])
        output[i,1] = y
        output[i,2] = x
    end
     output
end
