export Σcount, fentropia, fs

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
