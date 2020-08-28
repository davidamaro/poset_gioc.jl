# Please, keep this ref: 
# https://pixorblog.wordpress.com/2016/07/13/savitzky-golay-filters-julia/
using LinearAlgebra
using DelimitedFiles
using DSP

export SG, apply_filter

function vandermonde(halfWindow::Int, polyDeg::Int,T::Type=Float64)
    
    @assert halfWindow>=0
    @assert polyDeg>=0
    
    x=T[i for i in -halfWindow:halfWindow]

    n = polyDeg+1
    m = length(x)
    
    #V = Array{T}(m, n)
    V = zeros(T,(m, n))
    
    for i = 1:m
        V[i,1] = T(1)
    end
    for j = 2:n
        for i = 1:m
            V[i,j] = x[i] * V[i,j-1]
        end
    end

    return V
end

#________________________________________________________________

function SG_bug(halfWindow::Int, polyDeg::Int,T::Type=Float64)

    @assert 2*halfWindow>polyDeg
    
    println("good")
    V=vandermonde(halfWindow,polyDeg,T)
    println("good?")
    Q,R=qr(V)
    Q,R,V
end
function SG(halfWindow::Int, polyDeg::Int,T::Type=Float64)

    @assert 2*halfWindow>polyDeg
    
    V=vandermonde(halfWindow,polyDeg,T)
    Q,R=qr(V)
    SG=Base.inv(R)*Q'

    for i in 1:size(SG,1)
        SG[i,:]*=factorial(i-1)
    end
    
# CAVEAT: returns the transposed matrix

    return SG'
end

#________________________________________________________________

#function apply_filter{T}(filter::StridedVector{T},signal::StridedVector{T})
function apply_filter(filter::StridedVector{T},signal::StridedVector{T}) where T <: Real

    @assert isodd(length(filter))

    halfWindow = round(Int,(length(filter)-1)/2)
    
    padded_signal = 
	    [signal[1]*ones(halfWindow);
         signal;
         signal[end]*ones(halfWindow)]

    filter_cross_signal = conv(filter[end:-1:1], padded_signal)

    return filter_cross_signal[2*halfWindow+1:end-2*halfWindow]
end

##________________________________________________________________
#
#using Winston
#
#s=readdlm("signal.txt")[:,1]
#
#sg=SG(20,3) # halt-window, polynomal degree
#
##________________
#
#smoothed=apply_filter(sg[:,1],s)
#
#plot(s,"r")
#oplot(smoothed)
#title("Smoothed")
#savefig("smoothed.png")
#
##________________
#
#smoothed_d1=apply_filter(sg[:,2],s)
#
#plot(smoothed_d1)
#title("Smoothed derivative")
#savefig("smoothed_d1.png")
#
##________________
#
#smoothed_d2=apply_filter(sg[:,3],s)
#
#plot(smoothed_d2)
#title("Smoothed 2-derivative")
#savefig("smoothed_d2.png")
