using Plots

∇x2x2_F2(x1, x2) = 1/4 * exp(-x2^2/2) * (x2^4 - 2 * x1 * x2^3 + 6*x1 * x2 - 5*x2^2 + 2)
∇x2x1_F2(x1, x2) = -1/2 * exp(-x2^2/2) + x2^2/2 * exp(-x2^2/2)
Sx1x2(x1,x2) = ∇x2x1_F2(x1, x2) / (-1*∇x2x2_F2(x1, x2))
∇x2_F1(x1,x2) = 2*x2

Dx1_F1(x1, x2) = -x1 + Sx1x2(x1,x2) * ∇x2_F1(x1,x2)

∇x2_F2(x1, x2) = (x2/2 - x1/2) * exp(-x2^2/2) + (x2^2/4 - (x1 * x2)/2) * (-x2 * exp(-x2^2/2))

x1_kp1(x1, x2; τ = 1/4) = x1 - τ * Dx1_F1(x1, x2)

x2_kp1(x1, x2; τ = 1/4) = x2 - τ * ∇x2_F2(x1, x2) + Sx1x2(x1, x2) * (-τ * Dx1_F1(x1, x2))


function main(x1_0 = 1.0, x2_0 = 1.0, itermax = 40; τ=1/4)
    x1s = zeros(itermax); x2s = zeros(itermax)
    x1s[1] = x1_0; x2s[1] = x2_0
    for k = 2:itermax
        x1k = x1s[k-1]; x2k = x2s[k-1]
        x1s[k] = x1_kp1(x1k, x2k; τ=τ)
        x2s[k] = x2_kp1(x1k, x2k; τ=τ)
    end
    return x1s, x2s
end

x1s, x2s = main(0.6,0.6)  # 1,1 does not converge !? nor 0.7,0.7 and higher tenths places
plot(x1s, x2s, seriestype= :scatter)

x1s, x2s = main(0.7,0.7, 100; τ=0.01) # still diverges
plot!(x1s, x2s, seriestype= :scatter)
