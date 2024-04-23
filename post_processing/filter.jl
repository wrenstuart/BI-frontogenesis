using CairoMakie
using FFTW

f(x, y) = cos(x + 2y) - 3sin(5x) + 2cos(10y)
Δ = 0.1
x = [x for x in -5π:Δ:5π]
y = [y for y in 0:Δ:6π]
M = length(x)
N = length(y)
z = [f(x, y) for x in x, y in y]

L = 5   # Filter lengthscale

# Outline of how wavenumbers work is below
#=
k_min = 2π/(x[end] - x[1])
l_min = 2π/(y[end] - y[1])
k = [k_min * (m-1) for m in 1:M]
l = [l_min * (n-1) for n in 1:N]
# Last values of above arrays should be 2π/Δ ^
k_cut = 2π/L
=#

m_cut = Int(round((x[end] - x[1])/L))
n_cut = Int(round((y[end] - y[1])/L))

#=z_f = fft(z)
z_f[1+m_cut:M-m_cut, :] .= 0
z_f[:, 1+n_cut:N-n_cut] .= 0
z_smooth = real(ifft(z_f))=#

function gaussian_filter_2d(z, m_cut, n_cut)
    (M, N) = size(z)
    z_f = fft(z)
    for i in 1:M, j in 1:N
        m = minimum([mod(i-1,M), mod(1-i,M)])
        n = minimum([mod(j-1,N), mod(1-j,N)])
        @info exp(-(m^2/m_cut^2 + n^2/n_cut^2))
        z_f[i, j] *= exp(-(m^2/m_cut^2 + n^2/n_cut^2))
    end
    real(ifft(z_f))
end

heatmap(x, y, gaussian_filter_2d(z, m_cut, n_cut))
#heatmap(x, y, z)