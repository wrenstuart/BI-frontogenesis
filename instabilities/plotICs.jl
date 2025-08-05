using CairoMakie
include("modes.jl")
include("../QOL.jl")

Ri = 5
s = 1e4

p = get_scales(Ri, s)
Lx = 2 * 2π * p.L * 0.4^0.5
Ly = Lx
uᵢ, vᵢ, wᵢ, bᵢ = generate_ic(Ri, Lx, p.U, N = 50)

correctPeriodicIndex(i::Int64, n::Int64) :: Int64 = (n + (i-1)) % n + 1
function δx(perArr::Matrix{Float64}, i::Int64, j::Int64)
    n, ~ = size(perArr)
    return perArr[correctPeriodicIndex(i+1, n), j] - perArr[correctPeriodicIndex(i-1, n), j]
end
function δy(perArr::Matrix{Float64}, i::Int64, j::Int64)
    ~, m = size(perArr)
    return perArr[i, correctPeriodicIndex(j+1, m)] - perArr[i, correctPeriodicIndex(j-1, m)]
end
∂x(perArr::Matrix{Float64}, i, j) = δx(perArr::Matrix{Float64}, i, j) / Δx
∂y(perArr::Matrix{Float64}, i, j) = δy(perArr::Matrix{Float64}, i, j) / Δy

Δx = Lx/100
Δy = Ly/100
x = 0 : Δx : Lx
y = 0 : Δy : Ly
bSurf = reshape([uᵢ(x, y, 0) for y in y for x in x], 101, 101)
uSurf = reshape([uᵢ(x, y, 0) for y in y for x in x], 101, 101)
vSurf = reshape([vᵢ(x, y, 0) for y in y for x in x], 101, 101)
wSurf = reshape([wᵢ(x, y, 0) for y in y for x in x], 101, 101)
ζSurf = reshape([∂x(vSurf, i, j) - ∂y(uSurf, i, j) for j in 1 : 101 for i in 1 : 101], 101, 101)
δSurf = reshape([∂x(uSurf, i, j) + ∂y(vSurf, i, j) for j in 1 : 101 for i in 1 : 101], 101, 101)

fig = Figure()
bAx = Axis(fig[1, 1])
ζAx = Axis(fig[1, 2])
δAx = Axis(fig[1, 3])
uAx = Axis(fig[2, 1])
vAx = Axis(fig[2, 2])
wAx = Axis(fig[2, 3])
heatmap!(bAx, x, y, bSurf)
heatmap!(ζAx, x, y, ζSurf)
heatmap!(δAx, x, y, δSurf)
heatmap!(uAx, x, y, uSurf)
heatmap!(vAx, x, y, vSurf)
heatmap!(wAx, x, y, wSurf)
display(fig)

#b̅ = [sum([bSurf[i,j] for i = 1 : 101]) for j = 1 : 101]
#plot(y, b̅)