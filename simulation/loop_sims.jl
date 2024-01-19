# A function to loop over a set of values of Ri and α and run simulations for each one

include("../QOL.jl")
include("sim.jl")

function loop_sims(log_Ris, log_αs, resolution)
    for i₁ = 1 : length(log_Ris)
        for i₂ = 1 : length(log_αs)
            log_Ri = log_Ris[i₁]
            log_α = log_αs[i₂]
            Ri = 10 ^ log_Ri
            α = 10 ^ log_α
            label = get_label(log_Ri, log_α)
            run_sim(Ri, α, label, resolution)
        end
    end
end