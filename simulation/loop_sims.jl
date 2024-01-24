# A function to loop over a set of values of Ri and s and run simulations for each one

include("../QOL.jl")
include("sim.jl")

function loop_sims(log_Ris, log_ss)
    for i₁ = 1 : length(log_Ris)
        for i₂ = 1 : length(log_ss)
            log_Ri = log_Ris[i₁]
            log_s = log_ss[i₂]
            Ri = 10 ^ log_Ri
            s = 10 ^ log_s
            label = get_label(log_Ri, log_s)
            run_sim(Ri, s, label)
        end
    end
end