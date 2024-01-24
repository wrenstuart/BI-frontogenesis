# Define sets of Richardson numbers (Ri) and isopycnal slopes (α) over which to iterate
# These numbers are defined by their base-10 logarithms
log_Ris = [3]
log_ss = [4]

include("loop_sims.jl")
loop_sims(log_Ris, log_ss)