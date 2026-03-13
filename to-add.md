# Post-processing

- For aggregated drifter timeseries, consider moving from per-drifter time intervals to a fixed time interval across all drifters
- Also consider adding variance alongside the aggregated mean for the data

# Simulation

- Get checkpointing working and document the change in Oceananigans version
- Submit git issue about why it wasn't working
- Then get pickup working
- Look into whether I can make this all a whole lot more efficient using AbstractFields rather than Fields
- Start outputting xz profiles