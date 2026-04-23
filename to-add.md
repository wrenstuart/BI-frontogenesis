# Post-processing

- For aggregated drifter timeseries, consider moving from per-drifter time intervals to a fixed time interval across all drifters
- omega equation, pv
- Introduce a way to look select a section of phase space within a certain time window and isolate all drifters which pass through this region
- A little Fourier transform to estimate the scales of the fronts

# Simulation

- Submit git issue about why checkpointing wasn't working
- Also submit git issue about interpolation of complicated fields onto drifters
    - This is the reason that I can't use AbstractFields to save on memory usage
- See if can reduce viscosity, either without changing anything or by increasing the amplitudes of less unstable modes

# For the paper

- Add some stuff from 21st April meeting to here