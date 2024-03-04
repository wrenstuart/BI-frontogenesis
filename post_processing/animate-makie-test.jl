using GLMakie

function animate()

    x = [x for x in 0:0.1:10]           # x domain
    t_end = 10                          # Final time
    n_frames = 100                      # Number of frames to animate
    time = range(0, t_end, n_frames)    # Collection of times to sample

    y = sin.(x)                         # Initial condition
    fig, ax, lineplot = lines(x, y, color = :black) # Plot the initial condition and save the axes
    
    record(fig, "test_ani.mp4", time, framerate = 20) do t
        # define anonymous function for t, which will take its values from time
        y = sin.(x .- t)        # Calculate the new curve to draw
        delete!(ax, lineplot)   # Delete the previous curve from the axes
        lineplot = lines!(ax, x, y, color = :black) # Draw the new one on
    end

end

animate()

function animate2()
    
    t = Observable(0.0)                 # This will be updated to change the whole plot

    x = [x for x in 0:0.1:10]           # x domain
    t_end = 10                          # Final time
    n_frames = 100                      # Number of frames to animate
    time = range(0, t_end, n_frames)    # Collection of times to sample

    y = @lift sin.(x .+ $t)             # Define observable y that updates when t does
    # f(t) = sin.(x .- t)               # Another way
    # y = lift(f, t)                    # of defining y
    
    fig = lines(x, y, color = :black)
    record(τ -> t[] = τ, fig, "test_ani.mp4", time, framerate = 25)

end

animate2()