import numpy as np
from calculos import solve_wave_equation
import environment_engine
import utils

# valores finitos para solução
dx = 0.1
dt = 0.01
t_max = 20
x_max = 40
y_max = 15
freq = 2
decay = 2

array_t, X, Y, array_wave = utils.create_wave(x_max, y_max, t_max, dx, dt, freq, decay)
environ_params = {
    "base": {
        "constant": 1.
    },
    "borehole": {
        "center": [15, 20],
        "constant": 0.5,
        "x_distance": 0.5,
        "y_distance": 20
    }
}
constant = environment_engine.create_environment(array_wave, dx, environ_params)

# solução da edp
result = solve_wave_equation(
    array_wave, 
    [dx, dt], 
    t_max, 
    environment=constant
)

utils.plot_f_l_frames(result)
utils.plot_response(array_t, result, dt, dx)
utils.animate_simulation(array_t, result, 150, file_name="movie.mp4")