import numpy as np
from calculos import solve_wave_equation
import environment_engine
import utils

# valores finitos para solução
dx = 0.1
dt = 0.05
t_max = 25
x_max = 40
y_max = 25
freq = 2
decay = 2

array_t, X, Y, array_wave = utils.create_wave(x_max, y_max, t_max, dx, dt, freq, decay)
environ_params = [
    {
        "type": "base",
        "constant": 1.
    },
    {
        "type": "borehole",
        "constant": 0.7,
        "x_distance": 2
    },
    {
        "type": "solid_rectangle",
        "center": [20, 15],
        "constant": 0.3,
        "x_distance": 100,
        "y_distance": 5
    },
    {
        "type": "solid_circle",
        "center": [20, 15],
        "constant": 0.8,
        "radius": 4,
        "x_pos": X,
        "y_pos": Y
    },
    {
        "type": "solid_circle",
        "center": [15, 20],
        "constant": 0.5,
        "radius": 4,
        "x_pos": X,
        "y_pos": Y
    },
    {
        "type": "permeable",
        "center": [27, 12],
        "constant": 0.3,
        "x_distance": 100,
        "y_distance": 5,
        "matrix_constant": 1.,
        "liquid_constant": 0.3,
        "percent": 0.3
    },
    {
        "type": "stripes",
        "center": [7, 30],
        "constant": 0.3,
        "x_distance": 100,
        "y_distance": 10,
        "1_constant": 1.,
        "2_constant": 0.3,
        "height": 0.5,
    }
]
constant = environment_engine.create_environment(array_wave, dx, environ_params, True)

# solução da edp
result = solve_wave_equation(
    array_wave, 
    [dx, dt], 
    t_max, 
    environment=constant
)

utils.plot_f_l_frames(result)
utils.plot_response(array_t, result, dt, dx)
#utils.animate_simulation(array_t, result, 150, file_name="movie.mp4")