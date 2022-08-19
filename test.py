import numpy as np
from calculos import *
import environment_engine
import matplotlib.pyplot as plt
import utils

# valores finitos para solução
dx = 0.2
dt = 0.0001
t_max = 135
x_max = 20
y_max = 10
freq = 1

print("started")
array_t, X, Y, array_wave = utils.create_wave(x_max, y_max, t_max, dx, dt)
print("generated wave")
excited_wave = utils.generate_excited_wave(t_max, dt, freq, type="ricker", simu_type = "acoustic")
print("generated excited wave")
environ_params = [
    {
        "type": "base",
        "constant": 0.1
    },
    {
        "type": "borehole",
        "constant": 0.05,
        "x_distance": 1.# 0.1
    }
]
shear_speed = environment_engine.create_environment(array_wave, dx, environ_params)#, True)

environ_params = [
    {
        "type": "base",
        "constant": 0.1
    },
    {
        "type": "borehole",
        "constant": 0.05,
        "x_distance": 1.# 0.1
    }
]
"""
    {
        "type": "permeable",
        "center": [27, 12],
        "constant": 0.6,
        "x_distance": 100,
        "y_distance": 100,
        "matrix_constant": 1.,
        "liquid_constant": 0.5,
        "percent": 0.3
    },
    {
        "type": "solid_rectangle",
        "center": [20, 15],
        "constant": 0.4,
        "x_distance": 100,
        "y_distance": 5
    },
    #{
    #    "type": "permeable",
    #    "center": [27, 12],
    #    "constant": 0.3,
    #    "x_distance": 100,
    #    "y_distance": 100,
    #    "matrix_constant": 1.,
    #    "liquid_constant": 0.3,
    #    "percent": 0.3
    #},
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
        "constant": 0.7,
        "radius": 4,
        "x_pos": X,
        "y_pos": Y
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
"""
pressure_speed = environment_engine.create_environment(array_wave, dx, environ_params)#, True)

# solução da edp
#result = solve_wave_equation(
#    array_wave, 
#    excited_wave,
#    [dx, dt], 
#    t_max, 
#    environment=constant
#)

lambda_1, mu = utils.generate_mu_lambda(shear_speed, pressure_speed)

dt, array_t, result = solve_elastodynamic_equation(
    initial_condition=array_wave,
    excited_wave=excited_wave,
    steps=[dx, dt],
    mu = mu,
    lambda_1=lambda_1,
    rho=shear_speed,
    max_time=t_max
)

utils.plot_f_l_frames(result)
utils.plot_response(array_t, result, dt, dx, print_freqs=True)
# utils.animate_simulation(array_t, result, 150, file_name="movie.mp4")