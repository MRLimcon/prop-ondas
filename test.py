import numpy as np
from calculos import *
import environment_engine
import matplotlib.pyplot as plt
import utils

# valores finitos para solução
dx = 0.02
dt = 0.00000125
t_max = 0.004
x_max = 20
y_max = 10
freq = 8000

print("started")
array_t, X, Y, array_wave = utils.create_wave(x_max, y_max, t_max, dx, dt)
print("generated wave")
excited_wave = utils.generate_excited_wave(t_max, dt, dx, freq, type="ricker", simu_type = "acoustic")
print("generated excited wave")
# em metros/20us
environ_params = [
    {
        "type": "base",
        "constant": 2000
    },
    {
        "type": "borehole",
        "constant": 0.0,
        "x_distance": 0.11
    }
]
shear_speed, bore_params = environment_engine.create_environment(array_wave, dx, environ_params, True)

environ_params = [
    {
        "type": "base",
        "constant": 3500
    },
    {
        "type": "borehole",
        "constant": 1500,
        "x_distance": 0.11
    }
]

pressure_speed = environment_engine.create_environment(array_wave, dx, environ_params, True)

lambda_1, mu = utils.generate_mu_lambda(shear_speed, pressure_speed)#, True)

print("Starting simulation")

#result = solve_wave_equation(
#    array_wave, 
#    excited_wave,
#    [dx, dt], 
#    t_max, 
#    environment=constant
#)

dt, array_t, result = solve_elastodynamic_equation(
    initial_condition=array_wave,
    excited_wave=excited_wave,
    steps=[dx, dt],
    mu = mu,
    lambda_1=lambda_1,
    rho=shear_speed,
    max_time=t_max,
    borehole_params=bore_params,
    visu_steps=500
)
print("Simulation ended")

utils.plot_f_l_frames(result)
utils.plot_response(array_t, result, dt, dx, save_data=True)#, print_freqs=True)
# utils.animate_simulation(array_t, result, 150, file_name="movie.mp4")