import numpy as np
from calculos import solve_wave_equation
import utils

# valores finitos para solução
dx = 0.1
dt = 0.01
t_max = 25
x_max = 20
y_max = 15
freq = 2
decay = 2

array_t, X, Y, array_wave = utils.create_wave(x_max, y_max, t_max, dx, dt, freq, decay)

# solução da edp
result = solve_wave_equation(
    array_wave, 
    [dx, dt], 
    t_max, 
    constant=1., 
    velocity_ratio=0.1
)

utils.plot_f_l_frames(result)
utils.plot_response(array_t, result, dt, dx)
#utils.animate_simulation(array_t, X, Y, result, 150)