from calculos import *
import environment_engine
import utils

# valores finitos para solução
dx = 0.01
x_max = 0.6
y_max = 0.6
z_max = 0.6
freq = 0.05
t_max = 10/freq
dt = 0.0015909902497147803*20
voltage = 1
print("started")
X, Y, Z = utils.create_3d_space(x_max, y_max, z_max, dx)

ew_format, ew = utils.make_coil(X, Y, Z, dx, 0.0033, 0.04, [0,0,0.1], 0, voltage)
ew_format2, ew2 = utils.make_coil(X, Y, Z, dx, 0.0033, 0.04, [0,0,0], 0, voltage)
ew_format = np.logical_or(ew_format.astype(bool), ew_format2.astype(bool))
ew = ew - ew2

ew_format2 = []
ew2 = []

environ_params = [
    {
        "type": "base",
        "constant": 1
    },
    {
        "type": "solid_cilinder",
        "constant": 1,
        "radius": 0.05,
        "height": 0.4,
        "center": [0, 0, 0]
    }
]
magnetic_permittivity = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

environ_params = [
    {
        "type": "base",
        "constant": 4.5
    },
    {
        "type": "solid_cilinder",
        "constant": 80.2,
        "radius": 0.05,
        "height": 0.4,
        "center": [0, 0, 0]
    }
]
electric_permittivity = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

environ_params = [
    {
        "type": "base",
        "constant": 120
    },
    {
        "type": "solid_cilinder",
        "constant": 4.8,
        "radius": 0.05,
        "height": 0.4,
        "center": [0, 0, 0]
    }
]
conductivity = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

environ_params = [
    {
        "type": "base",
        "constant": 10
    },
    {
        "type": "solid_cuboid",
        "constant": 0,
        "x_distance": 0.40,
        "y_distance": 0.40,
        "z_distance": 0.40,
        "center": [0.30, 0.30, 0.30]
    }
]
pml = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)


print("Starting simulation")

dt, array_t, result = solve_electromagnetic_equation(
    excited_wave=ew,
    excited_wave_format=ew_format,
    dx=dx,
    mag_permi=magnetic_permittivity,
    conductivity=conductivity,
    elec_permi=electric_permittivity,
    max_time=t_max,
    freq=freq,
    pml_layer=pml,
    visu_steps=250,
    dt=dt
)
print("Simulation ended")

utils.plot_f_l_frames(result, X, Y, Z)

response_params = {
    "radius": 0.04,
    "ring_radius": 0.0033,
    "center": [0, 0, -0.10]
}

values1 = utils.get_electromagnetic_response(array_t, X, Y, Z, result, dx, response_params, show_response=True)
values1.to_csv(f"./results/voltage-{voltage}.csv")

#utils.plot_response(array_t, result, dt, dx, save_data=True)#, print_freqs=True)
# utils.animate_simulation(array_t, result, 150, file_name="movie.mp4")