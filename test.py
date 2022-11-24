from calculos import *
import environment_engine
import utils
import pip

# valores finitos para solução
dx = 2
x_max = 60
y_max = 60
z_max = 60
freq = 1*(10**11)
dt = 5e-16
print(dt)
t_max = 3/freq
print("started")
X, Y, Z = utils.create_3d_space(x_max, y_max, z_max, dx)

ew_format, ew = utils.make_coil(X, Y, Z, dx, 2, 12.5, [0,0,10], 0)
ew_format2, ew2 = utils.make_coil(X, Y, Z, dx, 2, 12.5, [0,0,0], 0)
ew_format = np.logical_or(ew_format.astype(bool), ew_format2.astype(bool))
ew = ew - ew2

ew_format2 = []
ew2 = []

environ_params = [
    {
        "type": "base",
        "constant": 0.000001256637061*100
    },
    {
        "type": "solid_cilinder",
        "constant": 6.3*(10**-3)*100,
        "radius": 10,
        "height": 40,
        "center": [0, 0, 0]
    }
]
magnetic_permittivity = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

environ_params = [
    {
        "type": "base",
        "constant": 80.2*8.85*(10**-12)*100
    },
    {
        "type": "solid_cilinder",
        "constant": 4.5*8.85*(10**-12)*100,
        "radius": 10,
        "height": 40,
        "center": [0, 0, 0]
    }
]
electric_permittivity = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

environ_params = [
    {
        "type": "base",
        "constant": 0
    },
    {
        "type": "solid_cilinder",
        "constant": 1/(96.1*100*(10**-9)),
        "radius": 10,
        "height": 40,
        "center": [0, 0, 0]
    }
]
conductivity = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

environ_params = [
    {
        "type": "base",
        "constant": 1
    },
    {
        "type": "solid_cuboid",
        "constant": 0,
        "x_distance": 40,
        "y_distance": 40,
        "z_distance": 40,
        "center": [30, 30, 30]
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
    "radius": 12.5,
    "ring_radius": 2,
    "center": [0, 0, -10]
}

utils.get_electromagnetic_response(array_t, X, Y, Z, result, dx, response_params, show_response=True)
#utils.plot_response(array_t, result, dt, dx, save_data=True)#, print_freqs=True)
# utils.animate_simulation(array_t, result, 150, file_name="movie.mp4")