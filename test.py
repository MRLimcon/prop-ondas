from calculos import *
import environment_engine
import utils

# valores finitos para solução
dx = 0.25
t_max = 20
x_max = 8
y_max = 8
z_max = 8
freq = 2

print("started")
X, Y, Z = utils.create_3d_space(x_max, y_max, z_max, dx)

environ_params = [
    {
        "type": "base",
        "constant": 1
    },
    {
        "type": "solid_cilinder",
        "constant": 1.2,
        "radius": 1,
        "height": 4,
        "center": [0, 0, 0]
    }
]
magnetic_permittivity = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

environ_params = [
    {
        "type": "base",
        "constant": 1
    },
    {
        "type": "solid_cilinder",
        "constant": 1.2,
        "radius": 1,
        "height": 4,
        "center": [0, 0, 0]
    }
]
electric_permittivity = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

environ_params = [
    {
        "type": "base",
        "constant": 0.75
    },
    {
        "type": "solid_cilinder",
        "constant": 0,
        "radius": 1,
        "height": 4,
        "center": [0, 0, 0]
    }
]
conductivity = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

environ_params = [
    {
        "type": "base",
        "constant": 0.75
    },
    {
        "type": "solid_cuboid",
        "constant": 0,
        "x_distance": 4,
        "y_distance": 4,
        "z_distance": 4,
        "center": [4, 4, 4]
    }
]
pml = environment_engine.create_3d_environment(X, Y, Z, dx, environ_params)

shape = list(electric_permittivity.shape)
shape.append(3)
ew = np.zeros(shape)
ew_format = np.zeros(shape, dtype=bool)
ew[int(shape[0]/2), int(shape[1]/2), int(shape[2]/2), 1] = 1
ew_format[int(shape[0]/2), int(shape[1]/2), int(shape[2]/2), 1] = True


print("Starting simulation")

dt, array_t, result = solve_electromagnetic_equation(
    excited_wave=ew,
    excited_wave_format=ew_format,
    dx=dx,
    mag_permi=magnetic_permittivity,
    conductivity=conductivity,
    elec_permi=electric_permittivity,
    max_time=t_max,
    freq=1,
    pml_layer=pml,
    visu_steps=150,
    dt=0.005
)
print("Simulation ended")

utils.plot_f_l_frames(result, X, Y, Z)
#utils.plot_response(array_t, result, dt, dx, save_data=True)#, print_freqs=True)
# utils.animate_simulation(array_t, result, 150, file_name="movie.mp4")