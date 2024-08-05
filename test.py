from calculos import *
import environment_engine
import utils

for cond in [0, 160]:
    for c_squared in [80.2, 8020]:
        # valores finitos para solução
        dx = 0.01
        x_max = 0.5
        y_max = 0.5
        z_max = 0.5
        freq = 0.05
        t_max = 10 / freq
        dt = 0.0015909902497147803 * 10
        voltage = 1
        condu = cond
        print(f"started for conductivity {condu} and c² {c_squared}")
        X, Y, Z = utils.create_3d_space(x_max, y_max, z_max, dx)

        ew_format, ew = utils.make_coil(
            X, Y, Z, dx, 0.0033, 0.06, [0, 0, 0], 0, voltage
        )
        # ew_format2, ew2 = utils.make_coil(X, Y, Z, dx, 0.0033, 0.04, [0,0,0], 0, voltage)
        # ew_format = np.logical_or(ew_format.astype(bool), ew_format2.astype(bool))
        # ew = ew - ew2

        ew_format2 = []
        ew2 = []

        environ_params = [
            {"type": "base", "constant": 1},
            {
                "type": "solid_cilinder",
                "constant": 1,
                "radius": 0.05,
                "height": 0.7,
                "center": [0, 0, -0.35],
            },
        ]
        magnetic_permittivity = environment_engine.create_3d_environment(
            X, Y, Z, dx, environ_params
        )

        environ_params = [
            {"type": "base", "constant": 1},
            {
                "type": "solid_cilinder",
                "constant": c_squared,
                "radius": 0.05,
                "height": 0.7,
                "center": [0, 0, -0.35],
            },
        ]
        electric_permittivity = environment_engine.create_3d_environment(
            X, Y, Z, dx, environ_params
        )

        environ_params = [
            {"type": "base", "constant": 0},
            {
                "type": "solid_cilinder",
                "constant": condu,
                "radius": 0.05,
                "height": 0.7,
                "center": [0, 0, -0.35],
            },
        ]
        conductivity = environment_engine.create_3d_environment(
            X, Y, Z, dx, environ_params
        )

        environ_params = [
            {"type": "base", "constant": 10},
            {
                "type": "solid_cuboid",
                "constant": 0,
                "x_distance": 0.40,
                "y_distance": 0.40,
                "z_distance": 0.40,
                "center": [0.25, 0.25, 0.25],
            },
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
            dt=dt,
        )
        print("Simulation ended")

        # utils.plot_f_l_frames(result, X, Y, Z)

        response_params = {
            "radius": 0.06,
            "ring_radius": 0.0033,
            "center": [0, 0, -0.10],
        }

        values1 = utils.get_electromagnetic_response(
            array_t, X, Y, Z, result, dx, response_params, show_response=False
        )
        values1.to_csv(f"./results/obj-condu-{condu}-{freq}-{c_squared}.csv")

        response_params = {
            "radius": 0.06,
            "ring_radius": 0.0033,
            "center": [0, 0, 0.10],
        }

        values1 = utils.get_electromagnetic_response(
            array_t, X, Y, Z, result, dx, response_params, show_response=False
        )
        values1.to_csv(f"./results/free-condu-{condu}-{freq}-{c_squared}.csv")

        result = []
        # utils.plot_response(array_t, result, dt, dx, save_data=True)#, print_freqs=True)
        # utils.animate_simulation(array_t, result, 150, file_name="movie.mp4")
