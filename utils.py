from matplotlib import cm
import fortran_bins.utils as utils
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def generate_excited_wave(t_max: float, dt: float, freq: float, type: str = "ricker", simu_type: str = "normal", show_env: bool = False)-> np.ndarray:
    if type == "ricker":
        t_vals = np.arange(-(1/freq)+dt, (1/freq), dt)
        exponent = np.exp(-(np.pi**2)*(freq**2)*(t_vals**2))
        result = (1-(2*(np.pi**2)*(freq**2)*(t_vals**2)))*exponent
        #ints = utils.utils.fix_ew(ew_len=len(result), ew=result)
        #t_vals = t_vals[ints[0]-1: ints[1]-1]
        #result = result[ints[0]-1: ints[1]-1]

    else:
        t_vals = np.arange(0, t_max, dt)
        result = np.sin(2*np.pi*freq*t_vals)

    if show_env:
        plt.plot(t_vals, result)
        plt.xlabel("Tempo (s)")
        plt.ylabel("Valor da oscilação")
        plt.show()

    if simu_type == "acoustic":
        final_result = result 
        result = np.zeros([final_result.shape[0], 2])
        result[:, 0] = final_result
        result[:, 1] = final_result

    return result

def generate_mu_lambda(shear_speed, pressure_speed, plot_env: bool = False): #, rho):
    mu = shear_speed**2 #*rho
    lambda_1 = (pressure_speed**2) - (2*mu) #*rho

    if plot_env:
        shw = plt.imshow(lambda_1.T, cmap = cm.coolwarm)
        plt.title("Lambda")
        plt.colorbar(shw, cmap = cm.coolwarm)
        plt.show()

        shw = plt.imshow(mu.T, cmap = cm.coolwarm)
        plt.title("Mu")
        plt.colorbar(shw, cmap = cm.coolwarm)
        plt.show()

    return lambda_1, mu


def make_fft(data: np.ndarray, timestep: float, dist: float) -> None:
    """
        Visualizes the frequency peaks for the values
    """
    new_data = (data-np.average(data))/np.max(data)
    starts = 0

    for i, value in enumerate(new_data):
        if all(np.abs(val) >= 0.3 for val in new_data[i:i+3]):
            starts = i
            break

    values = new_data[starts:]

    n1 = len(values)
    fft_results1 = np.abs(np.fft.fft(values))
    fft_results1 = fft_results1 / np.max(fft_results1)
    freq = np.fft.fftfreq(n1, d=timestep)

    peaks1 = find_peaks(fft_results1)[0]
    filters = [fft_results1[peak] > 0.3 and freq[peak] >= 0 for peak in peaks1]
    peaks1 = peaks1[filters]

    print(f"Starts at {starts*timestep} s and distance is {dist}")
    print(f"Frequencies are {freq[peaks1]}")

def detect_signals(data: np.ndarray):
    new_data = (data-np.average(data))/np.max(data)
    len_data = len(data)
    starts = []
    ends = []

    for i, value in enumerate(new_data):
        if all(np.abs(val) >= 0.3 for val in new_data[i:i+3]):
            if len(starts) <= len(ends):
                starts.append(i)
        else:
            if len(ends) < len(starts):
                ending_integers = len_data - i - 1
                if ending_integers > 15 and all(np.abs(value) < 0.2 for value in new_data[i:i+10]):
                    ends.append(i)

        if i == len_data-1 and len(ends) < len(starts):
            ends.append(i)

    if len(starts) >= 2 and len(ends) >= 2:
        return starts[:2], ends[:2]
    else:
        return None

def create_wave(
        x_max: float, y_max: float, t_max: float, dx: float, dt: float) -> tuple[np.ndarray]:
    """
        Create the initial conditions for the simulation,
        x_max is the height of the simulated borehole,
        y_max is the radius of the simulation,
        t_max is the maximum time for the simulation,
        dx and dt are the steps in space and time,
        freq is the frequency of the generated monopole wave,
        decay is the spatial decay of the generated wave
    """
    array_x = np.arange(-x_max/2, (x_max/2)+dx, dx)
    array_y = np.arange(-y_max, y_max+dx, dx)
    array_t = np.arange(0, t_max, dt)

    X, Y = np.meshgrid(array_x, array_y)
    array_wave = np.zeros(X.shape)
    return array_t, X, Y, array_wave
    
def plot_f_l_frames(array: np.ndarray) -> None:
    """
        Plot the initial condition and the last frame of the simulation
    """
    if len(array.shape) == 3:
        shw = plt.imshow(array[0], cmap = cm.coolwarm)
        plt.title("First frame")
        plt.colorbar(shw, cmap = cm.coolwarm)
        plt.show()

        shw = plt.imshow(array[-1], cmap = cm.coolwarm)
        plt.title("Last frame")
        plt.colorbar(shw, cmap = cm.coolwarm)
        plt.show()
    elif len(array.shape) == 4:
        shw = plt.imshow(array[0, :, :, 0], cmap = cm.coolwarm)
        plt.title("First frame - u")
        plt.colorbar(shw, cmap = cm.coolwarm)
        plt.show()

        shw = plt.imshow(array[-1, :, :, 0], cmap = cm.coolwarm)
        plt.title("Last frame - u")
        plt.colorbar(shw, cmap = cm.coolwarm)
        plt.show()

        shw = plt.imshow(array[0, :, :, 1], cmap = cm.coolwarm)
        plt.title("First frame - v")
        plt.colorbar(shw, cmap = cm.coolwarm)
        plt.show()

        shw = plt.imshow(array[-1, :, :, 1], cmap = cm.coolwarm)
        plt.title("Last frame - v")
        plt.colorbar(shw, cmap = cm.coolwarm)
        plt.show()

def plot_response(array_t: np.ndarray, array: np.ndarray, dt: float, 
    dx: float, distance:float = 0.12, print_freqs: bool = False) -> None:
    """
        Generate the receiver responses along the borehole
    """
    y_pos = int(array.shape[2]/2)+2
    steps = int(distance/dx)
    initial = int(array.shape[1]/2 - int(3.5/dx)) 
    starts = []
    ends = []
    signal_data = []
    dists = []

    
    if len(array.shape) == 3:
        fig, axes = plt.subplots(nrows=4, ncols=1)
        for i in range(4):
            x_pos = initial - (steps*i)
            dist = np.abs(x_pos - int(array.shape[1]/2))*dx
            new_values = detect_signals(array[:, x_pos, y_pos])

            if new_values != None:
                starts.append(new_values[0])
                ends.append(new_values[1])
                signal_data.append(array[:, x_pos, y_pos])
                dists.append(dist)

            if print_freqs:
                make_fft(array[:, x_pos, y_pos], dt, dist)

            axes[i].plot(array_t, array[:, x_pos, y_pos])
            axes[i].set_title(f"Distance from source: {dist} m")
        
        fig.tight_layout()
        plt.show()
    elif len(array.shape) == 4:
        for j in [0, 1]:
            fig, axes = plt.subplots(nrows=4, ncols=1)
            for i in range(4):
                x_pos = initial - (steps*i)
                dist = np.abs(x_pos - int(array.shape[1]/2))*dx
                new_values = detect_signals(array[:, x_pos, y_pos, j])

                if new_values != None:
                    starts.append(new_values[0])
                    ends.append(new_values[1])
                    signal_data.append(array[:, x_pos, y_pos, j])
                    dists.append(dist)

                if print_freqs:
                    make_fft(array[:, x_pos, y_pos, j], dt, dist)

                axes[i].plot(array_t, array[:, x_pos, y_pos, j])
                axes[i].set_title(f"Distance from source: {dist} m")

            fig.tight_layout()
            plt.show()
        
    #starts = np.array(starts)
    #ends = np.array(ends)
    #signal_data = np.array(signal_data)
    #dists = np.array(dists)

    #print("First wave slowness:")
    #print(np.average(np.abs((starts[:-1, 0]-starts[-1, 0])*dt / (dists[:-1]-dists[-1]))))
    #print("Stoneley wave slowness:")
    #print(np.average(np.abs((starts[:-1, 1]-starts[-1, 1])*dt / (dists[:-1]-dists[-1]))))

    #print("Timestamps:")
    #print(starts*dt)
    #print(ends*dt)

def animate_simulation(
        array_t: np.ndarray, result: np.ndarray, num_frames: float, file_name: str = None, 
        X: np.ndarray = None, Y: np.ndarray = None, ani_type: str = "2d") -> None:
    """
        Create an animation of the wave propagation
    """
    global ax, kind

    if ani_type == "surface":
        kind = True
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plot = ax.plot_surface(X, Y, result[0].T, cmap=cm.coolwarm,
                        linewidth=0)
    else:
        kind = False
        fig, ax = plt.subplots()
        plot = ax.imshow(result[0], cmap = cm.coolwarm)
    

    def change_plot(frame_number, zarray, plot):
        global ax, kind

        print(f"frame n: {frame_number}")
        ax.clear()
        if kind:
            plot = ax.plot_surface(X, Y, zarray[frame_number, :, :].T,
                cmap=cm.coolwarm, linewidth=0)
        else:
            plot = ax.imshow(zarray[frame_number, :, :], cmap = cm.coolwarm)
        return plot,

    fps = 30
    frame_step = int(len(array_t)/num_frames)
    plot_array = np.zeros([num_frames, result.shape[1], result.shape[2]])

    for i in range(num_frames):
        plot_array[i] = result[i*frame_step]

    ani = animation.FuncAnimation(fig, change_plot, frames=num_frames, fargs=(plot_array, plot), interval=1000 / fps, blit=False)

    if file_name != None:
        ani.save(f"./results/{file_name}")
    else:
        plt.show()
