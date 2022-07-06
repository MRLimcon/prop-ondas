from matplotlib import cm
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def make_fft(values: np.ndarray, timestep: float) -> None:
    """
        Visualizes the frequency peaks for the values
    """
    n1 = len(values)
    fft_results1 = abs(np.fft.fft(values))
    fft_results1 = fft_results1 / np.max(fft_results1)
    freq = np.fft.fftfreq(n1, d=timestep)

    peaks1 = find_peaks(fft_results1)[0]
    filters = [fft_results1[peak] > 0.3 and freq[peak] >= 0 for peak in peaks1]
    peaks1 = peaks1[filters]

    print(f"As frequências de pico para o sensor são {freq[peaks1]}")

def create_wave(
        x_max: float, y_max: float, t_max: float, 
        dx: float, dt: float, freq: float, decay: float) -> tuple[np.ndarray]:
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
    array_wave = np.cos(freq*np.sqrt(X**2 + Y**2))*np.exp(-(decay*(X**2 + Y**2)))
    return array_t, X, Y, array_wave
    
def plot_f_l_frames(array: np.ndarray) -> None:
    """
        Plot the initial condition and the last frame of the simulation
    """
    shw = plt.imshow(array[0], cmap = cm.coolwarm)
    plt.title("First frame")
    plt.colorbar(shw, cmap = cm.coolwarm)
    plt.show()

    shw = plt.imshow(array[-1], cmap = cm.coolwarm)
    plt.title("Last frame")
    plt.colorbar(shw, cmap = cm.coolwarm)
    plt.show()

def plot_response(array_t: np.ndarray, array: np.ndarray, dt: float, dx: float) -> None:
    """
        Generate the receiver responses along the borehole
    """
    x_pos = int(array.shape[1]/2.5)
    y_pos = int(array.shape[2]/2)
    steps = int(array.shape[1]/80)
    initial = int(array.shape[1]/2.5)
    # make_fft(array[:, x_pos, y_pos], dt)

    fig, axes = plt.subplots(nrows=5, ncols=1)
    for i in range(5):
        make_fft(array[:, x_pos, y_pos], dt)
        x_pos = initial + (steps*i)
        axes[i].plot(array_t, array[:, x_pos, y_pos])
        axes[i].set_title(f"Distance from source: {np.abs(x_pos - int(array.shape[1]/2))*dx} m")
    fig.tight_layout()
    plt.show()

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
