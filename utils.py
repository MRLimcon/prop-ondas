from matplotlib import cm
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def make_fft(values, timestep):
    n1 = len(values)
    fft_results1 = abs(np.fft.fft(values))
    fft_results1 = fft_results1 / np.max(fft_results1)
    freq = np.fft.fftfreq(n1, d=timestep)

    # PEGANDO OS PICOS DE FREQUÊNCIA
    peaks1 = find_peaks(fft_results1)[0]
    filters = [fft_results1[peak] > 0.3 and freq[peak] >= 0 for peak in peaks1]
    peaks1 = peaks1[filters]

    print(f"As frequências de pico para o sensor são {freq[peaks1]}")

def create_wave(x_max, y_max, t_max, dx, dt, freq, decay):
    # criação das condições iniciais
    array_x = np.arange(-x_max, x_max+dx, dx)
    array_y = np.arange(-y_max, y_max+dx, dx)
    array_t = np.arange(0, t_max, dt)

    X, Y = np.meshgrid(array_x, array_y)
    array_wave = np.cos(freq*np.sqrt(X**2 + Y**2))*np.exp(-((decay*X)**2))*np.exp(-((decay*Y)**2))
    return array_t, X, Y, array_wave
    
def plot_f_l_frames(array: np.ndarray):
    # plots
    plt.imshow(array[0], cmap = cm.coolwarm)
    plt.title("First frame")
    plt.show()

    plt.imshow(array[-1], cmap = cm.coolwarm)
    plt.title("Last frame")
    plt.show()

def plot_response(array_t, array, dt, dx):
    x_pos = int(array.shape[1]/2.5)
    y_pos = int(array.shape[2]/2)
    steps = int(array.shape[1]/80)
    initial = int(array.shape[1]/2.5)
    make_fft(array[:, x_pos, y_pos], dt)

    fig, axes = plt.subplots(nrows=5, ncols=1)
    for i in range(5):
        x_pos = initial + (steps*i)
        axes[i].plot(array_t, array[:, x_pos, y_pos])
        axes[i].set_title(f"distancia: {np.abs(x_pos - int(array.shape[1]/2))*dx}")
    #fig.legend()
    fig.tight_layout()
    plt.show()

def animate_simulation(array_t, X, Y, result, num_frames, ani_type = "2d"):
    global ax, kind

    if ani_type == "surface":
        kind = True
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        plot = ax.plot_surface(X, Y, result[0].T, cmap=cm.coolwarm,
                        linewidth=0)
    else:
        kind = False
        fig, ax = plt.subplots()
        plot = ax.imshow(result[0])
    

    def change_plot(frame_number, zarray, plot):
        global ax, kind

        print(f"frame n: {frame_number}")
        ax.clear()
        if kind:
            plot = ax.plot_surface(X, Y, zarray[frame_number, :, :].T,
                cmap=cm.coolwarm, linewidth=0)
        else:
            plot = ax.imshow(zarray[frame_number, :, :])
        return plot,

    fps = 30
    frame_step = int(len(array_t)/num_frames)
    plot_array = np.zeros([num_frames, result.shape[1], result.shape[2]])

    for i in range(num_frames):
        plot_array[i] = result[i*frame_step]

    ani = animation.FuncAnimation(fig, change_plot, frames=num_frames, fargs=(plot_array, plot), interval=1000 / fps, blit=False)

    ani.save("movie.mp4")
