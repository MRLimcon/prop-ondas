import fortran_bins.calculations as calculations
import numpy as np

# função para solução da edp de onda
def solve_wave_equation(
    initial_condition: np.ndarray,
    excited_wave: float,
    steps: tuple[float],
    max_time: float,
    environment: float = 1.
) -> np.ndarray:
    """
        Solution to 2d wave equation, initial_condition is an 2d array of floats,
        constant is the environment conditions for the speed constant,
        steps is an tuple of the format (dx, dt, ), remember to always put dt < c²/dx²,
        max_time is the max time for the solution, remember to scale with dt to allocate enough memory,
    """
    length = int(max_time/steps[1])
    ic_lenx = initial_condition.shape[0]
    ic_leny = initial_condition.shape[1]
    constant = (environment**2)/(steps[0]**2)

    array = calculations.calc.free_wave_equation_2d(
        ew=excited_wave,
        ew_len=len(excited_wave),
        c=constant,
        sol_len=length,
        dt=steps[1],
        lenx=ic_lenx,
        leny=ic_leny
    )
    array = array.T
    return array

def solve_elastodynamic_equation(
    initial_condition: np.ndarray,
    excited_wave: float,
    steps: tuple[float],
    borehole_params: tuple[int],
    mu,
    lambda_1,
    rho,
    max_time: float,
    visu_steps: int = 500
) -> np.ndarray:
    """
        Solution to 2d elastodynamic equation, initial_condition is an 2d array of floats,
        constant is the environment conditions for the speed constant,
        steps is an tuple of the format (dx, dt, ), remember to always put dt < c²/dx²,
        max_time is the max time for the solution, remember to scale with dt to allocate enough memory,
    """
    length = int(max_time/steps[1])
    ic_lenx = initial_condition.shape[0]
    ic_leny = initial_condition.shape[1]

    visu_dt = round(length/visu_steps)
    visu_steps = int(length/visu_dt)

    array = calculations.elastodynamic.elastodynamic_2d(
        lenx=ic_leny, 
        leny=ic_lenx, 
        sol_len=length, 
        ew=excited_wave, 
        ew_len=excited_wave.shape[0], 
        mu=mu.T, 
        ar_len=visu_steps,
        ar_steps=visu_dt,
        l=lambda_1.T, 
        rho=rho.T, 
        dt=steps[1],
        dx=steps[0],
        lower=borehole_params[0],
        upper=borehole_params[1]
    )

    times = np.array([visu_dt*steps[1]*i for i in range(array.shape[0])])
    return visu_dt*steps[1], times, array

def solve_electromagnetic_equation(
    excited_wave: np.ndarray,
    excited_wave_format: np.ndarray,
    dx: float,
    mag_permi: np.ndarray,
    conductivity: np.ndarray,
    elec_permi: np.ndarray,
    max_time: float,
    freq: float,
    pml_layer: np.ndarray = None,
    visu_steps: int = 500,
    dt: float = None
) -> np.ndarray:
    """
        Solution to 3d electromagnetic equation
    """
    if dt == None:
        dt = 0.577*dx/np.max(1/np.sqrt(mag_permi*elec_permi))

    length = int(max_time/dt)
    lenx = conductivity.shape[0]
    leny = conductivity.shape[1]
    lenz = conductivity.shape[2]

    visu_dt = round(length/visu_steps)
    visu_steps = int(length/visu_dt)

    if type(pml_layer) == type(None):
        pml_layer = np.zeros([lenx, leny, lenz])

    array = calculations.electromagnetic.electromagnetic_3d(
        lenx=lenx, 
        leny=leny, 
        lenz=lenz,
        sol_len=length, 
        ew=excited_wave,
        ew_format=excited_wave_format,
        ar_len=visu_steps,
        ar_steps=visu_dt,
        conductivity=conductivity, 
        elec_permi=elec_permi, 
        mag_permi=mag_permi,
        dt=dt,
        dx=dx,
        pml=pml_layer,
        freq=freq, 
    )

    times = np.array([visu_dt*dt*i for i in range(array.shape[0])])
    return visu_dt*dt, times, array