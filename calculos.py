import fortran_bins.calculations as calculations
import numpy as np

# função para solução da edp de onda
def solve_wave_equation(
    initial_condition: np.ndarray,
    steps: tuple[float],
    max_time: float,
    constant: float = 1.,
    velocity_ratio: float = 1.
) -> np.ndarray:
    """
        Solution to 1d wave equation, initial_condition is an 1d array of floats,
        initial_velocities is an 1d array of floats,
        constant is the speed constant (wave speed, light speed, etc),
        steps is an tuple of the format (dx, dt, ), remember to always put dt < 100*c²/dx²,
        max_time is the max time for the solution, remember to scale with dt to allocate enough memory,
        type is the type of the solution, "free" is the solution with two free edges, "closed" is the solution with two unmovable edges
    """
    length = int(max_time/steps[1])
    ic_lenx = initial_condition.shape[0]
    ic_leny = initial_condition.shape[1]

    array = calculations.calc.free_wave_equation_2d(
        ic=initial_condition,
        c=constant,
        sol_len=length,
        delta_x=steps[0],
        delta_t=steps[1],
        ic_lenx=ic_lenx,
        ic_leny=ic_leny,
        vr = velocity_ratio
    )
    return array.T