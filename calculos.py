import fortran_bins.calculations as calculations
import numpy as np

# função para solução da edp de onda
def solve_wave_equation(
    initial_condition: np.ndarray,
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
        ic=initial_condition,
        c=constant,
        sol_len=length,
        dt=steps[1],
        ic_lenx=ic_lenx,
        ic_leny=ic_leny
    )
    return array.T