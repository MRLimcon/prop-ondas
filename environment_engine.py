import numpy as np

def get_edges(dictionary: dict, steps: float) ->tuple:
    lenx = int(2*dictionary["x_distance"]/steps) 
    leny = int(dictionary["y_distance"]/steps)
    center = [int(dictionary["center"][0]/steps), int(dictionary["center"][1]/steps)]

    return lenx, leny, center

def correct_edges(lenx: float, leny: float, center: tuple[int], environment_shape: np.ndarray) -> list[int]:
    corrected_edges = np.zeros(4, dtype=np.int64)
    left_edge = int(center[0]-(lenx/2))
    right_edge = int(center[0]+(lenx/2))
    bottom_edge = int(center[1]-(leny/2))
    upper_edge = int(center[1]+(leny/2))

    if left_edge < 0:
        corrected_edges[0] = 0
    else:
        corrected_edges[0] = left_edge

    if right_edge > environment_shape[0]:
        corrected_edges[1] = environment_shape[0]
    else:
        corrected_edges[1] = right_edge

    if bottom_edge < 0:
        corrected_edges[2] = 0
    else:
        corrected_edges[2] = bottom_edge

    if upper_edge > environment_shape[1]:
        corrected_edges[3] = environment_shape[1]
    else:
        corrected_edges[3] = upper_edge

    return corrected_edges

def get_indexes(lenx: int, leny: int, steps: float) -> np.ndarray:
    return


def create_environment(init_array: np.ndarray, steps, params: dict[str, dict])-> np.ndarray:
    """
        Create an environment of velocity constants,
        the environment has the shape of init_array, with steps as spatial steps,
        params has the shapes and parameters for the intented area
    """
    environment = np.zeros(init_array.shape)

    base = params["base"]
    velocity_constant = base["constant"]
    environment[:,:] = velocity_constant

    for shape in params:
        shape_params = params[shape]

        if shape == "solid_square":
            lenx, leny, center = get_edges(shape_params, steps)
            velocity_constant = shape_params["constant"]

            end_shape = correct_edges(lenx, leny, center, environment.shape)
            environment[end_shape[0]:end_shape[1], end_shape[2]:end_shape[3]] = velocity_constant

        elif shape == "solid_circle":
            radius = shape_params["radius"]/steps
            velocity_constant = shape_params["constant"]
            environment[np.sqrt(((environment*steps)-center)**2) <= radius] = velocity_constant

        elif shape == "stripes":
            lenx, leny, center = get_edges(shape_params, steps)
            stripe1_constant = shape_params["1_constant"]
            stripe2_constant = shape_params["2_constant"]

        elif shape == "permeable":
            lenx, leny, center = get_edges(shape_params, steps)
            matrix_constant = shape_params["matrix_constant"]
            liquid_constant = shape_params["liquid_constant"]

        else:
            continue

    borehole_params = params["borehole"]
    lenx, leny, center = get_edges(borehole_params, steps)
    velocity_constant = borehole_params["constant"]

    end_shape = correct_edges(lenx, leny, center, environment.shape)
    environment[end_shape[0]:end_shape[1], :] = velocity_constant

    return environment