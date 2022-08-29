import numpy as np
import fortran_bins.utils as utils
import matplotlib.pyplot as plt
from matplotlib import cm

def get_edges(dictionary: dict, steps: float) ->tuple:
    lenx = int(2*dictionary["x_distance"]/steps)
    if "center" in dictionary:
        center = [int(dictionary["center"][0]/steps), int(dictionary["center"][1]/steps)]
        leny = int(dictionary["y_distance"]/steps)

        return lenx, leny, center
    elif "y_distance" in dictionary:
        leny = int(dictionary["y_distance"]/steps)

        return lenx, leny
    else:
        return lenx

def correct_edges(lenx: float, center: tuple[int], environment_shape: np.ndarray, leny: float = None) -> list[int]:
    if leny == None:
        leny = environment_shape[1]

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


def create_environment(init_array: np.ndarray, steps, params: list[dict], show_env: bool = False)-> np.ndarray:
    """
        Create an environment of velocity constants,
        the environment has the shape of init_array, with steps as spatial steps,
        params has the shapes and parameters for the intented area
    """
    environment = np.zeros(init_array.shape)

    base = next(param for param in params if param["type"] == "base")
    velocity_constant = base["constant"]
    environment[:,:] = velocity_constant

    for shape_params in params:
        shape = shape_params["type"]

        if shape == "solid_rectangle":
            lenx, leny, center = get_edges(shape_params, steps)
            velocity_constant = shape_params["constant"]

            end_shape = correct_edges(lenx, center, environment.shape, leny)

            if end_shape[0] == end_shape[1] or end_shape[2] == end_shape[3]:
                continue

            environment[end_shape[0]:end_shape[1], end_shape[2]:end_shape[3]] = velocity_constant

        elif shape == "solid_circle":
            radius = shape_params["radius"]
            center = [shape_params["center"][1], shape_params["center"][0]]
            velocity_constant = shape_params["constant"]
            X = shape_params["x_pos"]
            Y = shape_params["y_pos"]

            environment[
                utils.utils.make_circle(
                    lenx=environment.shape[0],
                    leny=environment.shape[1],
                    x=X-X[0, 0],
                    y=Y-Y[0, 0],
                    center=center,
                    radius=radius
                ).astype(dtype=np.bool)
            ] = velocity_constant

        elif shape == "stripes":
            lenx, leny, center = get_edges(shape_params, steps)
            end_shape = correct_edges(lenx, center, environment.shape, leny)
            height = int(shape_params["height"]/steps)
            stripe1_constant = shape_params["1_constant"]
            stripe2_constant = shape_params["2_constant"]
            layers = int(leny/height) + 1

            if end_shape[0] == end_shape[1] or end_shape[2] == end_shape[3]:
                continue
            
            environment[end_shape[0]:end_shape[1], 
                end_shape[2]:end_shape[3]] = utils.utils.make_stripes(
                lenx=end_shape[1]-end_shape[0],
                leny=end_shape[3]-end_shape[2],
                matrix_const=stripe1_constant,
                liquid_const=stripe2_constant,
                height=height,
                layers=layers
            )

        elif shape == "permeable":
            lenx, leny, center = get_edges(shape_params, steps)
            matrix_constant = shape_params["matrix_constant"]
            liquid_constant = shape_params["liquid_constant"]
            percent = shape_params["percent"]

            if percent >= 1.:
                percent = percent/100

            end_shape = correct_edges(lenx, center, environment.shape, leny)

            if end_shape[0] == end_shape[1] or end_shape[2] == end_shape[3]:
                continue
            
            environment[end_shape[0]:end_shape[1], 
                end_shape[2]:end_shape[3]] = utils.utils.make_permeable(
                lenx=end_shape[1]-end_shape[0],
                leny=end_shape[3]-end_shape[2],
                matrix_const=matrix_constant,
                liquid_const=liquid_constant,
                percent = percent
            )

        else:
            continue

    borehole_params = base = next(param for param in params if param["type"] == "borehole")
    lenx = get_edges(borehole_params, steps)
    center = [int(init_array.shape[0]/2), int(init_array.shape[1]/2)]
    velocity_constant = borehole_params["constant"]

    end_shape = correct_edges(lenx, center, environment.shape)
    environment[end_shape[0]:end_shape[1], :] = velocity_constant

    if show_env:
        shw = plt.imshow(environment.T, cmap = cm.coolwarm)
        plt.title("Velocity constant on the environment")
        bar = plt.colorbar(shw, cmap = cm.coolwarm)
        plt.tight_layout()
        plt.show()


    if np.abs(velocity_constant) < 0.0001:
        return environment, [end_shape[0], end_shape[1]]
    else:
        return environment