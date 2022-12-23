from Solvers.Composite_Wall_Solver import Simulation
import numpy as np


def image_from_frontiers(simulation: Simulation, ratio, step):
    x_max = len(simulation.materials_melt_distance_evolution) // step
    y_max = int(x_max * (ratio[1] / ratio[0]))

    image_array = np.zeros((x_max, y_max, 4), dtype=np.float64)

    for t in range(0, x_max):
        k = 0
        for m in simulation.materials_melt_distance_evolution[t * step]:
            for p in m:
                num_alloc_pixels = round(p[1] / simulation.experiment.wall.length * y_max)
                color = simulation.experiment.wall.material_repartition[p[0]].solid_color
                if p[2] == 1:
                    color = simulation.experiment.wall.material_repartition[p[0]].liquid_color

                for i in range(k, num_alloc_pixels + k):
                    image_array[t][clamp(i, 0, y_max - 1)] = color

                k += num_alloc_pixels

    return np.rot90(image_array)


def clamp(value, v_min, v_max):
    return min(max(v_min, value), v_max)
