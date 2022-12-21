from Solvers.Composite_Wall_Solver import Simulation
from Plot import tools_plot as tp
import numpy as np
import matplotlib.pyplot as plt


def plot_simulation(simulation: Simulation):
    plt.style.use('dark_background')

    # Create "x" arrays
    time_array = np.linspace(0, (simulation.experiment.time_steps - 1) * simulation.experiment.dt,
                             simulation.experiment.time_steps)
    step = len(simulation.frontier_evolution) // simulation.experiment.max_abscissa_values
    frontiers_arrays = prepare_frontiers_arrays(time_array, simulation.frontier_evolution, step)
    temperatures_array = prepare_double_arrays(time_array, simulation.wall_inside_temperature_evolution,
                                               simulation.outside_temperature_evolution, step)
    flux_in_array = prepare_simple_arrays(time_array, simulation.in_flux, step)
    # Create two subplots and unpack the output array immediately
    f, axs = plt.subplots(2, 2, figsize=(19, 11))

    # If there is a phase changing material: plot evolution of frontier
    if frontiers_arrays.shape[0] != 0:
        format_scale = tp.scale_formatting_unique(max(frontiers_arrays[1]))
        axs[1, 0].scatter(frontiers_arrays[0], frontiers_arrays[1] * format_scale[0], marker=".")
        axs[1, 0].set_title('Frontier evolution')
        axs[1, 0].set_ylabel(format_scale[1])

    axs[0, 1].plot(flux_in_array[0], flux_in_array[1], marker='.',
                   label=r'$\phi_{in}$', color='sandybrown')
    axs[0, 1].set_title(f'Surface heat flux in room over {tp.time_formatting(simulation.experiment.duration)}')
    axs[0, 1].hlines(y=0, xmin=flux_in_array[0][0], xmax=flux_in_array[0][-1], color='white', linestyle='-')
    axs[0, 1].fill_between(flux_in_array[0], flux_in_array[1], where=(flux_in_array[1] > 0), color='cyan', alpha=0.3)
    axs[0, 1].fill_between(flux_in_array[0], flux_in_array[1], where=(flux_in_array[1] < 0), color='coral', alpha=0.3)
    axs[0, 1].legend()

    axs[0, 0].plot(temperatures_array[0], temperatures_array[1] - 273.15, marker=".", label=r'$T_{w,int}$',
                   color='aqua')
    axs[0, 0].plot(temperatures_array[0], temperatures_array[2] - 273.15, marker=".", label=r'$T_{out}$', color='coral')
    axs[0, 0].legend()

    # Set to full-screen
    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    plt.show()


def prepare_frontiers_arrays(time_array, frontiers_evolutions: list, pas: int):
    frontiers_tracked = []
    for i in range(0, len(frontiers_evolutions), pas):
        for n in range(len(frontiers_evolutions[i])):
            frontiers_tracked.append([time_array[i], frontiers_evolutions[i][n]])
    return np.transpose(np.array(frontiers_tracked))


def prepare_double_arrays(time_array, in_temp: np.array, out_temp: np.array, pas: int):
    values_kept = []
    for i in range(0, len(in_temp), pas):
        values_kept.append([time_array[i], in_temp[i], out_temp[i]])
    return np.transpose(np.array(values_kept))


def prepare_simple_arrays(time_array, values_array: np.array, pas: int):
    values_kept = []
    for i in range(0, len(values_array), pas):
        values_kept.append([time_array[i], values_array[i]])
    return np.transpose(np.array(values_kept))
