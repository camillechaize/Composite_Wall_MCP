from Solvers.Composite_Wall_Solver import Simulation
from Plot import tools_plot as tp
from Plot import image_frontier_plot as ifp
import numpy as np
import matplotlib.pyplot as plt


def plot_simulation(simulation: Simulation):
    plt.style.use('dark_background')
    ratio = (19, 11)

    # Create "x" arrays
    time_array = np.linspace(0, (simulation.experiment.time_steps - 1) * simulation.experiment.dt,
                             simulation.experiment.time_steps)
    step = len(simulation.materials_melt_distance_evolution) // simulation.experiment.max_abscissa_values

    image_frontier = ifp.image_from_frontiers(simulation, ratio, step)

    temperatures_array = prepare_double_arrays(time_array, simulation.wall_inside_temperature_evolution,
                                               simulation.outside_temperature_evolution, step)
    flux_in_array = prepare_simple_arrays(time_array, simulation.in_flux, step)
    # Create two subplots and unpack the output array immediately
    f, axs = plt.subplots(2, 2, figsize=ratio)

    axs[0, 1].plot(flux_in_array[0], flux_in_array[1], linewidth=1,
                   label=r'$\phi_{in}$ in $W.m^{-2}$', color='sandybrown')
    axs[0, 1].set_title(f'Surface heat flux in room over ${tp.time_formatting(simulation.experiment.duration)}$')
    axs[0, 1].set_ylabel(f'Surface heat flux in $W.m^{-2}$')
    axs[0, 1].hlines(y=0, xmin=flux_in_array[0][0], xmax=flux_in_array[0][-1], color='white', linestyle='-')
    axs[0, 1].fill_between(flux_in_array[0], flux_in_array[1], where=(flux_in_array[1] > 0), color='cyan', alpha=0.3)
    axs[0, 1].fill_between(flux_in_array[0], flux_in_array[1], where=(flux_in_array[1] < 0), color='coral', alpha=0.3)
    axs[0, 1].legend()

    temp_tracker_array = prepare_simple_arrays(time_array, simulation.temp_tracker, step)
    axs[0, 0].plot(temperatures_array[0], temperatures_array[1] - 273.15, linewidth=1, label=r'$T_{w,int}$',
                   color='aqua')
    axs[0, 0].plot(temperatures_array[0], temperatures_array[2] - 273.15, linewidth=1, label=r'$T_{out}$', color='coral')
    axs[0, 0].plot(temp_tracker_array[0], temp_tracker_array[1] - 273.15, linewidth=1, label=r'$T_{middle}$', color='green')
    axs[0, 0].set_ylabel(f'Temperature in $°C$')
    axs[0, 0].set_title(f'Temperature evolution over ${tp.time_formatting(simulation.experiment.duration)}$')
    axs[0, 0].legend()

    scale_wall = tp.scale_formatting_unique(simulation.experiment.wall.length)
    wall_legend = tp.build_wall_legend(simulation.experiment.wall)
    axs[1, 0].imshow(image_frontier, extent=[0, simulation.experiment.duration, 0, simulation.experiment.wall.length * scale_wall[0]], aspect='auto')
    axs[1, 0].set_ylabel(f'Wall Depth in ${scale_wall[1]}$')
    axs[1, 0].set_xlabel(f'Time in $s$')
    axs[1, 0].set_title('Inside Wall Evolution')
    axs[1, 0].legend(wall_legend[0], wall_legend[1])
    axs[1, 1].imshow(image_frontier, extent=[0, simulation.experiment.duration, 0, simulation.experiment.wall.length * scale_wall[0]], aspect='auto')

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
