from Classes.Experiment import Experiment
from Classes.Material import PhaseChangeMaterial as Pcm
from Solvers.TDMA_Solver import tdma
from Tools import console_informations_tool as cit
import numpy as np


class Simulation:
    def __init__(self, experiment: Experiment):
        # Create Arrays that will store values during simulation
        self.experiment = experiment
        self.temp_distribution = np.empty(experiment.num_nodes, dtype=np.float64)
        self.sensible_enthalpy_distribution = np.empty(experiment.num_nodes, dtype=np.float64)
        self.full_enthalpy_distribution = np.empty(experiment.num_nodes,
                                                   dtype=np.float64)  # Total enthalpies of cells: H = h + r.L.f
        self.melt_fraction_distribution = np.empty(experiment.num_nodes, dtype=np.float64)
        self.frontier_evolution = []
        self.materials_melt_distance_evolution = []
        self.outside_temperature_evolution = np.empty(experiment.time_steps, dtype=np.float64)
        self.wall_inside_temperature_evolution = np.empty(experiment.time_steps, dtype=np.float64)
        self.temp_tracker = np.empty(experiment.time_steps, dtype=np.float64)
        self.in_flux = np.empty(experiment.time_steps, dtype=np.float64)

    def simulate(self):
        run_simulation(self.experiment, self.temp_distribution, self.sensible_enthalpy_distribution,
                       self.full_enthalpy_distribution, self.melt_fraction_distribution, self.frontier_evolution,
                       self.materials_melt_distance_evolution,
                       self.wall_inside_temperature_evolution, self.outside_temperature_evolution, self.in_flux, self.temp_tracker)


def run_simulation(experiment: Experiment,
                   temp_distribution, sensible_enthalpy_distribution, full_enthalpy_distribution,
                   melt_fraction_distribution, frontier_evolution, materials_melt_distance_evolution, wall_inside_temperature_evolution,
                   outside_temperature_evolution, in_flux, temp_tracker):
    initialize_experiment(temp_distribution, sensible_enthalpy_distribution, melt_fraction_distribution,
                          experiment)

    # Loop over iterations until end of experiment
    for t in range(experiment.time_steps):

        if t % 1000 == 0 or t == experiment.time_steps - 1:
            cit.loading_bar(t, experiment.time_steps - 1)

        update_conditions(t * experiment.dt, experiment)
        outside_temperature_evolution[t] = experiment.outside_temperature

        calculate_step_time(prepare_diagonals(experiment, melt_fraction_distribution), experiment,
                            sensible_enthalpy_distribution, melt_fraction_distribution, frontier_evolution, materials_melt_distance_evolution)
        # Save full enthalpy
        calculate_full_enthalpy(full_enthalpy_distribution, sensible_enthalpy_distribution, melt_fraction_distribution,
                                experiment)

        calculate_all_cells_temperature(temp_distribution, sensible_enthalpy_distribution, experiment)
        wall_inside_temperature_evolution[t] = temp_distribution[-1]
        in_flux[t] = (calculate_in_flux(experiment, temp_distribution[-1], melt_fraction_distribution[-1]))
        temp_tracker[t] = temp_distribution[60]


def prepare_diagonals(experiment, melt_fraction_distribution):
    def melt_parameter(parameter: str, index: int):
        material = experiment.wall.material_repartition[index]
        if parameter == 'alpha':
            p_l = material.alpha_l
            p_s = material.alpha_s
        elif parameter == 'lambda':
            p_l = material.lambda_l
            p_s = material.lambda_s
        elif parameter == 'rho':
            p_l = material.rho_l
            p_s = material.rho_s
        elif parameter == 'c':
            p_l = material.c_l
            p_s = material.c_s
        else:
            print('Parameter Not Defined')
            p_l = material.c_l
            p_s = material.c_s
        return melt_fraction_distribution[index] * p_l + (1 - melt_fraction_distribution[index]) * p_s

    def fusion_temperature(index: int):
        return experiment.wall.material_repartition[index].fusion_temp

    # Prepare diagonals for matrix
    # Create empty diagonals
    left_diagonal = np.empty(experiment.num_nodes - 1)
    right_diagonal = np.empty(experiment.num_nodes - 1)
    big_diagonal = np.empty(experiment.num_nodes)
    # Create the constant boundary matrix: useful for left boundary as it's driven by temperature
    constant_boundary_matrix = np.zeros(experiment.num_nodes)
    # Fill diagonals
    k = 0
    # Loop over composites in the wall
    for composite, cells_num_composite in enumerate(experiment.wall.cells_number_by_layer):
        for i in range(k, k + cells_num_composite):
            gamma = melt_parameter('alpha', i) * experiment.r
            if i == 0:
                # Fill left BC
                nu = (2 * melt_parameter('lambda', i) - experiment.h_ext * experiment.dx) / (
                        2 * melt_parameter('lambda', i) + experiment.h_ext * experiment.dx)
                rho, c = melt_parameter('rho', i), melt_parameter('c', i)
                big_diagonal[i] = 1 + gamma * (2 - nu)
                right_diagonal[i] = -gamma
                constant_boundary_matrix[0] = gamma * (1 - nu) * rho * c * (
                        experiment.outside_temperature - fusion_temperature(i))
            elif i == experiment.num_nodes - 1:
                # Fill right interface with room
                phi = (2 * melt_parameter('lambda', i) - experiment.h_int * experiment.dx) / (
                        2 * melt_parameter('lambda', i) + experiment.h_int * experiment.dx)
                rho, c = melt_parameter('rho', i), melt_parameter('c', i)
                left_diagonal[i - 1] = -gamma
                big_diagonal[i] = 1 + gamma * (2 - phi)
                constant_boundary_matrix[i] = gamma * (1 - phi) * rho * c * (
                        experiment.inside_temperature - fusion_temperature(i))
            elif i == k:
                # Fill interface of composites (POV: Composite from right)
                l_1, l_2 = melt_parameter('lambda', i - 1), melt_parameter('lambda', i)
                rho, c = melt_parameter('rho', i), melt_parameter('c', i)
                p = (rho * c) / (
                        melt_parameter('c', i - 1) * melt_parameter('rho', i - 1))
                beta = (l_2 - l_1) / (l_1 + l_2)
                left_diagonal[i - 1] = -gamma * p * (1 - beta)
                big_diagonal[i] = 1 + gamma * (2 - beta)
                right_diagonal[i] = -gamma
                constant_boundary_matrix[i] = gamma * (1 - beta) * (
                        fusion_temperature(i - 1) - fusion_temperature(i)) * rho * c
            elif i == k + cells_num_composite - 1:
                # Fill interface of composites (POV: Composite from left)
                l_1, l_2 = melt_parameter('lambda', i), melt_parameter('lambda', i + 1)
                rho, c = melt_parameter('rho', i), melt_parameter('c', i)
                p = (rho * c) / (
                        melt_parameter('c', i + 1) * melt_parameter('rho', i + 1))
                beta = (l_1 - l_2) / (l_1 + l_2)
                left_diagonal[i - 1] = -gamma
                big_diagonal[i] = 1 + gamma * (2 - beta)
                right_diagonal[i] = -gamma * p * (1 - beta)
                constant_boundary_matrix[i] = gamma * (1 - beta) * (
                        fusion_temperature(i + 1) - fusion_temperature(i)) * rho * c
            else:
                # Fill common part: everything but boundaries
                left_diagonal[i - 1] = -gamma
                big_diagonal[i] = 1 + 2 * gamma
                right_diagonal[i] = -gamma
        k += cells_num_composite

    return left_diagonal, big_diagonal, right_diagonal, constant_boundary_matrix


def calculate_step_time(matrices: tuple, experiment: Experiment, sensible_enthalpy_distribution,
                        melt_fraction_distribution, frontier_evolution, materials_melt_distance_evolution):
    # Set the (k=0)th fraction to initial fraction (i.e. the one calculated at the end of previous step_time)
    new_enthalpy_distribution = sensible_enthalpy_distribution.copy()
    new_melt_fraction_vector = melt_fraction_distribution.copy()
    # Iterate over one step time until convergence
    max_i = 0
    while max_i < 5:
        max_i += 1

        # Prepare coefficients in the case of a phase-changing node
        change_coefficients(matrices, new_melt_fraction_vector)

        # Calculate new enthalpy
        delta_fraction_vector = new_melt_fraction_vector - melt_fraction_distribution
        new_enthalpy_distribution_save = new_enthalpy_distribution.copy()  # Save it to check convergence
        new_enthalpy_distribution = calculate_enthalpy_iteration(matrices, delta_fraction_vector,
                                                                 sensible_enthalpy_distribution, experiment)

        # Recalculate liquid fraction
        new_melt_fraction_vector = calculate_fraction_iteration(new_enthalpy_distribution,
                                                                new_melt_fraction_vector, matrices, experiment)
        if has_converged(new_enthalpy_distribution_save, new_enthalpy_distribution, experiment.tolerance):
            read_frontier(new_melt_fraction_vector, frontier_evolution, materials_melt_distance_evolution, experiment)
            break

    # Save new fraction and new enthalpy distribution
    melt_fraction_distribution[:] = np.round(new_melt_fraction_vector,
                                             6)  # Round melt fraction to avoid cases like 1.02e-19, etc...
    sensible_enthalpy_distribution[:] = new_enthalpy_distribution


def calculate_enthalpy_iteration(matrices: tuple, delta_fraction_vector: np.array, sensible_enthalpy_distribution,
                                 experiment: Experiment) -> np.array:
    # Prepare "d" array for TDMA solver
    wall = experiment.wall
    d_array = np.empty(experiment.num_nodes)
    for i in range(experiment.num_nodes):
        d_array[i] = sensible_enthalpy_distribution[i] - wall.material_repartition[i].n * delta_fraction_vector[i] \
                     + matrices[3][i]

    # Calculate new enthalpy
    new_enthalpy_distribution = tdma(matrices[0], matrices[1], matrices[2], d_array)

    return new_enthalpy_distribution


def calculate_fraction_iteration(new_enthalpy_distribution_vector: np.array,
                                 melt_fraction_vector: np.array, matrices: tuple, experiment: Experiment) -> np.array:
    # Calculate new melt fraction
    new_melt_fraction_vector = np.zeros(experiment.num_nodes)
    for i in range(0, experiment.num_nodes):
        material = experiment.wall.material_repartition[i]
        if material.latent_heat != 0:  # Check that it is considered as a phase changing material
            new_melt_fraction_vector[i] = melt_fraction_vector[i] + (
                    matrices[1][i] * new_enthalpy_distribution_vector[i]) / (
                                                  material.rho_l * material.latent_heat)
            new_melt_fraction_vector[i] = clamp(new_melt_fraction_vector[i], 0, 1)
    return new_melt_fraction_vector


def change_coefficients(matrices, new_melt_fraction_vector: np.array):
    for i in range(new_melt_fraction_vector.shape[0]):
        # Change ap coefficient
        if 0 < new_melt_fraction_vector[i] < 1:
            matrices[1][i] = 10 ** 150


def initialize_experiment(temp_distribution, sensible_enthalpy_distribution, melt_fraction_distribution, experiment):
    # Temperature distribution
    for i in range(experiment.num_nodes):
        temp_distribution[i] = experiment.wall_initial_temperature

    # Melt fraction distribution
    for i in range(0, experiment.num_nodes):
        if experiment.wall.material_repartition[i].fusion_temp < temp_distribution[i]:
            melt_fraction_distribution[i] = 1.0
        else:
            melt_fraction_distribution[i] = 0.0

    # Enthalpy distribution
    calculate_all_cells_enthalpy_temp(temp_distribution, sensible_enthalpy_distribution, experiment)


def calculate_all_cells_enthalpy_temp(temp_distribution, sensible_enthalpy_distribution, experiment):
    for i in range(experiment.num_nodes):
        sensible_enthalpy_distribution[i] = calculate_cell_enthalpy_temp(temp_distribution[i],
                                                                         experiment.wall.material_repartition[i])


def calculate_cell_enthalpy_temp(temperature: float, material: Pcm) -> float:
    # Determine enthalpy based on temperature
    if temperature < material.fusion_temp:
        enthalpy = (temperature - material.fusion_temp) * material.rho_s * material.c_s
    elif temperature == material.fusion_temp:
        enthalpy = 0
    else:
        enthalpy = (temperature - material.fusion_temp) * material.rho_l * material.c_l

    return enthalpy


def calculate_full_enthalpy(full_enthalpy_distribution, sensible_enthalpy_distribution, melt_fraction_distribution,
                            experiment: Experiment):
    for i in range(experiment.num_nodes):
        material = experiment.wall.material_repartition[i]
        if sensible_enthalpy_distribution[i] < 0:
            full_enthalpy_distribution[i] = sensible_enthalpy_distribution[i]
        elif sensible_enthalpy_distribution[i] == 0:
            full_enthalpy_distribution[i] = material.rho_l * material.latent_heat * melt_fraction_distribution[i]
        else:
            full_enthalpy_distribution[i] = sensible_enthalpy_distribution[i] + material.rho_l * material.latent_heat


def calculate_all_cells_temperature(temp_distribution, sensible_enthalpy_distribution, experiment):
    for i in range(experiment.num_nodes):
        temp_distribution[i] = calculate_cell_temperature(sensible_enthalpy_distribution[i],
                                                          experiment.wall.material_repartition[i])


def calculate_cell_temperature(enthalpy: float, material: Pcm) -> float:
    if enthalpy > 0:
        return material.fusion_temp + enthalpy / (material.rho_l * material.c_l)
    else:
        return material.fusion_temp + enthalpy / (material.rho_s * material.c_s)


def calculate_in_flux(experiment: Experiment, last_node_temperature: float, last_node_melt_fraction: float) -> float:
    lamb = last_node_melt_fraction * experiment.wall.material_repartition[-1].lambda_l + \
           (1 - last_node_melt_fraction) * experiment.wall.material_repartition[-1].lambda_s
    phi = (2 * lamb - experiment.dx * experiment.h_int) / (2 * lamb + experiment.dx * experiment.h_int)
    return experiment.h_int * (1 + phi) * (last_node_temperature - experiment.inside_temperature) / 2


def update_conditions(t: float, experiment: Experiment):
    experiment.outside_temperature = experiment.outside_temperature_evolution(t)


def clamp(number: float, smallest, largest): return max(smallest, min(number, largest))


def has_converged(previous_iteration_enthalpy: np.array, new_iteration_enthalpy: np.array, tolerance) -> bool:
    enthalpy_difference = np.sum(np.abs(previous_iteration_enthalpy-new_iteration_enthalpy))
    return -tolerance < enthalpy_difference < tolerance


def read_frontier(melt_fraction, frontier_evolution_array: list, materials_melt_distance_evolution: list,
                  experiment: Experiment):
    frontier_evolution_array.append([])
    for i in range(experiment.num_nodes - 1):
        if 0 < melt_fraction[i] < 1:
            if i < experiment.num_nodes - 1 and melt_fraction[i + 1] == 1:
                frontier_evolution_array[-1].append((i + 1 - melt_fraction[i]) * experiment.dx)
            else:
                frontier_evolution_array[-1].append((i + melt_fraction[i]) * experiment.dx)
        elif (melt_fraction[i] == 1 and melt_fraction[i + 1] == 0) or (
                melt_fraction[i] == 0 and melt_fraction[i + 1] == 1):
            frontier_evolution_array[-1].append((i + 1) * experiment.dx)

    materials_melt_distance_evolution.append([])  # [[timestep 1], [[material 1], ..., [[start phase node, distance, phase: O solid], ...]]]
    # Fill distances of phases for each material
    k = 0
    # Loop over composites in the wall
    total_length = 0.0
    for composite, cells_num_composite in enumerate(experiment.wall.cells_number_by_layer):
        materials_melt_distance_evolution[-1].append([])  # Add element (for current step time) corresponding to material iterated
        for i in range(k, k + cells_num_composite):
            distance = [0]
            phase = [0]
            add_phase = False
            if i == k:
                if (melt_fraction[i] == int(melt_fraction[i])) and (melt_fraction[i + 1] == int(melt_fraction[i + 1])) and (
                        melt_fraction[i] == 1 - melt_fraction[i + 1]):
                    distance = [experiment.dx]
                    phase = [int(melt_fraction[i])]
                    add_phase = True
                # First node of the material: check if melting
                elif 0 < melt_fraction[i] < 1:
                    # by assuming the phase on the right is the same that the i+1 node: (if i+1 node melting, arbitrary solid)
                    distance = [experiment.dx * (
                            int(melt_fraction[i + 1]) * (1 - melt_fraction[i]) + (1 - int(melt_fraction[i + 1])) * melt_fraction[i])]
                    phase = [1 - int(melt_fraction[i + 1])]
                    add_phase = True
            elif i == k + cells_num_composite - 1:
                add_phase = True
                # Last node of the material: check if melting
                if 0 < melt_fraction[i] < 1:
                    # by assuming the phase on the right is the same that the i-1 node: (if i+1 node melting, arbitrary solid)
                    phase = [int(melt_fraction[i - 1]), 1 - int(melt_fraction[i - 1])]
                    distance = [experiment.dx * ((1 - phase[0]) * (1 - melt_fraction[i]) + phase[0] * melt_fraction[i]),
                                experiment.dx * (phase[0] * (1 - melt_fraction[i]) + (1 - phase[0]) * melt_fraction[i])]
                else:
                    if not materials_melt_distance_evolution[-1][-1]:
                        distance = [experiment.dx * (i + 1 - k)]
                    else:
                        distance = [experiment.dx * (i + 1) - total_length]
                    phase = [int(melt_fraction[i])]
            elif 0 < melt_fraction[i] < 1:
                add_phase = True
                # If phase is only in the node and not neighbours nodes
                if melt_fraction[i + 1] == melt_fraction[i - 1]:
                    phase = [1 - int(melt_fraction[i + 1])]
                    # Check if it's the first frontier in material
                    if not materials_melt_distance_evolution[-1][-1]:
                        phase = [1 - phase[0], phase[0]]
                        distance = [experiment.dx * (i - k),
                                    experiment.dx * ((phase[0]) * (1 - melt_fraction[i]) + (1 - phase[0]) * melt_fraction[i])]
                    else:
                        distance = [experiment.dx * (1 - phase[0]) * (1 - melt_fraction[i]) + (phase[0]) * melt_fraction[i]]
                # else distance is well-defined
                else:
                    if melt_fraction[i - 1] == 1 or melt_fraction[i + 1] == 0:
                        phase = [1]

                    # Check if it's the first frontier in material
                    if not materials_melt_distance_evolution[-1][-1]:
                        distance = [experiment.dx * (i + ((1 - phase[0]) * (1 - melt_fraction[i]) + phase[0] * melt_fraction[i]) - k)]
                    else:
                        distance = [experiment.dx * (i + ((1 - phase[0]) * (1 - melt_fraction[i]) + phase[0] * melt_fraction[i])) - total_length]
            elif (melt_fraction[i] == int(melt_fraction[i])) and (melt_fraction[i + 1] == int(melt_fraction[i + 1])) and (
                    melt_fraction[i] == 1 - melt_fraction[i + 1]):
                # Check if it's the first frontier in material
                phase = [int(melt_fraction[i])]
                if not materials_melt_distance_evolution[-1][-1]:
                    distance = [experiment.dx * (i + ((1 - phase[0]) * (1 - melt_fraction[i]) + phase[0] * melt_fraction[i]) - k)]
                else:
                    distance = [experiment.dx * (i + ((1 - phase[0]) * (1 - melt_fraction[i]) + phase[0] * melt_fraction[i])) - total_length]
                add_phase = True
            # Add distance of phase to array
            if add_phase:
                for p in range(len(phase)):
                    total_length += distance[p]
                    materials_melt_distance_evolution[-1][-1].append([i, distance[p], phase[p], total_length])

        k += cells_num_composite
