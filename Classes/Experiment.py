from Classes.Wall import Wall
from Config import temperature_config
from Tools.ini_tool import read_ini


class Experiment:
    def __init__(self, experience_ini_path: str, wall_csv_path: str, materials_csv_path: str):
        # Get values from config file
        exp_settings = get_experience_settings(experience_ini_path)

        self.num_nodes = exp_settings.getint('Precision', 'num_nodes')
        self.dt = exp_settings.getfloat('Precision', 'delta_time')
        self.duration = exp_settings.getfloat('Precision', 'duration')
        self.tolerance = exp_settings.getfloat('Precision', 'tolerance')
        self.time_steps = int(self.duration // self.dt)

        self.wall_surface = exp_settings.getfloat('Wall', 'wall_surface')
        self.wall_initial_temperature = exp_settings.getfloat('Wall', 'wall_initial_temperature') + 273.15

        self.h_int = exp_settings.getfloat('Room', 'h_int')
        self.inside_temperature = exp_settings.getfloat('Room', 'inside_temperature') + 273.15

        self.h_ext = exp_settings.getfloat('Outside', 'h_ext')
        self.outside_temperature = exp_settings.getfloat('Outside', 'outside_temperature') + 273.15
        self.outside_temperature_evolution = temperature_config.update_outside_temperature

        self.wall = Wall(wall_csv_path, materials_csv_path)
        self.num_nodes = self.wall.calculate_cells_number_by_material(
            self.num_nodes)  # Update number of nodes, so it fits layers ratios
        self.dx = self.wall.length / self.num_nodes
        self.r = self.dt / (self.dx ** 2)

        self.max_abscissa_values = exp_settings.getint('Plot', 'max_abscissa_values')
        self.number_of_apartments = exp_settings.getint('Plot', 'number_of_apartments')
        self.nuclear_power_plant_power = exp_settings.getint('Plot', 'nuclear_power_plant_power')
        self.heater_cooler_efficiency = exp_settings.getfloat('Plot', 'heater_cooler_efficiency')
        self.energy_cost = exp_settings.getfloat('Plot', 'energy_cost')


def get_experience_settings(experience_ini_path: str):
    experience_settings = read_ini(experience_ini_path)
    return experience_settings
