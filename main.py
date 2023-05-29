from Classes.Experiment import Experiment
import Plot.standard_plot as std_plot
import Solvers.Composite_Wall_Solver as Cws
from pathlib import Path

current_path = Path.cwd()

E = Experiment(str((current_path / 'Config' / 'experiment_config.ini').resolve()),
               str((current_path / 'Config' / 'wall_config.csv').resolve()),
               str((current_path / 'Config' / 'materials_config.csv').resolve()),
               str((current_path / 'Weather_Data' / 'Paris_2022-07-01_to_2022-07-31.csv').resolve()))

simulation = Cws.Simulation(E)
simulation.simulate()
std_plot.plot_simulation(simulation)
