import Classes.Material
from Tools.csv_tool import read_csv
from Classes.Material import create_pcm_list_from_csv as c_pcm


class Wall:
    def __init__(self, wall_csv_path: str, materials_csv_path: str):
        self.wall_materials, self.thicknesses = get_materials_in_wall(wall_csv_path, materials_csv_path)
        self.length = sum(self.thicknesses)
        self.cells_number = 0
        self.cells_number_by_layer = []
        self.material_repartition = []

    def calculate_cells_number_by_material(self, total_cells: int):
        self.cells_number_by_layer = []
        for thickness in self.thicknesses:
            self.cells_number_by_layer.append(int((thickness / self.length) * total_cells))
        self.cells_number = sum(self.cells_number_by_layer)
        self.fill_properties_repartition()
        return self.cells_number

    def fill_properties_repartition(self):
        for index, cells_in_material in enumerate(self.cells_number_by_layer):
            self.material_repartition += [self.wall_materials[index]] * cells_in_material


def get_materials_in_wall(wall_csv_path: str, materials_csv_path: str) -> (
        list[Classes.Material.PhaseChangeMaterial], list[float]):
    all_materials_available = c_pcm(materials_csv_path)
    materials_in_wall, thicknesses = [], []
    wall_layers = read_csv(wall_csv_path)
    for layer in wall_layers:
        # Get corresponding material from available materials
        for material in all_materials_available:
            if material.name == layer[0]:
                materials_in_wall.append(material)
                thicknesses.append(float(layer[1]))
                continue
    return materials_in_wall, thicknesses
