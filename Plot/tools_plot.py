import math
import numpy as np


def time_formatting(seconds: float) -> str:
    time_count = seconds
    partial_times = [1, 60, 60, 24, 365]
    partial_names = ["s", "m", "h", "d", "yrs"]
    partial_results = [0, 0, 0, 0, 0]
    phrase = f''
    for i in range(len(partial_times) - 1, -1, -1):
        product = math.prod(partial_times[:i + 1])
        partial_results[i] = int(time_count // product)
        time_count -= partial_results[i] * product
        if partial_results[i] != 0:
            phrase = phrase + f' {partial_results[i]} {partial_names[i]}'
    return phrase


def scale_formatting_unique(max_value_meters: float):
    final_product = 0.001
    product_names = ["km", "hm", "dam", "m", "dm", "cm", "mm"]
    for i in range(len(product_names)):
        control_value = max_value_meters * final_product
        if control_value < 1 or control_value >= 10:
            final_product *= 10
        else:
            return final_product, product_names[i]


def sample_of_array(array, n):
    result = []
    for i in range(n):
        result.append(array[array.shape[0] // n * i])
    return np.array(result)
