import math

period = 3600 * 24  # A day in seconds
w = 2 * math.pi / period


# This function is called on each time step to update outside temperature
def update_outside_temperature(t: float) -> float:
    return 273.15 + 27 + 5 * math.sin(w * t)
