import math
import numpy as np

def calculate_value1(ydb):
    return max(0.0, min(1.0, 3 * (2 * ydb - 1)**2 - 2 * (2 * ydb - 1)**3))

def calculate_value(ydb):    
    if ydb < 0.5:
        return 0.0
    elif ydb > 1.0:
        return 1.0
    else:
        return 3 * (2 * ydb - 1)**2 - 2 * (2 * ydb - 1)**3

# Example usage
# yd = np.linspace(0,2,10)

print(calculate_value1(0),calculate_value(0))
print(calculate_value1(0.1),calculate_value(0.1))
print(calculate_value1(0.6),calculate_value(0.6))
print(calculate_value1(1.5),calculate_value(1.5))
print(calculate_value1(10),calculate_value(10))
