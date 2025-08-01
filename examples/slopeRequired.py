"""
Required Slope Across Geodesic Algorithms
=========================================

Example of computing the required climb/descent slope and horizontal
distance between two WGS84 waypoints, using three different geodesic
algorithms.
"""

from pyBADA import conversions as conv
from pyBADA import geodesic as geo

# Define two waypoints with WGS84 lat/lon and pressure altitude (ft)
waypoint1 = {
    'latitude': 52.2367946579192,
    'longitude': 20.7129809016565,
    'altitude': 3500.0
}
waypoint2 = {
    'latitude': 52.1697191213371,
    'longitude': 20.9519554471793,
    'altitude': 412.0
}

# Iterate over each geodesic algorithm
for algo_name, Algo in [
    ("Vincenty",   geo.Vincenty),
    ("Haversine",  geo.Haversine),
    ("RhumbLine",  geo.RhumbLine),
]:
    print(f"\n=== {algo_name} ===")
    # Compute required slope (degrees) and horizontal distance (meters)
    slope_deg, dist_m = Algo.requiredSlope(waypoint1, waypoint2)
    # Convert distance to nautical miles for reporting
    dist_nm = conv.m2nm(dist_m)
    print(f"Required slope = {slope_deg:.8f}°")
    print(f"Horizontal distance = {dist_nm:.8f} NM")