"""
Destination Point Applying Slope For Distance and Bearing
=========================================================

Example of computing the destination waypoint (latitude, longitude, altitude)
starting from an initial WGS84 waypoint, flying a given slope and horizontal
distance on a specified bearing, using three different geodesic algorithms.

"""

from pyBADA import conversions as conv
from pyBADA import geodesic as geo

# Common inputs
initial_alt_ft = 3500.0         # feet
slope_deg      = 3.5            # degrees (positive = climb, negative = descent)
distance_nm    = 12.0           # nautical miles
waypoint1 = {
    'latitude': 52.2367946579192,
    'longitude': 20.7129809016565,
    'altitude': 3500
}
bearing_deg = 75.0  # degrees from true north

# Iterate over each geodesic algorithm
for algo_name, Algo in [
    ("Vincenty", geo.Vincenty),
    ("Haversine", geo.Haversine),
    ("RhumbLine", geo.RhumbLine)
]:
    print(f"\n=== {algo_name} ===")

    dest_wp = Algo.destinationPointApplyingSlopeForDistance(
        waypoint1,
        slope_deg,
        distance_nm,
        bearing_deg
    )
    print("Destination waypoint:")
    print(f"  Latitude  = {dest_wp['latitude']:.6f}°")
    print(f"  Longitude = {dest_wp['longitude']:.6f}°")
    print(f"  Altitude  = {dest_wp['altitude']:.2f} ft")