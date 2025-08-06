"""
Required Slope & Destination Point Applying Slope For Distance and Bearing
==========================================================================
1. Computes the required climb/descent slope and horizontal distance between two WGS84 waypoints.
2. Computes the destination waypoint (lat, lon, alt) from an initial waypoint, given a slope,
   horizontal distance and bearing.
Uses three different geodesic algorithms: Vincenty, Haversine, and Rhumb-Line.
"""

from pyBADA import conversions as conv
from pyBADA import geodesic as geo

# ——— Inputs for required–slope example ———
waypoint1 = {
    "latitude": 52.2367946579192,
    "longitude": 20.7129809016565,
    "altitude": 3500.0,  # ft
}
waypoint2 = {
    "latitude": 52.1697191213371,
    "longitude": 20.9519554471793,
    "altitude": 412.0,  # ft
}

# ——— Inputs for destination–point example ———
initial_waypoint = waypoint1.copy()
slope_deg = 3.5  # degrees (positive = climb)
distance_nm = 12.0  # nautical miles
bearing_deg = 75.0  # degrees from true north

# ——— Run both examples for each algorithm ———
for algo_name, Algo in [
    ("Vincenty", geo.Vincenty),
    ("Haversine", geo.Haversine),
    ("RhumbLine", geo.RhumbLine),
]:
    print(f"\n=== {algo_name} ===")

    # 1) Required slope & horizontal distance
    slope_req_deg, dist_m = Algo.requiredSlope(waypoint1, waypoint2)
    dist_nm = conv.m2nm(dist_m)
    print("Required slope from waypoint1 → waypoint2:")
    print(f"  Slope               = {slope_req_deg:.8f}°")
    print(f"  Horizontal distance = {dist_nm:.8f} NM")

    # 2) Destination point given a slope, distance & bearing
    dest_wp = Algo.destinationPointApplyingSlopeForDistance(
        initial_waypoint, slope_deg, distance_nm, bearing_deg
    )
    print("Destination waypoint applying slope & distance:")
    print(f"  Latitude  = {dest_wp['latitude']:.6f}°")
    print(f"  Longitude = {dest_wp['longitude']:.6f}°")
    print(f"  Altitude  = {dest_wp['altitude']:.2f} ft")
