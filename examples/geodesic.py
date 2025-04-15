"""
Geodesic Calculation
====================

Example of using the geodesic module and comparision
between Haversine and Vincenty implementation

"""

import folium

from pyBADA import geodesic as geo


# input parameters
LAT_init = 48.596289
LON_init = 2.351640
LAT_final = 48.613721
LON_final = 2.376616

distance = 54972.271  # [m]
bearing = 310

# ===========================================================
# calculation of distance between 2 latitude/longitude points
# ===========================================================

# Haversine formula
distance_haversine = geo.Haversine.distance(
    LAT_init=LAT_init,
    LON_init=LON_init,
    LAT_final=LAT_final,
    LON_final=LON_final,
)
print("Haversine distance:", distance_haversine)

# Vincenty formula
distance_vincenty = geo.Vincenty.distance(
    LAT_init=LAT_init,
    LON_init=LON_init,
    LAT_final=LAT_final,
    LON_final=LON_final,
)
print("Vincenty distance: ", distance_vincenty)
print("")

# ===========================================================
# calculation of bearing between 2 latitude/longitude points
# ===========================================================

# Haversine formula
bearing_haversine = geo.Haversine.bearing(
    LAT_init=LAT_init,
    LON_init=LON_init,
    LAT_final=LAT_final,
    LON_final=LON_final,
)
print("Haversine initial bearing:", bearing_haversine)

# Vincenty formula for initial bearing
bearing_vincenty = geo.Vincenty.bearing_initial(
    LAT_init=LAT_init,
    LON_init=LON_init,
    LAT_final=LAT_final,
    LON_final=LON_final,
)
print("Vincenty initial bearing: ", bearing_vincenty)

# Vincenty formula for final bearing
bearing_vincenty = geo.Vincenty.bearing_final(
    LAT_init=LAT_init,
    LON_init=LON_init,
    LAT_final=LAT_final,
    LON_final=LON_final,
)
# print("Vincenty final bearing:   ", bearing_vincenty)
print("")

# =====================================================================================
# calculation of destination point based in initial point, distance and initial bearing
# =====================================================================================

# Haversine formula
destPoint_haversine = geo.Haversine.destinationPoint(
    LAT_init=LAT_init, LON_init=LON_init, distance=distance, bearing=bearing
)
print("Haversine destination point:", destPoint_haversine)

# Vincenty formula
destPoint_vincenty = geo.Vincenty.destinationPoint(
    LAT_init=LAT_init, LON_init=LON_init, distance=distance, bearing=bearing
)
print("Vincenty destination point: ", destPoint_vincenty)
print("")

# =====================================================================================
# calculation of destination point based in initial point, distance and initial bearing
# =====================================================================================

# Haversine formula
destPoint_haversine = geo.Haversine.destinationPoint(
    LAT_init=LAT_init, LON_init=LON_init, distance=distance, bearing=bearing
)
print("Haversine destination point:", destPoint_haversine)

# Vincenty formula
destPoint_vincenty = geo.Vincenty.destinationPoint(
    LAT_init=LAT_init, LON_init=LON_init, distance=distance, bearing=bearing
)
print("Vincenty destination point: ", destPoint_vincenty)

# =====================================================================================
# drawing destination point on the 2D map
# =====================================================================================

initPoint = (LAT_init, LON_init)

# create a base map
myMap = folium.Map(location=initPoint, zoom_start=10)

# construct the lines
haversineLine = []
haversineLine.append(initPoint)
haversineLine.append(destPoint_haversine)

vincentyLine = []
vincentyLine.append(initPoint)
vincentyLine.append(destPoint_vincenty)

# draw the lines on the map
folium.PolyLine(haversineLine, color="red", weight=2).add_to(myMap)
folium.PolyLine(vincentyLine, color="blue", weight=2).add_to(myMap)

# add markers on the map
folium.Marker(initPoint).add_to(myMap)
folium.Marker(
    destPoint_haversine,
    popup="<i>Haversine</i>",
    tooltip="Haversine",
    icon=folium.Icon(color="red"),
).add_to(myMap)
folium.Marker(
    destPoint_vincenty,
    popup="<i>Vincenty</i>",
    tooltip="Vincenty",
    icon=folium.Icon(color="blue"),
).add_to(myMap)

# saving the map in the HTML file. This map then can be viewed in the browser and you can interact with it
# myMap.save("map_geodesic.html")

# display the map
myMap
