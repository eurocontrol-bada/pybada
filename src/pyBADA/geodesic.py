"""
Geodesic calculation module
"""

from math import (
    tan,
    atan2,
    sin,
    asin,
    cos,
    radians,
    degrees,
    sqrt,
    pi,
    log,
    log2,
)
from pyBADA.aircraft import Airplane as airplane
from pyBADA import constants as const


class Haversine:
    """This class implements the geodesic calculations on sherical earth
    (ignoring ellipsoidal effects).

    .. note::
            https://www.movable-type.co.uk/scripts/latlong.html
    """

    def __init__(self):
        pass

    @staticmethod
    def distance(LAT_init, LON_init, LAT_final, LON_final):
        """Calculate the great-circle distance between two points on the
        Earth's surface using the haversine formula.

        The great-circle distance is the shortest distance between two points
        over the Earth's surface, ignoring elevation changes (i.e., hills or
        mountains).

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: Great-circle distance between the two points in meters.
        :rtype: float
        """

        phi_init = radians(LAT_init)
        phi_final = radians(LAT_final)
        delta_phi = phi_final - phi_init
        lambda_init = radians(LON_init)
        lambda_final = radians(LON_final)
        delta_lambda = lambda_final - lambda_init

        a = pow(sin(delta_phi / 2), 2) + cos(phi_init) * cos(phi_final) * pow(
            sin(delta_lambda / 2), 2
        )
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
        d = const.AVG_EARTH_RADIUS_KM * 1000 * c

        return d

    @staticmethod
    def destinationPoint(LAT_init, LON_init, distance, bearing):
        """Calculate the destination point given an initial point, distance,
        and bearing.

        Given an initial latitude and longitude, this function calculates the
        destination point after traveling the specified distance along the
        given initial bearing (direction).

        Note that the bearing may vary along the path, but this calculation
        assumes a constant bearing.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param distance: Distance to travel from the initial point in meters.
        :param bearing: Initial bearing (direction) in degrees.
        :type LAT_init: float
        :type LON_init: float
        :type distance: float
        :type bearing: float
        :returns: Tuple containing the destination latitude and longitude in
            degrees.
        :rtype: (float, float)
        """

        delta = distance / (const.AVG_EARTH_RADIUS_KM * 1000)
        theta = radians(bearing)

        phi_init = radians(LAT_init)
        lambda_init = radians(LON_init)

        sinPhi_final = sin(phi_init) * cos(delta) + cos(phi_init) * sin(
            delta
        ) * cos(theta)
        phi_final = asin(sinPhi_final)
        y = sin(theta) * sin(delta) * cos(phi_init)
        x = cos(delta) - sin(phi_init) * sinPhi_final
        lambda_final = lambda_init + atan2(y, x)

        lat = degrees(phi_final)
        lon = degrees(lambda_final)

        return (lat, lon)

    @staticmethod
    def bearing(LAT_init, LON_init, LAT_final, LON_final):
        """Calculate the initial bearing between two points along a great-
        circle path.

        The initial bearing (forward azimuth) is the direction one would need
        to travel in a straight line along the great-circle arc from the start
        point to the end point.

        This bearing is measured clockwise from true north.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: Initial bearing in degrees (0° to 360°).
        :rtype: float
        """

        bearing = atan2(
            sin(radians(LON_final) - radians(LON_init))
            * cos(radians(LAT_final)),
            cos(radians(LAT_init)) * sin(radians(LAT_final))
            - sin(radians(LAT_init))
            * cos(radians(LAT_final))
            * cos(radians(LON_final) - radians(LON_init)),
        )
        bearing = (degrees(bearing) + 360) % 360

        return bearing


class Vincenty:
    """This class implements the vincenty calculations of geodesics on the
    ellipsoid-model earth.

    .. note::
            https://www.movable-type.co.uk/scripts/latlong-vincenty.html
    """

    @staticmethod
    def distance_bearing(LAT_init, LON_init, LAT_final, LON_final):
        """Calculate the geodesic distance, initial bearing, and final bearing
        between two points on the Earth's surface.

        This method uses the Vincenty formula to account for the Earth's
        ellipsoidal shape, providing accurate calculations for long distances.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: Tuple containing distance in meters, initial bearing in
            degrees, and final bearing in degrees.
        :rtype: (float, float, float)
        """

        LON2 = radians(LON_final)
        LON1 = radians(LON_init)
        LAT2 = radians(LAT_final)
        LAT1 = radians(LAT_init)

        L = LON2 - LON1
        tanU1 = (1 - const.f) * tan(LAT1)
        cosU1 = 1 / sqrt(1 + tanU1 * tanU1)
        sinU1 = tanU1 * cosU1
        tanU2 = (1 - const.f) * tan(LAT2)
        cosU2 = 1 / sqrt(1 + tanU2 * tanU2)
        sinU2 = tanU2 * cosU2

        antipodal = False
        if abs(L) > pi / 2 or abs(LAT2 - LAT1) > pi / 2:
            antipodal = True

        lambd = L
        lambd_new = 0.0
        iterations = 0
        while iterations == 0 or (
            abs(lambd - lambd_new) > 1e-12 and iterations < 1000
        ):
            iterations += 1
            sinlambda = sin(lambd)
            coslambda = cos(lambd)
            sinSqsigma = pow((cosU2 * sinlambda), 2) + pow(
                (cosU1 * sinU2 - sinU1 * cosU2 * coslambda), 2
            )
            sinsigma = sqrt(sinSqsigma)
            cossigma = sinU1 * sinU2 + cosU1 * cosU2 * coslambda
            sigma = atan2(sinsigma, cossigma)
            sinalpha = cosU1 * cosU2 * sinlambda / sinsigma
            cosSqalpha = 1 - pow(sinalpha, 2)

            if cosSqalpha != 0.0:
                cos2sigmam = cossigma - 2 * sinU1 * sinU2 / cosSqalpha
            else:
                cos2sigmam = 0.0

            C = (
                const.f
                / 16
                * cosSqalpha
                * (4 + const.f * (4 - 3 * cosSqalpha))
            )
            lambd_new = lambd
            lambd = L + (1 - C) * const.f * sinalpha * (
                sigma
                + C
                * sinsigma
                * (
                    cos2sigmam
                    + C * cossigma * (-1 + 2 * cos2sigmam * cos2sigmam)
                )
            )

            if antipodal:
                iterationcheck = abs(lambd) - pi
            else:
                iterationcheck = abs(lambd)

            if iterationcheck > pi:
                return [None, None, None]

        # vincenty formula failed to converge
        if iterations >= 1000:
            return [None, None, None]

        uSq = (
            cosSqalpha * (pow(const.a, 2) - pow(const.b, 2)) / pow(const.b, 2)
        )
        A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
        B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
        deltaSigma = (
            B
            * sinsigma
            * (
                cos2sigmam
                + B
                / 4
                * (
                    cossigma * (-1 + 2 * pow(cos2sigmam, 2))
                    - B
                    / 6
                    * cos2sigmam
                    * (-3 + 4 * pow(sinsigma, 2))
                    * (-3 + 4 * pow(cos2sigmam, 2))
                )
            )
        )

        s = const.b * A * (sigma - deltaSigma)

        # initial bearing
        alpha1 = atan2(
            cosU2 * sinlambda, cosU1 * sinU2 - sinU1 * cosU2 * coslambda
        )
        alpha1 = (degrees(alpha1) + 360) % 360

        # final bearing
        alpha2 = atan2(
            cosU1 * sinlambda, -sinU1 * cosU2 + cosU1 * sinU2 * coslambda
        )
        alpha2 = (degrees(alpha2) + 360) % 360

        return (s, alpha1, alpha2)

    @staticmethod
    def distance(LAT_init, LON_init, LAT_final, LON_final):
        """Calculate the geodesic distance between two latitude/longitude
        points on the Earth's surface.

        This method uses an accurate ellipsoidal model (Vincenty's formula)
        for calculating the distance, which is particularly useful for long
        distances across the globe.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: The geodesic distance in meters.
        :rtype: float
        """

        dist_bearing = Vincenty.distance_bearing(
            LAT_init, LON_init, LAT_final, LON_final
        )
        return dist_bearing[0]

    @staticmethod
    def bearing_initial(LAT_init, LON_init, LAT_final, LON_final):
        """Calculate the initial bearing (forward azimuth) from the initial
        point to the final point.

        This function returns the initial bearing that, if followed in a
        straight line along a great-circle path, will take you from the start
        point to the end point.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: The initial bearing in degrees.
        :rtype: float
        """

        b_initial = Vincenty.distance_bearing(
            LAT_init, LON_init, LAT_final, LON_final
        )
        return b_initial[1]

    @staticmethod
    def bearing_final(LAT_init, LON_init, LAT_final, LON_final):
        """Calculate the final bearing (reverse azimuth) from the final point
        to the initial point.

        This function calculates the final bearing at the destination point,
        which is the direction one would need to take to return to the initial
        point along the great-circle path.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: The final bearing in degrees.
        :rtype: float
        """

        b_final = Vincenty.distance_bearing(
            LAT_init, LON_init, LAT_final, LON_final
        )
        return b_final[2]

    @staticmethod
    def destinationPoint_finalBearing(LAT_init, LON_init, distance, bearing):
        """Calculate the destination point and final bearing given an initial
        point, distance, and bearing.

        This method calculates the latitude and longitude of the destination point after traveling a specified
        distance along the given bearing from the starting point. It also returns the final bearing at the
        destination point.

        Note: The bearing normally varies along the path due to the Earth's curvature.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param distance: Distance traveled from the initial point in meters.
        :param bearing: Initial bearing (direction) in degrees.
        :type LAT_init: float
        :type LON_init: float
        :type distance: float
        :type bearing: float
        :returns: Tuple containing the destination latitude, destination longitude, and final bearing (degrees).
        :rtype: (float, float, float)
        """

        LON1 = radians(LON_init)
        LAT1 = radians(LAT_init)

        sinalpha1 = sin(radians(bearing))
        cosalpha1 = cos(radians(bearing))

        tanU1 = (1 - const.f) * tan(LAT1)
        cosU1 = 1 / sqrt(1 + tanU1 * tanU1)
        sinU1 = tanU1 * cosU1

        sigma1 = atan2(tanU1, cosalpha1)
        sinalpha = cosU1 * sinalpha1
        cosSqalpha = 1 - pow(sinalpha, 2)
        uSq = (
            cosSqalpha * (pow(const.a, 2) - pow(const.b, 2)) / pow(const.b, 2)
        )
        A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
        B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))

        sigma = distance / (const.b * A)

        sigma_new = 0.0
        iterations = 0
        while iterations == 0 or (
            abs(sigma - sigma_new) > 1e-12 and iterations < 1000
        ):
            iterations += 1
            cos2sigmam = cos(2 * sigma1 + sigma)
            sinsigma = sin(sigma)
            cossigma = cos(sigma)
            deltaSigma = (
                B
                * sinsigma
                * (
                    cos2sigmam
                    + B
                    / 4
                    * (
                        cossigma * (-1 + 2 * cos2sigmam * cos2sigmam)
                        - B
                        / 6
                        * cos2sigmam
                        * (-3 + 4 * sinsigma * sinsigma)
                        * (-3 + 4 * cos2sigmam * cos2sigmam)
                    )
                )
            )
            sigma_new = sigma
            sigma = distance / (const.b * A) + deltaSigma

        # vincenty formula failed to converge
        if iterations >= 1000:
            return [None, None, None]
        # print(distance, sigma,sigma_new,abs(sigma-sigma_new))
        # print(sinsigma)
        x = sinU1 * sinsigma - cosU1 * cossigma * cosalpha1
        LAT2 = atan2(
            sinU1 * cossigma + cosU1 * sinsigma * cosalpha1,
            (1 - const.f) * sqrt(sinalpha * sinalpha + x * x),
        )
        lambd = atan2(
            sinsigma * sinalpha1,
            cosU1 * cossigma - sinU1 * sinsigma * cosalpha1,
        )
        C = const.f / 16 * cosSqalpha * (4 + const.f * (4 - 3 * cosSqalpha))
        L = lambd - (1 - C) * const.f * sinalpha * (
            sigma
            + C
            * sinsigma
            * (cos2sigmam + C + cossigma * (-1 + 2 * cos2sigmam * cos2sigmam))
        )
        LON2 = LON1 + L

        alpha2 = atan2(sinalpha, -x)
        finalBearing = (degrees(alpha2) + 360) % 360

        return (degrees(LAT2), degrees(LON2), finalBearing)

    @staticmethod
    def destinationPoint(LAT_init, LON_init, distance, bearing):
        """Calculate the destination point after traveling a specified
        distance on a given bearing.

        This method returns the latitude and longitude of the destination
        point after traveling the given distance on the specified initial
        bearing, following a great-circle path.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param distance: Distance to be traveled from the initial point in
            meters.
        :param bearing: Initial bearing (direction) in degrees.
        :type LAT_init: float
        :type LON_init: float
        :type distance: float
        :type bearing: float
        :returns: Tuple containing the destination latitude and longitude in
            degrees.
        :rtype: (float, float)
        """

        dest = Vincenty.destinationPoint_finalBearing(
            LAT_init, LON_init, distance, bearing
        )

        return (dest[0], dest[1])


class RhumbLine:
    """This class implements the rhumb line (loxodrome) calculations of
    geodesics on the ellipsoid-model earth.

    .. note::
            https://github.com/SpyrosMouselinos/distancly/blob/master/distancly/rhumbline.py
    """

    @staticmethod
    def simple_project(latitiude: float) -> float:
        """Applies a projection to the latitude for use in rhumbline
        calculations.

        The projection is based on the Mercator projection, where latitudes
        are projected to account for the curvature of the Earth. This formula
        ensures that the calculations along the rhumbline are accurate.

        :param latitiude: Latitude in radians.
        :return: The projected latitude in radians.
        """

        return tan(pi / 4 + latitiude / 2)

    @staticmethod
    def distance(LAT_init, LON_init, LAT_final, LON_final) -> float:
        """Calculates the rhumbline distance between two geographical points
        in meters.

        The rhumbline is a path of constant bearing that crosses all meridians
        at the same angle, unlike a great-circle route which is the shortest
        distance between two points on the Earth's surface.

        This method adjusts for longitudes that span more than half of the
        globe.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :return: The rhumbline distance in meters.
        """

        lat_a = radians(LAT_init)
        lon_a = radians(LON_init)

        lat_b = radians(LAT_final)
        lon_b = radians(LON_final)

        delta_phi = lat_b - lat_a
        delta_psi = log(
            RhumbLine.simple_project(lat_b) / RhumbLine.simple_project(lat_a)
        )
        delta_lambda = lon_b - lon_a

        if abs(delta_psi) > 10e-12:
            q = delta_phi / delta_psi
        else:
            q = cos(lat_a)

        if abs(delta_lambda) > pi:
            if delta_lambda > 0:
                delta_lambda = -(2 * pi - delta_lambda)
            else:
                delta_lambda = 2 * pi + delta_lambda

        dist = (
            sqrt(delta_phi * delta_phi + q * q * delta_lambda * delta_lambda)
            * const.AVG_EARTH_RADIUS_KM
        )
        return dist * 1000

    @staticmethod
    def bearing(LAT_init, LON_init, LAT_final, LON_final) -> float:
        """Calculates the rhumbline bearing from the initial point to the
        final point.

        This returns the constant bearing (direction) required to travel along
        a rhumbline between the two points. The bearing is adjusted for
        longitudes that cross the 180-degree meridian.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :return: The rhumbline bearing in degrees.
        """

        lat_a = radians(LAT_init)
        lon_a = radians(LON_init)

        lat_b = radians(LAT_final)
        lon_b = radians(LON_final)

        delta_psi = log(
            RhumbLine.simple_project(lat_b) / RhumbLine.simple_project(lat_a)
        )
        delta_lambda = lon_b - lon_a

        if abs(delta_lambda) > pi:
            if delta_lambda > 0:
                delta_lambda = -(2 * pi - delta_lambda)
            else:
                delta_lambda = 2 * pi + delta_lambda

        return degrees(atan2(delta_lambda, delta_psi)) % 360

    @staticmethod
    def destinationPoint(LAT_init, LON_init, bearing, distance) -> tuple:
        """Calculates the destination point given an initial point, a bearing,
        and a distance traveled.

        This method computes the final latitude and longitude after traveling
        along a rhumbline for a given distance in meters from the initial
        point at a constant bearing.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param bearing: The constant bearing in degrees.
        :param distance: The distance to travel from the initial point in
            meters.
        :return: A tuple containing the destination latitude and longitude in
            degrees.
        """

        lat_a = radians(LAT_init)
        lon_a = radians(LON_init)
        theta = radians(bearing)
        delta = (distance / 1000) / const.AVG_EARTH_RADIUS_KM
        delta_phi = delta * cos(theta)
        lat_b = lat_a + delta_phi
        delta_psi = log(
            RhumbLine.simple_project(lat_b) / RhumbLine.simple_project(lat_a)
        )

        if abs(delta_psi) > 10e-12:
            q = delta_phi / delta_psi
        else:
            q = cos(lat_a)

        delta_lambda = delta * sin(theta) / q
        lon_b = lon_a + delta_lambda

        # Normalise latitude
        if abs(lat_b) > pi / 2:
            if lat_b > 0:
                lat_b = pi - lat_b
            else:
                lat_b = -pi - lat_b

        lat_b = degrees(lat_b)
        lon_b = degrees(lon_b)
        # Normalize longitude
        lon_b = (540 + lon_b) % 360 - 180
        return (lat_b, lon_b)

    @staticmethod
    def loxodromic_mid_point(
        LAT_init, LON_init, LAT_final, LON_final
    ) -> tuple:
        """Calculates the midpoint along a rhumbline between two geographical
        points.

        The midpoint is calculated using the rhumbline path between the
        initial and final points. This takes into account the Earth's
        curvature by projecting the latitudes.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :return: A tuple representing the midpoint's latitude and longitude in
            degrees.
        """

        lat_a = radians(LAT_init)
        lon_a = radians(LON_init)

        lat_b = radians(LAT_final)
        lon_b = radians(LON_final)

        # Anti - Meridian Crossing
        if abs(lon_b - lon_a) > pi:
            lon_a += 2 * pi

        lat_mid = (lat_a + lat_b) / 2
        f1 = RhumbLine.simple_project(lat_a)
        f2 = RhumbLine.simple_project(lat_b)
        f3 = RhumbLine.simple_project(lat_mid)
        if abs(f2 - f1) < 1e-6:
            lon_mid = lon_a + lon_b / 2
        else:
            lon_mid = (
                (lon_b - lon_a) * log(f3) + lon_a * log(f2) - lon_b * log(f1)
            ) / log(f2 / f1)

        lat_mid = degrees(lat_mid)
        lon_mid = degrees(lon_mid)
        # Normalize longitude
        lon_mid = (540 + lon_mid) % 360 - 180
        return lat_mid, lon_mid

    @staticmethod
    def loxodromic_power_interpolation(
        LAT_init, LON_init, LAT_final, LON_final, n_points: int
    ) -> list:
        """Generates a specified number of points between two geographical
        locations along a rhumbline path.

        This method recursively calculates intermediate points between two
        points on the Earth's surface, following a constant bearing rhumbline
        path. The number of points should be a power of 2 minus 1.

        :param LAT_init: Initial latitude in degrees.
        :param LON_init: Initial longitude in degrees.
        :param LAT_final: Final latitude in degrees.
        :param LON_final: Final longitude in degrees.
        :param n_points: Number of intermediate points to generate. Must be a
            power of 2 minus 1.
        :return: A list of tuples, where each tuple represents an interpolated
            point's latitude and longitude in degrees.
        """

        n_points = int(n_points)
        if not log2(n_points + 1).is_integer():
            print(
                "N_Points must be an power of 2 minus 1 Number! e.g. 1,3,7,15,..."
            )
            return []

        lmp = RhumbLine.loxodromic_mid_point

        # Recursive Solution #
        def solution(a, b, idx):
            if idx == 1:
                return lmp(a[0], a[1], b[0], b[1])
            else:
                return (
                    solution(a, solution(a, b, 1), (idx - 1) / 2),
                    solution(a, b, 1),
                    solution(solution(a, b, 1), b, (idx - 1) / 2),
                )

        points = solution(
            (LAT_init, LON_init), (LAT_final, LON_final), n_points
        )

        # Decouple points
        decoupled_points = []
        if len(points) == 2:
            decoupled_points.append(points)
        else:
            for midpoint in points:
                decoupled_points.append(midpoint)
        return decoupled_points


class Turn:
    """This class implements the calculations of geodesics turns."""

    @staticmethod
    def destinationPoint_finalBearing(
        LAT_init,
        LON_init,
        bearingInit,
        TAS,
        rateOfTurn,
        timeOfTurn,
        directionOfTurn,
        centerPoint=None,
    ):
        """Calculates the destination point and final bearing after traveling
        for a given time with a specified turn.

        This function computes the aircraft's final position and bearing
        after making a turn at a specified rate of turn, direction, and
        true airspeed (TAS). If TAS is zero, the aircraft rotates in
        place. The calculation accounts for turning radius and bank
        angle.

        :param LAT_init: Initial latitude [deg].
        :param LON_init: Initial longitude [deg].
        :param timeOfTurn: Time spent in turn [s].
        :param bearingInit: Initial bearing [deg].
        :param TAS: True Airspeed (TAS) [m/s].
        :param rateOfTurn: Rate of turn [deg/s].
        :param directionOfTurn: Direction of turn ('LEFT' or 'RIGHT').
        :param centerPoint: Optional latitude and longitude of the
            rotation center (defaults to None) [deg, deg].
        :type LAT_init: float.
        :type LON_init: float.
        :type timeOfTurn: float.
        :type bearingInit: float.
        :type TAS: float.
        :type rateOfTurn: float.
        :type directionOfTurn: str.
        :type centerPoint: tuple(float, float).
        :returns: Destination point's latitude, longitude, and final
            bearing [deg, deg, deg].
        :rtype: tuple(float, float, float).
        """

        if TAS == 0:
            arcLength = (
                rateOfTurn * timeOfTurn
            )  # amount of degrees to do the rotation

            if directionOfTurn == "RIGHT":
                bearing_final = (bearingInit + arcLength) % 360
            elif directionOfTurn == "LEFT":
                bearing_final = (bearingInit - arcLength) % 360

            return (LAT_init, LON_init, bearing_final)

        else:
            bankAngle = airplane.bankAngle(
                rateOfTurn=rateOfTurn, v=TAS
            )  # [degrees]

            arcLength = (
                rateOfTurn * timeOfTurn
            )  # amount of degrees to do the rotation
            turnRadius = airplane.turnRadius_bankAngle(
                v=TAS, ba=bankAngle
            )  # [m]

            # find center of rotation, which is at (bearingInit + 90 degrees) and distance of turnRadius
            if directionOfTurn == "RIGHT":
                centerAngle = bearingInit + 90
            elif directionOfTurn == "LEFT":
                centerAngle = bearingInit - 90

            if centerPoint is None:
                centerPoint = RhumbLine.destinationPoint(
                    LAT_init=LAT_init,
                    LON_init=LON_init,
                    distance=turnRadius,
                    bearing=centerAngle,
                )

            # calcualte new angle after the rotation from the center point to new destination point
            if directionOfTurn == "RIGHT":
                newAngle = (centerAngle + arcLength + 180) % 360
            elif directionOfTurn == "LEFT":
                newAngle = (centerAngle - arcLength + 180) % 360

            # calcualte the new destination point after the rotation from the center point, using the same distance
            finalPoint = RhumbLine.destinationPoint(
                LAT_init=centerPoint[0],
                LON_init=centerPoint[1],
                distance=turnRadius,
                bearing=newAngle,
            )

            if directionOfTurn == "RIGHT":
                bearing_final = (bearingInit + arcLength) % 360
            elif directionOfTurn == "LEFT":
                bearing_final = (bearingInit - arcLength) % 360

            dist = RhumbLine.distance(
                LAT_init=LAT_init,
                LON_init=LON_init,
                LAT_final=centerPoint[0],
                LON_final=centerPoint[1],
            )

            return (finalPoint[0], finalPoint[1], bearing_final)

    @staticmethod
    def distance(rateOfTurn, TAS, timeOfTurn):
        """Calculates the distance traveled during a turn based on the rate of
        turn, true airspeed, and time.

        This function computes the total distance traveled during a constant
        turn, based on the aircraft's rate of turn, true airspeed, and the
        duration of the turn.

        :param rateOfTurn: Rate of turn [deg/s].
        :param TAS: True Airspeed (TAS) [m/s].
        :param timeOfTurn: Duration of the turn [s].
        :type rateOfTurn: float.
        :type TAS: float.
        :type timeOfTurn: float.
        :returns: Distance traveled during the turn [m].
        :rtype: float.
        """

        if TAS == 0:
            return 0
        else:
            bankAngle = airplane.bankAngle(rateOfTurn=rateOfTurn, v=TAS)
            arcLengthDegrees = (
                rateOfTurn * timeOfTurn
            )  # amount of degrees to do the rotation
            turnRadius = airplane.turnRadius_bankAngle(
                v=TAS, ba=bankAngle
            )  # [m]
            distance = radians(arcLengthDegrees) * turnRadius  # arcLength [m]

        return distance
