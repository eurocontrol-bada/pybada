# -*- coding: utf-8 -*-
"""
pyBADA
Geodesic calculation module
Developped @EUROCONTROL (EIH)
2024
"""

__author__ = "Henrich Glaser-Opitz"
__copyright__ = "Copyright 2024, EUROCONTROL (EIH)"
__license__ = "BADA Eurocontrol"
__version__ = "1.0.0"
__maintainer__ = "Henrich Glaser-Opitz"
__email__ = "henrich.glaser-opitz@eurocontrol.int"
__status__ = "Development"
__docformat__ = "reStructuredText"


from math import tan, atan2, sin, asin, cos, radians, degrees, sqrt, pi, log, log2, acos
from pyBADA.aircraft import Airplane as airplane
from pyBADA import conversions as conv
from pyBADA import constants as const


class Haversine:
    """This class implements the geodesic calculations on sherical earth (ignoring ellipsoidal effects).

    .. note::
            https://www.movable-type.co.uk/scripts/latlong.html
    """

    def __init__(self):
        pass

    @staticmethod
    def distance(LAT_init, LON_init, LAT_final, LON_final):
        """This function returns the great-circle distance between two point using haversine formula.
        That is the shortest distance over the earth's surface (ignoring any hills on earth surface)

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param LAT_final: final Latitude [deg]
        :param LON_final: final Longitude [deg]
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: distance [NM]
        :rtype: float.
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

        return conv.m2nm(d)

    @staticmethod
    def destinationPoint(LAT_init, LON_init, distance, bearing):
        """This function returns the destination point having travelled the given distance on the initial bearing

        .. note::
                bearing normally varies around path followed

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param distance: distance travelled from initial point at defined bearing [NM]
        :param bearing:  [deg]
        :type LAT_init: float
        :type LON_init: float
        :type distance: float
        :type bearing: float
        :returns: detination point ([deg],[deg])
        :rtype: (float, float).
        """

        delta = conv.nm2m(distance) / (const.AVG_EARTH_RADIUS_KM * 1000)
        theta = radians(bearing)

        phi_init = radians(LAT_init)
        lambda_init = radians(LON_init)

        sinPhi_final = sin(phi_init) * cos(delta) + cos(phi_init) * sin(delta) * cos(
            theta
        )
        phi_final = asin(sinPhi_final)
        y = sin(theta) * sin(delta) * cos(phi_init)
        x = cos(delta) - sin(phi_init) * sinPhi_final
        lambda_final = lambda_init + atan2(y, x)

        lat = degrees(phi_final)
        lon = degrees(lambda_final)

        return (lat, lon)

    @staticmethod
    def bearing(LAT_init, LON_init, LAT_final, LON_final):
        """This function returns the initial bearing (sometimes referred to as forward azimuth), which if followed
        in a straight line along the great-circle arc will take you from start point to the end point

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param LAT_final: final Latitude [deg]
        :param LON_final: final Longitude [deg]
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: initial bearing [deg]
        :rtype: float.
        """

        bearing = atan2(
            sin(radians(LON_final) - radians(LON_init)) * cos(radians(LAT_final)),
            cos(radians(LAT_init)) * sin(radians(LAT_final))
            - sin(radians(LAT_init))
            * cos(radians(LAT_final))
            * cos(radians(LON_final) - radians(LON_init)),
        )
        bearing = (degrees(bearing) + 360) % 360

        return bearing


class Vincenty(object):
    """This class implements the vincenty calculations of geodesics on the ellipsoid-model earth

    .. note::
            https://www.movable-type.co.uk/scripts/latlong-vincenty.html

    """

    @staticmethod
    def distance_bearing(LAT_init, LON_init, LAT_final, LON_final):
        """This function returns the geodesic distance, initial and final bearing between
        a pair of latitude/longitude points on the earth's surface using an accurate ellipsoidal
        model of earth

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param LAT_final: final Latitude [deg]
        :param LON_final: final Longitude [deg]
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: [distance, initial bearing, final bearing] [m, deg, deg]
        :rtype: (float, float, float).
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
        while iterations == 0 or (abs(lambd - lambd_new) > 1e-12 and iterations < 1000):
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

            C = const.f / 16 * cosSqalpha * (4 + const.f * (4 - 3 * cosSqalpha))
            lambd_new = lambd
            lambd = L + (1 - C) * const.f * sinalpha * (
                sigma
                + C
                * sinsigma
                * (cos2sigmam + C * cossigma * (-1 + 2 * cos2sigmam * cos2sigmam))
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

        uSq = cosSqalpha * (pow(const.a, 2) - pow(const.b, 2)) / pow(const.b, 2)
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
        alpha1 = atan2(cosU2 * sinlambda, cosU1 * sinU2 - sinU1 * cosU2 * coslambda)
        alpha1 = (degrees(alpha1) + 360) % 360

        # final bearing
        alpha2 = atan2(cosU1 * sinlambda, -sinU1 * cosU2 + cosU1 * sinU2 * coslambda)
        alpha2 = (degrees(alpha2) + 360) % 360

        return (s, alpha1, alpha2)

    @staticmethod
    def distance(LAT_init, LON_init, LAT_final, LON_final):
        """This function returns the geodesic distance between a pair of latitude/longitude points
        on the earth's surface using an accurate ellipsoidal model of earth

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param LAT_final: final Latitude [deg]
        :param LON_final: final Longitude [deg]
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: distance [m]
        :rtype: float.
        """

        dist = Vincenty.distance_bearing(LAT_init, LON_init, LAT_final, LON_final)
        return dist[0]

    @staticmethod
    def bearing_initial(LAT_init, LON_init, LAT_final, LON_final):
        """This function returns the initial bearing between a pair of latitude/longitude points
        on the earth's surface using an accurate ellipsoidal model of earth

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param LAT_final: final Latitude [deg]
        :param LON_final: final Longitude [deg]
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: initial bearing [deg]
        :rtype: float.
        """

        b_initial = Vincenty.distance_bearing(LAT_init, LON_init, LAT_final, LON_final)
        return b_initial[1]

    @staticmethod
    def bearing_final(LAT_init, LON_init, LAT_final, LON_final):
        """This function returns the final bearing between a pair of latitude/longitude points
        on the earth's surface using an accurate ellipsoidal model of earth

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param LAT_final: final Latitude [deg]
        :param LON_final: final Longitude [deg]
        :type LAT_init: float
        :type LON_init: float
        :type LAT_final: float
        :type LON_final: float
        :returns: final bearing [deg]
        :rtype: float.
        """

        b_final = Vincenty.distance_bearing(LAT_init, LON_init, LAT_final, LON_final)
        return b_final[2]

    @staticmethod
    def destinationPoint_finalBearing(LAT_init, LON_init, distance, bearing):
        """This function returns the destination point and final bearing having travelled the given distance on the initial bearing
        from the initial point

        .. note::
                bearing normally varies around path followed

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param distance: distance travelled from initial point at defined bearing [m]
        :param bearing:  initial bearing [deg]
        :type LAT_init: float
        :type LON_init: float
        :type distance: float
        :type bearing: float
        :returns: [detination point and final bearing] ([deg],[deg],[deg])
        :rtype: (float, float, float).
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
        uSq = cosSqalpha * (pow(const.a, 2) - pow(const.b, 2)) / pow(const.b, 2)
        A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
        B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))

        sigma = distance / (const.b * A)

        sigma_new = 0.0
        iterations = 0
        while iterations == 0 or (abs(sigma - sigma_new) > 1e-12 and iterations < 1000):
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
            sinsigma * sinalpha1, cosU1 * cossigma - sinU1 * sinsigma * cosalpha1
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
        """This function returns the destination point having travelled the given distance on the initial bearing
        from the initial point

        .. note::
                bearing normally varies around path followed

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param distance: distance travelled from initial point at defined bearing [m]
        :param bearing:  initial bearing [deg]
        :type LAT_init: float
        :type LON_init: float
        :type distance: float
        :type bearing: float
        :returns: detination point ([deg],[deg])
        :rtype: (float, float).
        """

        dest = Vincenty.destinationPoint_finalBearing(
            LAT_init, LON_init, distance, bearing
        )

        return (dest[0], dest[1])


class RhumbLine(object):
    """This class implements the rhumb line (loxodrome) calculations of geodesics on the ellipsoid-model earth

    .. note::
            https://github.com/SpyrosMouselinos/distancly/blob/master/distancly/rhumbline.py

    """

    @staticmethod
    def simple_project(latitiude: float) -> float:
        """
        Projects a point to its corrected latitude for the rhumbline calculations.
        :param latitiude: A float in radians.
        :return: The projected value in radians.
        """

        return tan(pi / 4 + latitiude / 2)

    @staticmethod
    def distance(LAT_init, LON_init, LAT_final, LON_final) -> float:
        """
        Returns Rhumbline distance in [m] between two points.
        :param point_a: Start point. Tuple of degrees.
        :param point_b: End point. Tuple of degrees.
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
        """
        Returns bearing between two points in degrees
        :param point_a: Start point. Tuple of degrees.
        :param point_b: End point. Tuple of degrees.
        :return: The bearing in degrees.
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
        """
        Returns point B from point A, travelling at constant bearing θ in deg, and distance d in [m].
        :param point_a: Start point. Tuple of degrees.
        :param bearing: The Bearing in degrees.
        :param distance: The Distance in m.
        :return: Point B coordinates.
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
    def loxodromic_mid_point(LAT_init, LON_init, LAT_final, LON_final) -> tuple:
        """
        Finds the rhumbline mid point between 2 points
        :param point_a: Start point. Tuple of degrees.
        :param point_b: End point. Tuple of degrees.
        :return: The rhumbline midpoint tuple.
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
        """
        Returns n_points points between point_a and point_b
        according to the rhumbline loxodromic interpolation,
        using recursive programming.
        :param point_a: Start point. Tuple of degrees.
        :param point_b: End point. Tuple of degrees.
        :param n_points: Number of midpoints. Needs to be a Power of 2 minus 1.
        :returns The list of interpolation points from start to end.
        """
        n_points = int(n_points)
        if not log2(n_points + 1).is_integer():
            print("N_Points must be an power of 2 minus 1 Number! e.g. 1,3,7,15,...")
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

        points = solution((LAT_init, LON_init), (LAT_final, LON_final), n_points)

        # Decouple points
        decoupled_points = []
        if len(points) == 2:
            decoupled_points.append(points)
        else:
            for midpoint in points:
                decoupled_points.append(midpoint)
        return decoupled_points


class Turn(object):
    """This class implements the calculations of geodesics turns

    .. note::

    """

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
        """This function returns the destination point and final bearing having travelled the given time Of Turn from the initial bearing
        from the initial point while turning with Rate of Turn

        :param LAT_init: initial Latitude [deg]
        :param LON_init: initial Longitude [deg]
        :param timeOfTurn: time travelled from initial point from defined initial bearing [s]
        :param bearingInit: initial bearing [deg]
        :param TAS: aircraft True Airspeed TAS [m/s]
        :param rateOfTurn: rate of turn [deg/s]
        :param directionOfTurn: direction of turn {LEFT/RIGHT}
        :param centerPoint: (LAT, LON) of point of rotation [deg, deg]
        :type LAT_init: float
        :type LON_init: float
        :type timeOfTurn: float
        :type bearing: float
        :type TAS: float
        :type rateOfTurn: float
        :type directionOfTurn: string
        :type centerPoint: tuple(float, float)
        :returns: [detination point and final bearing] ([deg],[deg],[deg])
        :rtype: (float, float, float).
        """

        if TAS == 0:
            arcLength = rateOfTurn * timeOfTurn  # amount of degrees to do the rotation

            if directionOfTurn == "RIGHT":
                bearing_final = (bearingInit + arcLength) % 360
            elif directionOfTurn == "LEFT":
                bearing_final = (bearingInit - arcLength) % 360

            return (LAT_init, LON_init, bearing_final)

        else:
            bankAngle = airplane.bankAngle(rateOfTurn=rateOfTurn, v=TAS)  # [degrees]

            arcLength = rateOfTurn * timeOfTurn  # amount of degrees to do the rotation
            turnRadius = airplane.turnRadius_bankAngle(v=TAS, ba=bankAngle)  # [m]

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
        """This function returns the distance travelled the given time of turn
        while turning with Rate of Turn and TAS

        :param timeOfTurn: time travelled from initial point from defined initial bearing [s]
        :param TAS: aircraft True Airspeed TAS [m/s]
        :param rateOfTurn: rate of turn [deg/s]
        :type timeOfTurn: float
        :type TAS: float
        :type rateOfTurn: float
        :returns: distance [m]
        :rtype: float.
        """

        if TAS == 0:
            return 0
        else:
            bankAngle = airplane.bankAngle(rateOfTurn=rateOfTurn, v=TAS)
            arcLengthDegrees = (
                rateOfTurn * timeOfTurn
            )  # amount of degrees to do the rotation
            turnRadius = airplane.turnRadius_bankAngle(v=TAS, ba=bankAngle)  # [m]
            distance = radians(arcLengthDegrees) * turnRadius  # arcLength [m]

        return distance
