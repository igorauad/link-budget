"""Satellite look angle computations

References:

- [1] https://www.ngs.noaa.gov/CORS/Articles/SolerEisemannJSE.pdf.
- [2] https://en.wikipedia.org/wiki/Earth_radius.
- [3] Maral, Gerard., Sun, Zhili., Bousquet, Michel. Satellite Communications
  Systems: Systems, Techniques and Technology. 3rd ed.

"""
import logging
import os
from datetime import datetime
from math import sqrt, sin, cos, tan, atan, atan2, degrees, radians
from urllib.parse import quote

import numpy as np
from skyfield.api import load, wgs84

from . import util

CELESTRAK_GROUPS = [
    'active', 'stations', 'geo', 'intelsat', 'ses', 'iridium', 'iridium-NEXT',
    'starlink', 'oneweb', 'orbcomm', 'globalstar', 'swarm', 'amateur',
    'x-comm', 'other-comm', 'satnogs', 'gorizont', 'raduga', 'molniya',
    'weather', 'noaa', 'goes', 'resource', 'sarsat', 'dmc', 'tdrss', 'argos',
    'planet', 'spire', 'gnss', 'gps-ops', 'glo-ops', 'galileo', 'beidou',
    'sbas', 'nnss', 'musson', 'military', 'radar', 'cubesat', 'other'
]
CELESTRAK_BASE_URL = "https://celestrak.org/NORAD/elements/gp.php"
skyfield_ts = load.timescale()


def get_default_tle_dataset_dir():
    """Get the default directory for saving TLE datasets"""
    return os.path.join(util.get_default_lb_dir(), "tle")


def get_sat_pos_by_tle(name,
                       group=None,
                       obs_time=None,
                       save_dir=None,
                       no_save=False):
    """Get satellite position from Two-Line Element (TLE) predictions

    First, searches for the satellite on CelesTrak's TLE database. Then,
    computes the expected satellite position on the given observation time
    using the implementation from the Skyfield package.

    The TLE dataset is queried using the format described in
    https://celestrak.org/NORAD/documentation/gp-data-formats.php. More
    specifically, it is either requested with:

    https://celestrak.org/NORAD/elements/gp.php?NAME=VALUE&FORMAT=tle

    or

    https://celestrak.org/NORAD/elements/gp.php?GROUP=VALUE&FORMAT=tle

    The former is used when the satellite is specified by name only (group set
    to None). The second format is used when the specific group on CelesTrak's
    database is also informed. See the groups at
    https://celestrak.org/NORAD/elements/. For example, the 'active' group
    holds a long list of active satellites, whereas the 'geo' group focuses on
    active geosynchronous satellites, and so on.

    This function downloads the TLE dataset before processing, and the
    downloads are saved as txt files at `~/.link-budget/tle/` by default. When
    the dataset is already available locally (downloaded previously), this
    function does not need to download it again. This behavior can be
    customized using the `save_dir` and `no_save` arguments.

    Args:
        name (str): Satellite name on the TLE database.
        group (str, optional): Group to which the satellite belongs on the
            CelesTrak database. Defaults to None, in which case the satellite
            is searched by name over the entire database.
        obs_time (datetime, optional): Observation time. Defaults to
            None, in which case the current time is used.
        save_dir (str, optional): Directory on which the downloaded TLE dataset
            should be saved. Defaults to None, in which case the
            `~/.link-budget/tle/` directory is used.
        no_save (bool, optional): Whether to save the downloaded TLE dataset
            permanently on the save directory (save_dir) to speed up future
            queries.

    Returns:
        tuple: Satellite longitude (degrees), latitude (degrees), and altitude
        (meters).
    """
    name = name.upper()  # ensure upper case

    if save_dir is None:
        save_dir = get_default_tle_dataset_dir()

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    if group is not None:
        url = CELESTRAK_BASE_URL + "?GROUP={}&FORMAT=tle".format(group)
        filename = 'tle-group-{}.txt'.format(group)
    else:
        url = CELESTRAK_BASE_URL + "?NAME={}&FORMAT=tle".format(quote(name))
        filename = 'tle-name-{}.txt'.format(name.replace(" ", "-"))

    filename = os.path.join(save_dir, filename)
    previously_downloaded = os.path.exists(filename)
    satellites = {
        sat.name: sat
        for sat in load.tle_file(url, filename=filename)
    }

    if name not in satellites:
        logging.error(
            "Satellite \"{}\" not found on the TLE database.".format(name))
        closest_options = []
        for key in satellites.keys():
            if name in key:
                closest_options.append(key)
        if closest_options:
            logging.info("Perhaps you mean one of the following options?")
            for option in closest_options:
                logging.info("- " + option)
        if no_save and not previously_downloaded:  # just downloaded
            os.remove(filename)
        raise ValueError("Invalid satellite name")

    satellite = satellites[name]

    obs_time = datetime.utcnow() if obs_time is None else obs_time
    t = skyfield_ts.utc(obs_time.year, obs_time.month, obs_time.day,
                        obs_time.hour, obs_time.minute, obs_time.second)
    geocentric = satellite.at(t)
    lat, lon = wgs84.latlon_of(geocentric)
    height = wgs84.height_of(geocentric)

    util.log_result(
        "Satellite Pos (Lat, Lon, Alt)",
        "{:.2f}°, {:.2f}°, {:.2f} km".format(lat.degrees, lon.degrees,
                                             height.km))

    if no_save and not previously_downloaded:  # just downloaded
        os.remove(filename)

    return lon.degrees, lat.degrees, height.m


def look_angles(sat_long,
                sat_lat,
                sat_alt,
                rx_long,
                rx_lat,
                rx_height=0,
                ellipsoid="WGS84"):
    """Calculate look angles (elevation, azimuth) and slant range

    Computes the angles relative to a reflector, either active (satellite) or
    passive (radar object). Compute using the rigorous ellipsoidal approach
    discussed in [1].

    Args:
        sat_long       : Subsatellite point's longitude in degrees.
        sat_lat:       : Subsatellite point's geodetic latitude in degrees.
        sat_alt        : Satellite/reflector altitude in meters above the
                         reference ellipsoid.
        rx_long        : Longitude of the receiver station in degrees.
        rx_lat         : Geodetic latitude of the receiver station in degrees.
        rx_height      : Receiver station's orthometric height (height above
                         sea-level) in meters. Defaults to zero.
        ellipsoid      : Reference ellipsoid for the computation, GRS80 or
                         WGS84. Defaults to WGS84.

    Note:
        Positive longitudes are east of the prime meridian and negative
        longitudes are west of the prime meridian. Positive latitudes are north
        of the equator and negative latitudes are south of the equator.

    Returns:
        Tuple with elevation (degrees), azimuth (degrees) and slant range (m).

    """

    # Step 1: PROGRAM INPUTS

    # Convert to radians
    sat_long = radians(sat_long)
    sat_lat = radians(sat_lat)
    rx_long = radians(rx_long)
    rx_lat = radians(rx_lat)

    # Step 2: TRANSFORMATION CURVILINEAR TO CARTESIAN COORDINATES

    # Ellipsoid parameters from GRS80
    f_inv = {
        'GRS80': 298.257222100882711,
        'WGS84': 298.257223563
    }  # reciprocal flattening
    f = 1 / f_inv[ellipsoid]  # flattening
    e_sq = 2 * f - f**2  # eccentricity squared

    # Earth parameters
    R_eq = 6378.137e3  # equatorial radius in meters (see [2])
    r = R_eq + sat_alt  # from the earth's center to the spacecraft

    # Principal radius of curvature in the prime vertical (see the the
    # discussion below Eq 12 in [1]). The computation is also discussed in [2].
    N = R_eq / sqrt(1 - e_sq * sin(rx_lat)**2)

    # Ellipsoidal (geodetic) height of the antenna location:
    Ng = 0  # undulation of the geoid or geoid height
    h = Ng + rx_height

    # Rectangular coordinates of the antenna location, using Eq. 12 from [1]:
    x_p = (N + h) * cos(rx_long) * cos(rx_lat)
    y_p = (N + h) * sin(rx_long) * cos(rx_lat)
    z_p = (N * (1 - e_sq) + h) * sin(rx_lat)

    # Rectangular coordinates of the satellite. See Fig. 5 in [1]:
    x_s = r * cos(sat_long) * cos(sat_lat)
    y_s = r * sin(sat_long) * cos(sat_lat)
    z_s = (N * (1 - e_sq) + sat_alt) * sin(sat_lat)

    # Step 3: SATELLITE COMPONENTS ON LOCAL (x, y, z) COORDINATES
    rect_coor = np.array([x_s, y_s, z_s]) - np.array([x_p, y_p, z_p])

    # rect_coor is a vector starting on the receiver position P and ending on
    # the satellite S (i.e., the topocentric range PS). In other words, it
    # represents the satellite rectangular coordinates referenced to the
    # receiver position. The Euclidean norm of the vector is the slant (or
    # topocentric) range:
    slant_range = np.linalg.norm(rect_coor)

    # Step 4: SATELLITE COMPONENTS ON LOCAL (e, n, u)
    #
    # e-axis points to (geodetic) east; n to (geodetic) north; and u to
    # (geodetic) zenith.

    # Rotation matrix - Eq. 9b from [1]:
    rot_mtx = np.array(
        [[-sin(rx_long), cos(rx_long), 0],
         [
             -sin(rx_lat) * cos(rx_long), -sin(rx_lat) * sin(rx_long),
             cos(rx_lat)
         ],
         [cos(rx_lat) * cos(rx_long),
          cos(rx_lat) * sin(rx_long),
          sin(rx_lat)]])
    # Conversion using Eq. 10 [1]:
    geodetic_coor = np.dot(rot_mtx, rect_coor)

    # Step 5: GEODETIC AZIMUTH AND GEODETIC VERTICAL ANGLE
    e, n, u = geodetic_coor
    azimuth = atan2(e, n)
    vert_angle = atan(u / sqrt(e**2 + n**2))  # elevation

    azimuth_degrees = degrees(azimuth) % 360
    elevation_degrees = degrees(vert_angle)

    util.log_result("Elevation", "{:.2f} degrees".format(elevation_degrees))
    util.log_result("Azimuth", "{:.2f} degrees".format(azimuth_degrees))
    util.log_result("Distance", "{:.2f} km".format(slant_range / 1e3))

    return elevation_degrees, azimuth_degrees, slant_range


def polarization_angle(sat_long, sat_lat, rx_long, rx_lat):
    """Compute the polarization angle (skew) at a given location

    A linearly polarized satellite transmission has the electric field oriented
    at a constant angle relative to a reference plane. For a satellite antenna,
    the reference plane is the equatorial plane. If the polarization is
    horizontal, the electric field is parallel to the equatorial
    plane. Otherwise, when the linear polarization is vertical, it is
    perpendicular to the equatorial plane.

    The earth station antenna feed must have its polarization aligned with the
    polarization plane of the received wave. However, the local vertical
    (normal) plane does not match the equatorial plane due to the curvature of
    the earth. Hence, there is an angle difference between the polarization of
    the signal transmitted by the satellite and the apparent polarization of
    the received signal, which is known as the polarization angle or
    polarization skew (sometimes also referred to as "polarity", "polarity
    skew", "LNB skew", or just "skew").

    When pointing an linearly-polarized earth station antenna, this angle has
    to be taken into account. In contrast, if pointing a circularly-polarized
    antenna, there is no need to compensate for the polarization skew.

    For geostationary satellites, [3] presents the skew formula that follows:

    .. math::

        \\psi = \\tan^{-1} \\left( \\frac{\\sin(L)}{\\tan(l)} \\right)

    where :math:`L` is the relative longitude (difference between the earth
    station's longitude and the satellite longitude) and :math:`l` is the earth
    station's latitude.

    As indicated by the equation, the skew is 0 when the earth station and the
    satellite are at the same longitude (when the relative longitude is zero).

    Note that, for a positive skew value, the LNB must be rotated clockwise,
    whereas, for a negative polarization angle, the LNB must be rotated
    counterclockwise. Furthermore, note that the clockwise and counterclockwise
    directions are relative to the front face (the feed horn) of the LNB. In
    other words, for a positive skew value, someone standing behind the dish
    and facing the satellite in the sky would rotate the LNB clockwise.

    Args:
        sat_long (float): Subsatellite point's longitude in degrees.
        sat_lat  (float): Subsatellite point's geodetic latitude in degrees.
        rx_long  (float): Earth station's longitude in degrees.
        rx_lat   (float): Earth station's geodetic latitude in degrees.

    Returns:
        float: Polarization angle in degrees.

    """
    # Convert to radians
    sat_long = radians(sat_long)
    sat_lat = radians(sat_lat)
    rx_long = radians(rx_long)
    rx_lat = radians(rx_lat)

    # FIXME support computation for satellites in non-geostationary orbit
    if sat_lat != 0:
        logging.warning(
            "The polarization angle is inaccurate for non-geostationary "
            "satellites")

    # Computation based on Eq. (8.22c) from [3]:
    delta_long = rx_long - sat_long
    pol_angle_rad = atan(sin(delta_long) / tan(rx_lat))

    # Return the polarization angle in degrees
    pol_angle = degrees(pol_angle_rad)
    util.log_result("Polarizationa angle", "{:.2f} degrees".format(pol_angle))
    return pol_angle
