#!/usr/bin/env python3
#
# Satellite Link Budget
#
# References:
#
# [1] Couch, Leon W.. Digital & Analog Communication Systems.
import argparse, math, logging
from math import log10, sin, asin, cos, acos, tan, atan, pi, degrees, radians


speed_of_light = 299792458 # in m/s


def calc_look_angles(sat_long, rx_long, rx_lat):
    """ Calculate look angles (elevation, azimuth) and slant range

    Args:
        sat_long : Longitude of the geosynchronous satellite in degrees
        rx_long  : Longitude of the receiver station in degrees
        rx_lat   : Latitute of the receiver station in degrees

    Returns:
        Tuple with elevation (degrees), azimuth (degrees) and slant range (m)

    """
    # Convert to radians
    sat_long  = radians(sat_long)
    rx_long   = radians(rx_long)
    rx_lat    = radians(rx_lat)

    long_diff = sat_long - rx_long
    R         = 3963   # Earth radius in statute miles
    h         = 22242  # Geosync satellite height in statute miles
    # Eq. 8-48b from [1]:
    beta = acos(cos(rx_lat) * cos(long_diff))
    # Elevation: Eq. 8-48a from [1]:
    E = math.atan((1/math.tan(beta)) - (R/((R + h)*sin(beta))))
    # Distance between the satellite and the receiver (slant range), from
    # Eq. 8-49 in [1]:
    d = math.sqrt((R+h)**2 + R**2 - 2*R*(R+h)*math.cos(beta))
    # Azimuth, from Eq. 8-50 in [1]:
    A = 2*pi - math.acos(-math.tan(rx_lat)/math.tan(beta))

    # Convert back to degrees and m
    E = degrees(E)
    A = degrees(A)
    d = d * 1609.3472

    logging.info("Elevation:  {:6.2f} degrees".format(E))
    logging.info("Azimuth:    {:6.2f} degrees".format(A))
    logging.info("Distance:   {:8.2f} km".format(d/1e3))

    return E, A, d


def calc_path_loss(d, freq):
    """ Calculate free-space path loss

    Args:
        d    : Distance in meters
        freq : Carrier frequency in Hz

    Returns:
        Path loss in dB

    """
    wavelength = speed_of_light / freq

    # Eq. 8-11 from [1]:
    Lfs_db = 20*log10(4*pi*d/wavelength)

    logging.info("Path loss:  {:6.2f} dB".format(Lfs_db))
    return Lfs_db


def calc_dish_gain(diameter, freq):
    """ Calculate parabolic dish gain

    Args:
        diameter : Diameter in m
        freq     : Frequency of interest in Hz

    Returns:
        Gain in dB
    """
    radius     = diameter / 2
    face_area  = pi * (radius**2) # assume circle
    wavelength = speed_of_light / freq

    # See Table 8-4 in [1]:
    gain    = 7*face_area/(wavelength**2)
    gain_db = 10*log10(gain)

    logging.info("Dish gain:  {:6.2f} dB".format(gain_db))
    return gain_db


def parser():
    parser = argparse.ArgumentParser(description="Link Budget")
    parser.add_argument('--eirp',
                        default=36,
                        type=float,
                        help='EIRP in dBw.')
    parser.add_argument('--f-dl',
                        default=4e9,
                        type=float,
                        help='Downlink carrier frequency in Hz.')
    parser.add_argument('--dish-size',
                        default=3.05,
                        type=float,
                        help='Parabolic antenna (dish) diameter in m.')
    parser.add_argument('--sat-long',
                        default=134,
                        type=float,
                        help='Satellite\'s longitude.')
    parser.add_argument('--sat-lat',
                        default=0,
                        type=float,
                        help='Satellite\'s latitude.')
    parser.add_argument('--rx-long',
                        default=77,
                        type=float,
                        help='Receive station\'s longitude.')
    parser.add_argument('--rx-lat',
                        default=38.8,
                        type=float,
                        help='Receive station\'s latitude.')
    args = parser.parse_args()
    return args


def main():
    logging.basicConfig(level=logging.INFO)
    args = parser()
    elevation, azimuth, slant_range = calc_look_angles(args.sat_long,
                                                       args.rx_long,
                                                       args.rx_lat)
    path_loss_db = calc_path_loss(slant_range, args.f_dl)
    dish_gain_db = calc_dish_gain(args.dish_size, args.f_dl)


if __name__ == '__main__':
    main()
