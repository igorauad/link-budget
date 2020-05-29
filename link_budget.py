#!/usr/bin/env python3
#
# Satellite Link Budget
#
# To reproduce Example SA8-1 from [1]:
#
# ./link_budget.py --eirp 52 \
#  --f-dl 12.45e9 \
#  --if-bw 24e6 \
#  --dish-size 0.46 \
#  --antenna-noise-temp 20 \
#  --lnb-noise-fig 0.6 \
#  --lnb-gain 40 \
#  --coax-length 110 \
#  --rx-noise-fig 10 \
#  --sat-long -101 \
#  --rx-long -82.43 \
#  --rx-lat 29.71
#
# References:
#
# [1] Couch, Leon W.. Digital & Analog Communication Systems.
#
# [2] https://www.ngs.noaa.gov/CORS/Articles/SolerEisemannJSE.pdf
#
# [3] Lindgren, M. (2015). A 1296 MHz Earth–Moon–Earth Communication System
# (Master's thesis).
#
# [4] https://en.wikipedia.org/wiki/Earth_radius
#
import argparse, math, logging
from math import log10, sqrt, sin, asin, cos, acos, tan, atan, pi, degrees, \
    radians


speed_of_light = 299792458 # in m/s


def _db_to_abs(val_db):
    return 10**(val_db/10)


def calc_look_angles(sat_long, rx_long, rx_lat):
    """Calculate look angles (elevation, azimuth) and slant range

    Args:
        sat_long : Longitude of the geosynchronous satellite in degrees
        rx_long  : Longitude of the receiver station in degrees
        rx_lat   : Latitute of the receiver station in degrees

    Note:
        - Positive longitudes are east, whereas negative longitudes are to the
          west.

    Returns:
        Tuple with elevation (degrees), azimuth (degrees) and slant range (m)

    """
    # Convert to radians
    sat_long  = radians(sat_long)
    rx_long   = radians(rx_long)
    rx_lat    = radians(rx_lat)

    # Constants
    R         = 6371e3           # mean radius of the earth in meters
    R_eq      = 6378.137e3       # equatorial radius in meters (see [4])
    geo_orbit = 35786e3
    r         = R_eq + geo_orbit # from the earth's center to the spacecraft

    # Eq. (1) from [2]:
    cos_gamma = cos(rx_lat) * cos(sat_long - rx_long)
    gamma     = acos(cos_gamma)
    # gamma is the angle between the radius vectors to the Rx location and the
    # sub-satellite point (intersection with the earth's surface of the
    # geocentric radius vector to the satellite). Equation (1) is the cosine of
    # the this angle.

    # Distance between the satellite and the receiver (a.k.a. slant range), from
    # Equation (2) of [2]:
    d = r * sqrt(1 + (R/r)**2 - 2*(R/r)*cos_gamma)

    # Zenith distance, Equation (4) from [2]:
    z = asin((r/d)*sin(gamma))

    # Elevation:
    v = 90 - degrees(z)

    # Angle of Equation (6) from [2]:
    beta = degrees(acos(tan(rx_lat)/tan(gamma)))

    # Azimuth:
    if (rx_lat > 0):
        # Rx is north of the satellite
        if (sat_long < rx_long):
            # Satellite to SW
            alpha = 180 + beta
        else:
            # Satellite to SE
            alpha = 180 - beta
    else:
        # Rx is south of the satellite
        if (sat_long < rx_long):
            # Satellite to NW
            alpha = 360 - beta
        else:
            # Satellite to NE
            alpha = beta

    logging.info("Elevation:          {:6.2f} degrees".format(v))
    logging.info("Azimuth:            {:6.2f} degrees".format(alpha))
    logging.info("Distance:           {:8.2f} km".format(d/1e3))

    return v, alpha, d


def calc_path_loss(d, freq):
    """Calculate free-space path loss

    Args:
        d    : Distance in meters
        freq : Carrier frequency in Hz

    Returns:
        Path loss in dB

    """
    wavelength = speed_of_light / freq

    # Eq. 8-11 from [1]:
    Lfs_db = 20*log10(4*pi*d/wavelength)

    logging.info("Path loss:          {:6.2f} dB".format(Lfs_db))
    return Lfs_db


def calc_dish_gain(diameter, freq):
    """Calculate parabolic dish gain

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

    logging.info("Dish gain:          {:6.2f} dB".format(gain_db))
    return gain_db


def calc_coax_gain_nf(length_ft):
    """Compute coaxial RG6 transmission line gain and noise figure

    Args:
        length_ft : Line length in feet

    Returns:
        Tuple with line loss (dB) and noise figure (dB).

    """
    loss_db_per_ft = 8/100
    loss_db        = length_ft * loss_db_per_ft

    # The noise figure (dB) of the coaxial line is equal to the loss in dB. See
    # Example 8-2 in [1]:
    noise_fig = loss_db

    logging.info("Coax loss:          {:6.2f} dB".format(loss_db))
    logging.info("Coax noise figure:  {:6.2f} dB".format(noise_fig))

    return loss_db, noise_fig


def calc_total_noise_figure(nfs, gains):
    """Calculate the overall noise figure of the receive system

    Args:
        nfs   : List with noise figures (dB) of cascaded linear devices, in
                order.
        gains : List with gains (dB) of cascaded linear devices, in order.

    Note: The list of gains should not include the gain of the last device in
    the chain, as it is irrelevant for the overall noise figure computation.

    Returns:
        The overall noise figure in dB

    """

    assert(len(nfs) > 0)
    assert(len(gains) > 0)
    assert(len(gains) == len(nfs) - 1)

    if (len(nfs) == 1):
        return nfs[0]

    # Implement Equation 8-34 from [1]:
    F      = _db_to_abs(nfs[0])
    G_prod = 1
    for i, nf in enumerate(nfs[1:]):
        nf_abs = _db_to_abs(nf)
        G_prod *= _db_to_abs(gains[i])
        F      += (nf_abs - 1) / G_prod

    F_db = 10*log10(F)
    logging.info("Noise figure:       {:6.2f} dB".format(F_db))
    return F_db


def noise_fig_to_noise_temp(nf):
    """Convert noise figure to the effective input-noise temperature

    Args:
        nf : Noise figure in dB

    Returns:
        Noise temperature in K

    """
    T0     = 290 # standard room temperature in K
    nf_abs = _db_to_abs(nf)

    # Using Equation 8-30b:
    Te = T0 * (nf_abs - 1)
    logging.info("Input-noise temp:   {:6.2f} K".format(Te))
    return Te


def calc_rx_sys_noise_temp(Tar, Te):
    """Compute the receiver system noise temperature in dB

    The receiver noise temperature is the sum of the effective input-noise
    temperature (Te) and the antenna noise temperature (Tar). The Te component
    is the noise introduced by the cascaded linear components (e.g., the LNB,
    the coax line and the receiver). The Tar component, in turn, is the noise
    captured by the antenna due to received cosmic noise and Earth blackbody
    radiation. The simplified model is as follows:

    Rx Antenna (Tar) -----> Sum ----> Noise-free Gain Stage ---> Detector
                       ^
                       |
                       |
                Receiver Noise (Te)

    Note that this is peculiar because we don't combine the antenna as another
    cascaded device to the receiver. Instead, we sum the cascaded devices with
    the antenna, based on this model. See Figure 8-24 in [1].

    Args:
        Tar : Antenna noise temperature in K
        Te  : Effective input-noise temperature in K

    Returns:
        The receiver system noise temperature in dB

    Note: here we return the noise temperature directly in dB. The dB form will
    be directly subtracted in the C/N computation. See Equation 8-43 in [1].

    """

    # Equation 8-41 from [1]:
    Tsyst = Tar + Te

    Tsyst_db = 10*log10(Tsyst)

    logging.info("System noise temp:  {:6.2f} dB".format(Tsyst_db))

    return Tsyst_db


def calc_c_to_n(eirp_db, path_loss_db, rx_ant_gain_db, T_sys_db, bw):
    """Compute the carrier-to-noise (C/N) ratio in dB

    Args:
        eirp_db        : EIRP in dBw
        path_loss_db   : Free-space path loss in dB
        rx_ant_gain_db : Receiver antenna gain in dB
        T_sys_db       : Receiver system noise temperature in dB
        bw             :

    """
    # According to Equation 8-40 from [1], the noise power is given by N =
    # k*Tsyst*bw, where k is Boltzmann’s constant, Tsyst is the receiver system
    # noise temperature (in absolute units) and bw is the IF equivalent
    # bandwidth. On the C/N computation in dB, N is in the denominator, and
    # hence we can simply subtract k_db and B_db to compute C/N, see Equation
    # 8-43 from [1].
    k_db  = -228.6 # Negative in dB because the Boltzmann’s constant is 1.38e-23
    bw_db = 10*log10(bw)

    # The ration between the Rx antenna gain and the receiver noise temperature
    # is a metric of interest, so log it:
    g_over_t_db = rx_ant_gain_db - T_sys_db
    logging.info("(Gar/Tsyst):        {:6.2f} dB/K".format(g_over_t_db))

    # Equation 8-43 from [1]:
    C_N_db = eirp_db - path_loss_db + g_over_t_db - k_db - bw_db
    logging.info("(C/N):              {:6.2f} dB".format(C_N_db))


def parser():
    parser = argparse.ArgumentParser(description="Link Budget")
    parser.add_argument('--eirp',
                        default=52,
                        type=float,
                        help='EIRP in dBw.')
    parser.add_argument('--f-dl',
                        default=12.45e9,
                        type=float,
                        help='Downlink carrier frequency in Hz.')
    parser.add_argument('--if-bw',
                        default=24e6,
                        type=float,
                        help='IF bandwidth in Hz.')
    parser.add_argument('--dish-size',
                        default=0.46,
                        type=float,
                        help='Parabolic antenna (dish) diameter in m.')
    parser.add_argument('--antenna-noise-temp',
                        default=20,
                        type=float,
                        help='Receive antenna\'s noise temperature in K.')
    parser.add_argument('--lnb-noise-fig',
                        default=0.6,
                        type=float,
                        help='LNB\'s noise figure in dB.')
    parser.add_argument('--lnb-gain',
                        default=40,
                        type=float,
                        help='LNB\'s gain.')
    parser.add_argument('--coax-length',
                        default=110,
                        type=float,
                        help='Length of the coaxial transmission line between '
                        'the LNB and the receiver in ft.')
    parser.add_argument('--rx-noise-fig',
                        default=10,
                        type=float,
                        help='Receiver\'s noise figure in dB.')
    parser.add_argument('--sat-long',
                        default=101,
                        type=float,
                        help='Satellite\'s longitude. Negative to the West and '
                        'positive to the East')
    parser.add_argument('--sat-lat',
                        default=0,
                        type=float,
                        help='Satellite\'s latitude. Positive to the North '
                        'and negative to the South')
    parser.add_argument('--rx-long',
                        default=82.43,
                        type=float,
                        help='Receive station\'s longitude. Negative to the '
                        'West and positive to the East')
    parser.add_argument('--rx-lat',
                        default=29.71,
                        type=float,
                        help='Receive station\'s latitude. Positive to the '
                        'North and negative to the South')
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

    coax_loss_db, coax_noise_fig_db = calc_coax_gain_nf(args.coax_length)

    noise_fig_db = calc_total_noise_figure(
        [args.lnb_noise_fig, coax_noise_fig_db, args.rx_noise_fig],
        [args.lnb_gain, -coax_loss_db]
    )

    effective_input_noise_temp = noise_fig_to_noise_temp(noise_fig_db)

    T_syst_db = calc_rx_sys_noise_temp(args.antenna_noise_temp,
                                       effective_input_noise_temp)

    calc_c_to_n(args.eirp, path_loss_db, dish_gain_db, T_syst_db, args.if_bw)


if __name__ == '__main__':
    main()
