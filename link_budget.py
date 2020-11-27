#!/usr/bin/env python3
#
# Satellite Link Budget
#
# To reproduce Example SA8-1 from [1]:
#
# ./link_budget.py --eirp 52 \
#  --freq 12.45e9 \
#  --if-bw 24e6 \
#  --rx-dish-size 0.46 \
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
import argparse
import logging
from math import log10, sqrt, sin, asin, cos, acos, tan, pi, degrees, \
    radians, log2
import util


speed_of_light = 299792458  # in m/s
T0 = 290  # standard room temperature in Kelvin


def calc_eirp(tx_power, tx_dish_gain, tx_dish_size, freq):
    """Compute the EIRP

    Given that EIRP (dB) = Tx Power (dB) + Tx Antenna Gain (dB), we only need
    the two terms. The antenna gain could either be provided directly or
    indirectly by providing the Tx dish size and frequency.

    Args:
        tx_power     : Transmit power feeding the antenna (dBW)
        tx_dish_gain : Transmit antenna gain (dB)
        tx_dish_size : Transmit antenna diameter (m)
        freq         : Transmit signal frequency (Hz)

    """

    if (tx_dish_gain is None):
        tx_dish_gain_db = calc_dish_gain(tx_dish_size, freq)
    else:
        tx_dish_gain_db = tx_dish_gain

    eirp = tx_power + tx_dish_gain_db

    return eirp


def calc_look_angles(sat_long, rx_long, rx_lat, sat_alt=35786e3):
    """Calculate look angles (elevation, azimuth) and slant range

    Computes the angles relative to a reflector, either an active reflector
    (satellite) or a passive reflector (radar object). Assumes the reflector is
    located in the equator (latitude 0).

    Args:
        sat_long   : Longitude of the satellite/reflector in degrees
        rx_long    : Longitude of the receiver station in degrees
        rx_lat     : Latitute of the receiver station in degrees
        sat_alt    : Satellite/reflector altitude (default to geosynchronous
                     altitude)

    Note:
        - Positive longitudes are east, whereas negative longitudes are to the
          west.

    Returns:
        Tuple with elevation (degrees), azimuth (degrees) and slant range (m)

    """
    # Convert to radians
    sat_long = radians(sat_long)
    rx_long = radians(rx_long)
    rx_lat = radians(rx_lat)

    # Constants
    R = 6371e3          # mean radius of the earth in meters
    R_eq = 6378.137e3   # equatorial radius in meters (see [4])
    r = R_eq + sat_alt  # from the earth's center to the spacecraft

    # Eq. (1) from [2]:
    cos_gamma = cos(rx_lat) * cos(sat_long - rx_long)
    gamma = acos(cos_gamma)
    # gamma is the angle between the radius vectors to the Rx location and the
    # sub-satellite point (intersection with the earth's surface of the
    # geocentric radius vector to the satellite). Equation (1) is the cosine of
    # the this angle.

    # Distance between the satellite and the receiver (a.k.a. slant range),
    # from Equation (2) of [2]:
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


def calc_path_loss(d, freq, radar=False, rcs=None, bistatic=None, d_rx=None):
    """Calculate free-space path loss

    Also known as free-space transmission loss.

    This function supports radar mode, in which case the path loss is computed
    by considering both forward and reverse paths to/from the radar object.

    Args:
        d        : Distance in meters between transmitter and receiver
                   or between transmitter and radar object
        freq     : Carrier frequency in Hz
        radar    : Radar mode
        bistatic : Bistatic radar
        rcs      : Radar cross section (RCS)
        d_rx     : Bistatic radar mode only: distance between radar object and
                   receiver that is not collocated with the transmitter

    Notes:

        - The RCS definition repeated in [3] is the following: "the RCS of a
          radar object is the hypothetical area intercepting that amount of
          power which, when scattered isotropically, produces a power density
          at the receiver equal to that from the actual object."

    Returns:
        Path loss in dB

    """
    wavelength = speed_of_light / freq

    # Eq. 8-11 from [1], or Eq. 3.16 from [3]:
    Lfs_one_way_db = 20*log10(4*pi*d/wavelength)

    if (radar):
        assert(rcs is not None)
        assert(bistatic is not None)

        # Radar object gain in dB, equation 3.23 in [3]:
        G_obj_db = 10*log10(4*pi*rcs/(wavelength**2))

        if (bistatic):
            assert(d_rx is not None)
            Lfs_tx_db = Lfs_one_way_db
            Lfs_rx_db = 20*log10(4*pi*d_rx/wavelength)
            # Bistatic radar transmission loss in dB, equation 3.24 in [3]:
            Lfs_db = Lfs_tx_db + Lfs_rx_db - G_obj_db
        else:
            # Monostatic radar transmission loss in dB, equation 3.26 in [3]:
            Lfs_db = 2*Lfs_one_way_db - G_obj_db
    else:
        Lfs_db = Lfs_one_way_db

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
    radius = diameter / 2
    face_area = pi * (radius**2)  # assume circle
    wavelength = speed_of_light / freq

    # See Table 8-4 in [1]:
    gain = 7*face_area/(wavelength**2)
    gain_db = 10*log10(gain)

    logging.info("Dish gain:          {:6.2f} dB".format(gain_db))
    return gain_db


def calc_coax_gain_nf(length_ft, Tl=T0):
    """Compute coaxial RG6 transmission line gain and noise figure

    Args:
        length_ft : Line length in feet
        Tl        : temperature of the line in Kelvin

    Returns:
        Tuple with line loss (dB) and noise figure (dB).

    """
    loss_db_per_ft = 8/100
    loss_db = length_ft * loss_db_per_ft
    loss = util.db_to_abs(loss_db)

    # The noise figure (dB) of a coaxial line is equal to the loss in dB if the
    # physical temperature of the line is equal to T0=290 K. See Equation 8.32a
    # on Example 8-2 in [1]. More generally, any passive two-port element (or
    # attenuator) at room temperature will have this property (noise figure =
    # attenuation in dB), see Equation 4.22 in [3].
    noise_factor = 1 + (Tl/T0)*(loss - 1)
    noise_fig = 10*log10(noise_factor)

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
    F = util.db_to_abs(nfs[0])
    G_prod = 1
    for i, nf in enumerate(nfs[1:]):
        nf_abs = util.db_to_abs(nf)
        G_prod *= util.db_to_abs(gains[i])
        F += (nf_abs - 1) / G_prod

    F_db = 10*log10(F)
    logging.info("Rx noise figure:    {:6.2f} dB".format(F_db))
    return F_db


def noise_fig_to_noise_temp(nf):
    """Convert noise figure to the effective input-noise temperature

    Note that the noise figure is always referenced to a noise source at the
    standard noise temperature of T0 = 290 K. In contrast, the noise
    temperature is independent of the temperature of the noise source.

    Args:
        nf : Noise figure in dB

    Returns:
        Noise temperature in K

    """
    nf_abs = util.db_to_abs(nf)

    # Using Equation 8-30b in [1]:
    Te = T0 * (nf_abs - 1)

    return Te


def noise_temp_to_noise_fig(Te):
    """Convert an effective input-noise temperature to a noise figure in dB

    Args:
        Te : Noise temperature in K

    Returns:
        Noise figure in dB

    """
    # Noise factor
    nf_abs = 1 + Te/T0
    # Return the noise figure
    return 10*log10(nf_abs)


def calc_rx_sys_noise_temp(Tar, Te):
    """Compute the receiver system noise temperature in dB

    The receiver noise temperature is the sum of the effective input-noise
    temperature (Te) of the entire receiver seen as a blackbox and the antenna
    noise temperature (Tar). The Te term represents the noise introduced by the
    cascaded linear components (e.g., the LNB, the coax line and e.g., the
    radio interface) of the receiver. The Tar component, in turn, is the noise
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

    As explained in [3], around equation 4.39, this is just a convenient choice
    in terms of where the effective input-noise temperature is observed. Note
    that an equivalent (or effective) input noise temperature represents the
    thermodynamic temperature of a noisy resistor, connected to the input of a
    noiseless two-port element, which gives the same output noise power as the
    noisy but otherwise equivalent two-port element, with an ideal noiseless
    source at its input [3]. Hence, if we group the entire receiver into a
    single equivalent two-port element, the Te term represents the noise
    generated by the entire receiver. This noise has power k*Te*B. When
    combined to the noise collected by the antenna, one obtaines the total
    system noise temperature.

    Args:
        Tar : Antenna noise temperature in K
        Te  : Effective input-noise temperature in K

    Returns:
        The receiver system noise temperature in dB

    Note: here we return the noise temperature directly in dB. The dB form will
    be directly subtracted in the C/N computation. See Equation 8-43 in [1].

    """

    # Equation 8-41 from [1], or 4.39 from [3]:
    Tsyst = Tar + Te
    logging.info("System noise temp:  {:6.2f} K".format(Tsyst))

    Tsyst_db = 10*log10(Tsyst)

    return Tsyst_db


def calc_cnr(eirp_db, path_loss_db, rx_ant_gain_db, T_sys_db, bw):
    """Compute the carrier-to-noise ratio (CNR) in dB

    Args:
        eirp_db        : EIRP in dBW
        path_loss_db   : Free-space path loss in dB
        rx_ant_gain_db : Receiver antenna gain in dB
        T_sys_db       : Receiver system noise temperature in dB
        bw             : Nominal signal bandwidth

    Returns:
        CNR (also known as C/N) in dB

    """
    # According to Equation 8-40 from [1], the noise power is given by N =
    # k*Tsyst*bw, where k is Boltzmann’s constant, Tsyst is the receiver system
    # noise temperature (in absolute units) and bw is the IF equivalent
    # bandwidth. On the C/N computation in dB, N is in the denominator, and
    # hence we can simply subtract k_db and B_db to compute C/N, see Equation
    # 8-43 from [1].
    k_db = -228.6  # Boltzmann’s constant (of 1.38e-23) in dB
    bw_db = 10*log10(bw)

    # The received power level at the antenna terminals is of interest, so
    # print it it out:
    P_rx_dbw = eirp_db - path_loss_db + rx_ant_gain_db
    P_rx_dbm = P_rx_dbw + 30
    logging.info("Rx Power:           {:6.2f} dBm".format(P_rx_dbm))

    # The ratio between the Rx antenna gain and the receiver noise temperature,
    # usually known as G/T, is also a metric of interest, so log it:
    g_over_t_db = rx_ant_gain_db - T_sys_db
    logging.info("(G/T):              {:6.2f} dB/K".format(g_over_t_db))

    # C/N, as computed in Equation 8-43 from [1]:
    cnr_db = eirp_db - path_loss_db + g_over_t_db - k_db - bw_db
    logging.info("(C/N):              {:6.2f} dB".format(cnr_db))

    return cnr_db


def calc_capacity(snr_db, bw):
    """Compute the channel capacity in bps

    Args:
        snr_db : signal-to-noise ratio in dB
        bw     : nominal bandwidth

    """
    snr = util.db_to_abs(snr_db)
    c = bw * log2(1 + snr)
    logging.info("Capacity:           {}".format(util.format_rate(c)))


def parser():
    parser = argparse.ArgumentParser(description="Link Budget")

    tx_pwr_group = parser.add_mutually_exclusive_group(required=True)
    tx_pwr_group.add_argument(
        '--eirp',
        type=float,
        help='EIRP in dBW.'
    )
    tx_pwr_group.add_argument(
        '--tx-power',
        type=float,
        help='Tx power feeding the antenna in dBW.'
    )
    tx_dish_group = parser.add_mutually_exclusive_group()
    tx_dish_group.add_argument(
        '--tx-dish-size',
        type=float,
        help='Diameter in meters of the parabolic antenna used for '
        'transmission. Used when the power is specified through option '
        '--tx-power'
    )
    tx_dish_group.add_argument(
        '--tx-dish-gain',
        type=float,
        help='Gain in dB of the parabolic antenna used for transmission. Used '
        'when the power is specified through option --tx-power'
    )
    parser.add_argument(
        '--freq',
        required=True,
        type=float,
        help='Downlink carrier frequency in Hz for satellite signals or '
        'simply the signal frequency in Hz for radar (passively reflected) '
        'signals.'
    )
    parser.add_argument(
        '--if-bw',
        required=True,
        type=float,
        help='IF bandwidth in Hz.'
    )
    rx_dish_group = parser.add_mutually_exclusive_group(required=True)
    rx_dish_group.add_argument(
        '--rx-dish-size',
        type=float,
        help='Parabolic antenna (dish) diameter in m.'
    )
    rx_dish_group.add_argument(
        '--rx-dish-gain',
        type=float,
        help='Parabolic antenna (dish) gain in dBi.'
    )
    parser.add_argument(
        '--antenna-noise-temp',
        required=True,
        type=float,
        help='Receive antenna\'s noise temperature in K.'
    )
    lnb_noise_group = parser.add_mutually_exclusive_group(required=True)
    lnb_noise_group.add_argument(
        '--lnb-noise-fig',
        type=float,
        help='LNB\'s noise figure in dB.'
    )
    lnb_noise_group.add_argument(
        '--lnb-noise-temp',
        type=float,
        help='LNB\'s noise temperature in K.'
    )
    parser.add_argument(
        '--lnb-gain',
        required=True,
        type=float,
        help='LNB\'s gain.'
    )
    parser.add_argument(
        '--coax-length',
        required=True,
        type=float,
        help='Length of the coaxial transmission line between the LNB and the '
        'receiver in ft.'
    )
    parser.add_argument(
        '--rx-noise-fig',
        required=True,
        type=float,
        help='Receiver\'s noise figure in dB.'
    )
    parser.add_argument(
        '--sat-long',
        required=True,
        type=float,
        help='Satellite\'s longitude. Negative to the West and positive to '
        'the East'
    )
    parser.add_argument(
        '--sat-lat',
        type=float,
        help='Satellite\'s latitude. Positive to the North and negative to '
        'the South'
    )
    parser.add_argument(
        '--rx-long',
        required=True,
        type=float,
        help='Receive station\'s longitude. Negative to the West and positive '
        'to the East'
    )
    parser.add_argument(
        '--rx-lat',
        required=True,
        type=float,
        help='Receive station\'s latitude. Positive to the North and negative '
        'to the South'
    )
    radar_p = parser.add_argument_group('radar options')
    radar_p.add_argument(
        '--radar',
        default=False,
        action='store_true',
        help='Activate radar mode, so that the link budget considers the '
        'pathloss to and back from object'
    )
    radar_p.add_argument(
        '--radar-alt',
        type=float,
        help='Altitude of the radar object'
    )
    radar_p.add_argument(
        '--radar-cross-section',
        type=float,
        help='Radar cross section of the radar object'
    )
    radar_p.add_argument(
        '--radar-bistatic',
        default=False,
        action='store_true',
        help='Bistatic radar scenario, i.e., radar transmitter and receiver '
        'are not collocated'
    )
    args = parser.parse_args()

    if (args.tx_power and args.tx_dish_size is None and
            args.tx_dish_gain is None):
        parser.error("Define either --tx-dish-size or --tx-dish-gain  "
                     "using option --tx-power")

    if (args.radar):
        if (args.radar_alt is None):
            parser.error("Argument --radar-alt is required in radar mode "
                         "(--radar)")
        if (args.radar_cross_section is None):
            parser.error("Argument --radar-cross-section is required in radar "
                         "mode (--radar)")
    return args


def main():
    logging.basicConfig(level=logging.INFO)
    args = parser()

    sat_alt = 35786e3 if not args.radar else args.radar_alt

    elevation, azimuth, slant_range = calc_look_angles(args.sat_long,
                                                       args.rx_long,
                                                       args.rx_lat,
                                                       sat_alt)

    # Compute the EIRP
    if (args.eirp is None):
        eirp = calc_eirp(args.tx_power, args.tx_dish_gain, args.tx_dish_size,
                         args.freq)
        logging.info("Tx Power:           {:6.2f} kW".format(
            util.db_to_abs(args.tx_power)/1e3))
    else:
        eirp = args.eirp

    logging.info("EIRP:               {:6.2f} dBW ({:6.2f} kW)".format(
        eirp, util.db_to_abs(eirp)/1e3))

    path_loss_db = calc_path_loss(slant_range, args.freq, args.radar,
                                  args.radar_cross_section,
                                  args.radar_bistatic)
    # TODO support bistatic radar. Add distance from radar object to rx
    # station.

    if (args.rx_dish_gain is None):
        dish_gain_db = calc_dish_gain(args.rx_dish_size, args.freq)
    else:
        dish_gain_db = args.rx_dish_gain

    coax_loss_db, coax_noise_fig_db = calc_coax_gain_nf(args.coax_length)

    if (args.lnb_noise_fig is None):
        lnb_noise_fig = noise_temp_to_noise_fig(args.lnb_noise_temp)
    else:
        lnb_noise_fig = args.lnb_noise_fig

    logging.info("LNB noise figure:   {:6.2f} dB".format(lnb_noise_fig))

    noise_fig_db = calc_total_noise_figure(
        [lnb_noise_fig, coax_noise_fig_db, args.rx_noise_fig],
        [args.lnb_gain, -coax_loss_db]
    )

    effective_input_noise_temp = noise_fig_to_noise_temp(noise_fig_db)

    logging.info("Antenna noise temp: {:6.2f} K".format(
        args.antenna_noise_temp))
    logging.info("Input-noise temp:   {:6.2f} K".format(
        effective_input_noise_temp))

    T_syst_db = calc_rx_sys_noise_temp(args.antenna_noise_temp,
                                       effective_input_noise_temp)

    cnr = calc_cnr(eirp, path_loss_db, dish_gain_db, T_syst_db, args.if_bw)

    calc_capacity(cnr, args.if_bw)


if __name__ == '__main__':
    main()
