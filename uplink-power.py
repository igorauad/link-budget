#!/usr/bin/env python3
import argparse
import logging

from linkbudget import calc, constants, pointing, util
from linkbudget.antenna import Antenna


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--freq',
                   required=True,
                   type=float,
                   help='Uplink carrier frequency in Hz')
    p.add_argument('--bw',
                   required=True,
                   type=float,
                   help='Noise bandwidth in Hz.')
    p.add_argument('--ul-cnr',
                   required=True,
                   type=float,
                   help='Target UL C/N ratio')
    p.add_argument('--ul-dish-size',
                   required=True,
                   type=float,
                   help='Diameter in meters of the uplink antenna.')
    p.add_argument('--ul-dish-efficiency',
                   type=float,
                   default=0.68,
                   help='Aperture efficiency of the parabolic antenna.')
    p.add_argument('--ul-losses',
                   type=float,
                   default=0,
                   help='Miscellaneous UL losses in dB')
    p.add_argument('--ul-long',
                   required=True,
                   type=float,
                   help='UL station\'s longitude.')
    p.add_argument('--ul-lat',
                   required=True,
                   type=float,
                   help='UL station\'s latitude.')
    p.add_argument('--sat-long',
                   required=True,
                   type=float,
                   help='Satellite\'s longitude.')
    tp_info = p.add_mutually_exclusive_group(required=True)
    tp_info.add_argument('--tp-noise-temp',
                         type=float,
                         help='Transponder system noise temperature.')
    tp_info.add_argument('--tp-gt',
                         type=float,
                         help='Tranponder\'s G/T ratio in dB/K')
    p.add_argument('--upc-range',
                   required=True,
                   type=float,
                   help='Target UPC dynamic range in dB.')
    p.add_argument(
        '--ul-contour',
        type=float,
        default=0.0,
        help='Contour line where the UL station is located. A negative value '
        'in dB representing the off-axis gain dropoff of the satellite '
        'antenna in the direction of the UL station relative to on-axis '
        '(boresight) gain. Taken into account only when option '
        '--tp-noise-temp is provided instead of --tp-gt.')

    args = p.parse_args()

    if (args.ul_contour > 0):
        p.error("--ul-contour must be a negative value in dB.")

    return args


def configure_logging():
    """Configure the logging format"""
    logging_fmt = "%(levelname)s %(message)s"
    logging.basicConfig(level=logging.INFO, format=logging_fmt)
    # Keep the level name except if it's an INFO message
    logging.addLevelName(logging.INFO, '')


def main():
    args = parse_args()
    configure_logging()

    _, _, slant_range_m = pointing.look_angles(args.sat_long, args.ul_long,
                                               args.ul_lat,
                                               constants.GEOSYNC_ORBIT)
    ul_antenna = Antenna(freq=args.freq,
                         diameter=args.ul_dish_size,
                         efficiency=args.ul_dish_efficiency,
                         label="UL Ant")

    path_loss_db = calc.path_loss(slant_range_m, args.freq)
    loss_db = path_loss_db + args.ul_losses

    # If the transponder's G/T ratio is known, compute the exact power
    # requirement. Otherwise, at least the system noise temperature is
    # required, and then several values of G can be assumed.
    if (args.tp_gt is not None):
        k_db = -228.6  # Boltzmann’s constant (of 1.38e-23) in dB
        bw_db = util.lin_to_db(args.bw)
        # Assume:
        # CNR = Pt + Gt + Gr/T - k - B - L
        #
        # Hence,
        # Pt = CNR - (Gt + Gr/T - k - B - L)
        #    = CNR - Gt - Gr/T + k + B + L
        Pt_dbw = args.ul_cnr - ul_antenna.gain_db - args.tp_gt + \
            k_db + bw_db + loss_db
        # Include also some extra available power for UPC
        Pt_dbw += args.upc_range

        util.log_result(
            "Required Tx Power",
            "{:.1f} dBW ({:.2f} W)".format(Pt_dbw, util.db_to_lin(Pt_dbw)))
        return

    tp_noise_temp_db = util.lin_to_db(args.tp_noise_temp)
    N_dbw = calc.noise_power(tp_noise_temp_db, args.bw)

    P_rx_dbw = N_dbw + args.ul_cnr
    util.log_result("Rx Power at TP Input", "{:.2f} dBW".format(P_rx_dbw))

    sat_antennas = [0.45, 0.6, 0.9]
    for idx, dish_size in enumerate(sat_antennas):
        util.log_result("---", "---")
        util.log_result("Satellite Antenna {}".format(idx + 1),
                        "{:.2f} m".format(dish_size))
        sat_antenna = Antenna(freq=args.freq,
                              diameter=dish_size,
                              efficiency=0.68,
                              label="Sat Ant")
        sat_ant_gain_off_axis = sat_antenna.gain_db + args.ul_contour
        calc.g_over_t(sat_ant_gain_off_axis, tp_noise_temp_db)
        Pt_dbw = P_rx_dbw + loss_db - ul_antenna.gain_db - \
            sat_ant_gain_off_axis + args.upc_range
        util.log_result(
            "Required Tx Power",
            "{:.1f} dBW ({:.2f} W)".format(Pt_dbw, util.db_to_lin(Pt_dbw)))


if __name__ == '__main__':
    main()
