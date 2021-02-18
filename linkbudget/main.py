"""Link budget analysis"""
import json
import logging
import argparse
from . import calc, pointing, util


__version__ = "0.1.2"


def get_parser():
    """Command-line arguments"""
    parser = argparse.ArgumentParser(
        description="Link Budget Calculator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--json',
        action='store_true',
        help='Print results in JSON format.'
    )
    tx_pwr_group = parser.add_mutually_exclusive_group(required=True)
    tx_pwr_group.add_argument(
        '--eirp',
        type=float,
        help='EIRP in dBW.'
    )
    tx_pwr_group.add_argument(
        '--tx-power',
        type=float,
        help='Power feeding the Tx antenna in dBW.'
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
        help='Gain in dBi of the parabolic antenna used for transmission. '
        'Used when the power is specified through option --tx-power'
    )
    parser.add_argument(
        '--tx-dish-efficiency',
        type=float,
        default=0.56,
        help='Aperture efficiency of the parabolic antenna used for '
        'transmission. Considered when the dish is specified by size.'
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
        '--rx-dish-efficiency',
        type=float,
        default=0.56,
        help='Aperture efficiency of the parabolic antenna used for '
        'reception. Considered when the dish is specified by size.'
    )
    sky_noise_group = parser.add_mutually_exclusive_group(required=True)
    sky_noise_group.add_argument(
        '--antenna-noise-temp',
        type=float,
        help='Receive antenna\'s noise temperature in K.'
    )
    sky_noise_group.add_argument(
        '--atmospheric-loss',
        type=float,
        help='Attenuation in dB experienced through the atmosphere. It should '
        'always include the clear air attenuation, and it could also include '
        'other effects such as rain attenuation'
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
        '--mispointing-loss',
        type=float,
        default=0,
        help='Loss in dB due to antenna mispointing'
    )
    pos_p = parser.add_argument_group('sat/rx positioning options')
    pos_p.add_argument(
        '--sat-long',
        type=float,
        help='Satellite\'s longitude. Negative to the West and positive to '
        'the East'
    )
    pos_p.add_argument(
        '--rx-long',
        type=float,
        help='Rx station\'s longitude. Negative to the West and positive '
        'to the East'
    )
    pos_p.add_argument(
        '--rx-lat',
        type=float,
        help='Rx station\'s latitude. Positive to the North and negative '
        'to the South'
    )
    pos_p.add_argument(
        '--slant-range',
        type=float,
        help='Slant path length in km between the Rx station and the '
        'satellite or reflector'
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
    return parser


def validate(parser, args):
    """Validate command-line arguments"""
    if (args.tx_power and args.tx_dish_size is None and
            args.tx_dish_gain is None):
        parser.error("Define either --tx-dish-size or --tx-dish-gain  "
                     "using option --tx-power")

    if (args.radar):
        if (args.slant_range is None and args.radar_alt is None):
            parser.error("Argument --radar-alt is required in radar mode "
                         "if the --slant-range argument is not provided")
        if (args.radar_cross_section is None):
            parser.error("Argument --radar-cross-section is required in radar "
                         "mode (--radar)")

    # Mutual exclusion between "--slant-range" and the rx/sat positioning args
    pos_args = [args.sat_long, args.rx_long, args.rx_lat]
    pos_arg_labels = ['--sat-long', '--rx-long', '--rx-lat']
    missing_pos_args = []
    defined_pos_args = []
    for arg, label in zip(pos_args, pos_arg_labels):
        if (arg is None):
            missing_pos_args.append(label)
        else:
            defined_pos_args.append(label)
    if (args.slant_range is None and len(missing_pos_args) > 0):
        parser.error("the following arguments are required: {}".format(
            ", ".join(missing_pos_args)
        ))
    elif (args.slant_range is not None and len(defined_pos_args) > 0):
        parser.error("argument{} {}: not allowed with argument "
                     "--slant-range".format(
                         "s" if len(defined_pos_args) > 1 else "",
                         ", ".join(defined_pos_args))
                     )


def analyze(args, verbose=False):
    """Main link budget analysis

    Args:
        args : Populated argparse namespace object.
        verbose : Verbose mode.

    Returns:
        Dictionary with the main link budget results.

    """
    if (verbose and not args.json):
        logging.basicConfig(level=logging.INFO)

    # -------- Look angles --------
    if (args.slant_range is None):
        # Satellite altitude
        sat_alt = 35786e3 if not args.radar else args.radar_alt

        # Look angles
        elevation, azimuth, slant_range = pointing.look_angles(
            args.sat_long, args.rx_long, args.rx_lat, sat_alt)
    else:
        elevation = azimuth = None
        slant_range = args.slant_range * 1e3  # km to m

    # -------- EIRP --------
    if (args.eirp is None):
        if args.tx_dish_gain is None:
            tx_gain = calc.dish_gain(args.tx_dish_size, args.freq,
                                     args.tx_dish_efficiency)
            logging.info("Tx dish gain:       {:6.2f} dB".format(tx_gain))
        else:
            tx_gain = args.tx_dish_gain
        eirp = calc.eirp(args.tx_power, tx_gain)
        logging.info("Tx Power:           {:6.2f} kW".format(
            util.db_to_abs(args.tx_power)/1e3))
    else:
        eirp = args.eirp

    logging.info("EIRP:               {:6.2f} dBW".format(eirp))

    # -------- Path loss --------
    path_loss_db = calc.path_loss(slant_range, args.freq, args.radar,
                                  args.radar_cross_section,
                                  args.radar_bistatic)
    # TODO support bistatic radar. Add distance from radar object to rx
    # station.

    # When defined, the atmospheric loss defines the antenna noise temperature
    # on reception. It also adds to the free-space path loss. On radar systems,
    # assume the atmospheric loss is experienced twice.
    if (args.atmospheric_loss is None):
        atmospheric_loss_db = 0
    elif (args.radar):
        atmospheric_loss_db = 2 * args.atmospheric_loss
    else:
        atmospheric_loss_db = args.atmospheric_loss

    logging.info("Atmospheric loss:   {:6.2f} dB".format(atmospheric_loss_db))

    # -------- Rx dish gain --------
    if (args.rx_dish_gain is None):
        dish_gain_db = calc.dish_gain(args.rx_dish_size, args.freq,
                                      args.rx_dish_efficiency)
        logging.info("Rx dish gain:       {:6.2f} dB".format(dish_gain_db))
    else:
        dish_gain_db = args.rx_dish_gain

    # -------- Noise figure --------
    coax_loss_db, coax_noise_fig_db = calc.coax_loss_nf(args.coax_length)

    if (args.lnb_noise_fig is None):
        lnb_noise_fig = calc.noise_temp_to_noise_fig(args.lnb_noise_temp)
    else:
        lnb_noise_fig = args.lnb_noise_fig

    logging.info("LNB noise figure:   {:6.2f} dB".format(lnb_noise_fig))

    noise_fig_db = calc.total_noise_figure(
        [lnb_noise_fig, coax_noise_fig_db, args.rx_noise_fig],
        [args.lnb_gain, -coax_loss_db]
    )

    # -------- System noise temperature --------
    effective_input_noise_temp = calc.noise_fig_to_noise_temp(noise_fig_db)

    logging.info("Input-noise temp:   {:6.2f} K".format(
        effective_input_noise_temp))

    if (args.antenna_noise_temp is None):
        antenna_noise_temp = calc.antenna_noise_temp(args.atmospheric_loss)
        # NOTE: consider the atmospheric loss only once here even if analyzing
        # a radar system. This call determines the sky's contribution to the Rx
        # antenna noise temperature.
    else:
        antenna_noise_temp = args.antenna_noise_temp

    logging.info("Antenna noise temp: {:6.2f} K".format(antenna_noise_temp))

    T_syst = calc.rx_sys_noise_temp(antenna_noise_temp,
                                    effective_input_noise_temp)
    T_syst_db = util.abs_to_db(T_syst)  # in dBK (for T_syst in K)

    # -------- Received Power, Noise Power, CNR, and other metrics --------
    P_rx_dbw = calc.rx_power(eirp, path_loss_db, dish_gain_db,
                             atmospheric_loss_db, args.mispointing_loss)

    N_dbw = calc.noise_power(T_syst_db, args.if_bw,)

    g_over_t_db = calc.g_over_t(dish_gain_db, T_syst_db)

    cnr = calc.cnr(P_rx_dbw,  N_dbw)

    # -------- Capacity --------
    capacity = calc.capacity(cnr, args.if_bw)

    # Results
    res = {
        'pointing': {
            'elevation': elevation,
            'azimuth': azimuth,
            'slant_range': slant_range
        },
        'eirp_db': eirp,
        'path_loss_db': path_loss_db,
        'atmospheric_loss_db': atmospheric_loss_db,
        'rx_dish_gain_db': dish_gain_db,
        'noise_fig_db': {
            'lnb': lnb_noise_fig,
            'coax': coax_noise_fig_db,
            'total': noise_fig_db
        },
        'noise_temp_k': {
            'effective_input': effective_input_noise_temp,
            'system': T_syst
        },
        'power_dbw': {
            'carrier': P_rx_dbw,
            'noise':  N_dbw
        },
        'g_over_t_db': g_over_t_db,
        'cnr_db': cnr,
        'capacity_bps': capacity

    }

    if (verbose and args.json):
        print(json.dumps(res, indent=4, ensure_ascii=True))

    return res


def main():
    parser = get_parser()
    args = parser.parse_args()
    validate(parser, args)
    analyze(args, verbose=True)