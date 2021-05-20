"""Link budget analysis"""
import json
import logging
import argparse
from . import calc, pointing, util, constants
from .antenna import Antenna

__version__ = "0.1.2"


def get_parser():
    """Command-line arguments"""
    parser = argparse.ArgumentParser(
        description="Link Budget Calculator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--json',
                        action='store_true',
                        help='Print results in JSON format.')
    tx_pwr_group = parser.add_mutually_exclusive_group(required=True)
    tx_pwr_group.add_argument('--eirp', type=float, help='EIRP in dBW.')
    tx_pwr_group.add_argument('--tx-power',
                              type=float,
                              help='Power feeding the Tx antenna in dBW.')
    tx_dish_group = parser.add_mutually_exclusive_group()
    tx_dish_group.add_argument(
        '--tx-dish-size',
        type=float,
        help='Diameter in meters of the parabolic antenna used for '
        'transmission. Used when the power is specified through option '
        '--tx-power')
    tx_dish_group.add_argument(
        '--tx-dish-gain',
        type=float,
        help='Gain in dBi of the parabolic antenna used for transmission. '
        'Used when the power is specified through option --tx-power')
    parser.add_argument(
        '--tx-dish-efficiency',
        type=float,
        default=0.56,
        help='Aperture efficiency of the parabolic antenna used for '
        'transmission. Considered when the dish is specified by size.')
    parser.add_argument(
        '--freq',
        required=True,
        type=float,
        help='Downlink carrier frequency in Hz for satellite signals or '
        'simply the signal frequency in Hz for radar (passively reflected) '
        'signals.')
    parser.add_argument('--if-bw',
                        required=True,
                        type=float,
                        help='IF bandwidth in Hz.')
    rx_dish_group = parser.add_mutually_exclusive_group(required=True)
    rx_dish_group.add_argument('--rx-dish-size',
                               type=float,
                               help='Parabolic antenna (dish) diameter in m.')
    rx_dish_group.add_argument('--rx-dish-gain',
                               type=float,
                               help='Parabolic antenna (dish) gain in dBi.')
    parser.add_argument(
        '--rx-dish-efficiency',
        type=float,
        default=0.56,
        help='Aperture efficiency of the parabolic antenna used for '
        'reception. Considered when the dish is specified by size.')
    sky_noise_group = parser.add_mutually_exclusive_group(required=True)
    sky_noise_group.add_argument(
        '--antenna-noise-temp',
        type=float,
        help='Receive antenna\'s noise temperature in K.')
    sky_noise_group.add_argument(
        '--atmospheric-loss',
        type=float,
        help='Attenuation in dB experienced through the atmosphere. It should '
        'always include the clear air attenuation, and it could also include '
        'other effects such as rain attenuation')
    lnb_noise_group = parser.add_mutually_exclusive_group(required=True)
    lnb_noise_group.add_argument('--lnb-noise-fig',
                                 type=float,
                                 help='LNB\'s noise figure in dB.')
    lnb_noise_group.add_argument('--lnb-noise-temp',
                                 type=float,
                                 help='LNB\'s noise temperature in K.')
    parser.add_argument('--lnb-gain',
                        required=True,
                        type=float,
                        help='LNB\'s gain.')
    parser.add_argument(
        '--coax-length',
        required=True,
        type=float,
        help='Length of the coaxial transmission line between the LNB and the '
        'receiver in ft.')
    parser.add_argument('--rx-noise-fig',
                        required=True,
                        type=float,
                        help='Receiver\'s noise figure in dB.')
    parser.add_argument('--mispointing-loss',
                        type=float,
                        default=0,
                        help='Loss in dB due to antenna mispointing')
    pos_p = parser.add_argument_group('sat/rx positioning options')
    pos_p.add_argument(
        '--sat-long',
        type=float,
        help='Satellite\'s longitude. Negative to the West and positive to '
        'the East')
    pos_p.add_argument(
        '--rx-long',
        type=float,
        help='Rx station\'s longitude. Negative to the West and positive '
        'to the East')
    pos_p.add_argument(
        '--rx-lat',
        type=float,
        help='Rx station\'s latitude. Positive to the North and negative '
        'to the South')
    pos_p.add_argument(
        '--slant-range',
        type=float,
        help='Slant path length in km between the Rx station and the '
        'satellite or reflector')
    radar_p = parser.add_argument_group('radar options')
    radar_p.add_argument(
        '--radar',
        default=False,
        action='store_true',
        help='Activate radar mode, so that the link budget considers the '
        'pathloss to and back from object')
    radar_p.add_argument('--radar-alt',
                         type=float,
                         help='Altitude of the radar object')
    radar_p.add_argument('--radar-cross-section',
                         type=float,
                         help='Radar cross section of the radar object')
    radar_p.add_argument(
        '--radar-bistatic',
        default=False,
        action='store_true',
        help='Bistatic radar scenario, i.e., radar transmitter and receiver '
        'are not collocated')
    return parser


def validate(parser, args):
    """Validate command-line arguments"""
    if (args.tx_power and args.tx_dish_size is None
            and args.tx_dish_gain is None):
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
            ", ".join(missing_pos_args)))
    elif (args.slant_range is not None and len(defined_pos_args) > 0):
        parser.error("argument{} {}: not allowed with argument "
                     "--slant-range".format(
                         "s" if len(defined_pos_args) > 1 else "",
                         ", ".join(defined_pos_args)))


def analyze(args, verbose=False):
    """Main link budget analysis

    Args:
        args : Populated argparse namespace object.
        verbose : Verbose mode.

    Returns:
        Dictionary with the main link budget results.

    """
    if (verbose and not args.json):
        logging_fmt = "%(message)s"
        logging.basicConfig(level=logging.INFO, format=logging_fmt)
        util.log_header()

    # -------- Look angles --------
    if (args.slant_range is None):
        # Satellite altitude
        sat_alt = constants.GEOSYNC_ORBIT if not args.radar else args.radar_alt

        # Look angles
        elevation, azimuth, slant_range_m = pointing.look_angles(
            args.sat_long, args.rx_long, args.rx_lat, sat_alt)
    else:
        elevation = azimuth = None
        slant_range_m = args.slant_range * 1e3  # km to m

    # -------- EIRP --------
    if (args.eirp is None):
        if args.tx_dish_gain is None:
            tx_dish = Antenna(freq=args.freq,
                              diameter=args.tx_dish_size,
                              efficiency=args.tx_dish_efficiency,
                              label="Tx dish")
        else:
            tx_dish = Antenna(freq=args.freq,
                              gain=args.tx_dish_gain,
                              label="Tx dish")

        eirp_dbw = calc.eirp(args.tx_power, tx_dish.gain_db)
        util.log_result(
            "Tx Power",
            "{:.2f} kW".format(util.db_to_abs(args.tx_power) / 1e3))
    else:
        eirp_dbw = args.eirp

    util.log_result("EIRP", "{:.2f} dBW".format(eirp_dbw))

    # -------- Path loss --------
    if (args.radar):
        radar_obj_gain_db = calc.radar_obj_gain(args.freq,
                                                args.radar_cross_section)
    else:
        radar_obj_gain_db = None
    path_loss_db = calc.path_loss(slant_range_m, args.freq, args.radar,
                                  radar_obj_gain_db, args.radar_bistatic)
    # TODO support bistatic radar. Add distance from radar object to rx
    # station.

    # -------- Atmospheric loss --------
    # When defined, the atmospheric loss defines the antenna noise temperature
    # on reception. It also adds to the free-space path loss. On radar systems,
    # assume the atmospheric loss is experienced twice.
    atmospheric_loss_db = args.atmospheric_loss or 0
    if (args.radar):
        util.log_result("One-way atmospheric loss",
                        "{:.2f} dB".format(atmospheric_loss_db))
        atmospheric_loss_db = 2 * atmospheric_loss_db
        util.log_result("Total atmospheric loss",
                        "{:.2f} dB".format(atmospheric_loss_db))
    else:
        util.log_result("Atmospheric loss",
                        "{:.2f} dB".format(atmospheric_loss_db))

    # -------- Reflected EIRP (radar mode only) --------
    # Compute the equivalent EIRP reflected off the radar object
    if (args.radar):
        one_way_path_loss_db = calc._path_loss(slant_range_m, args.freq)
        reflected_eirp_dbw = eirp_dbw \
            - one_way_path_loss_db \
            - (atmospheric_loss_db/2) \
            + radar_obj_gain_db
        util.log_result("Reflected EIRP",
                        "{:.2f} dBW".format(reflected_eirp_dbw))

    # -------- Rx dish gain --------
    if (args.rx_dish_gain is None):
        rx_dish = Antenna(freq=args.freq,
                          diameter=args.rx_dish_size,
                          efficiency=args.rx_dish_efficiency,
                          label="Rx dish")
    else:
        rx_dish = Antenna(freq=args.freq,
                          gain=args.rx_dish_gain,
                          label="Rx dish")

    # -------- Noise figure --------
    coax_loss_db, coax_noise_fig_db = calc.coax_loss_nf(args.coax_length)

    if (args.lnb_noise_fig is None):
        lnb_noise_fig = calc.noise_temp_to_noise_fig(args.lnb_noise_temp)
    else:
        lnb_noise_fig = args.lnb_noise_fig

    util.log_result("LNB noise figure", "{:.2f} dB".format(lnb_noise_fig))

    noise_fig_db = calc.total_noise_figure(
        [lnb_noise_fig, coax_noise_fig_db, args.rx_noise_fig],
        [args.lnb_gain, -coax_loss_db])

    # -------- System noise temperature --------
    effective_input_noise_temp = calc.noise_fig_to_noise_temp(noise_fig_db)

    util.log_result("Input-noise temp",
                    "{:.2f} K".format(effective_input_noise_temp))

    if (args.antenna_noise_temp is None):
        antenna_noise_temp = calc.antenna_noise_temp(args.atmospheric_loss)
        # NOTE: consider the atmospheric loss only once here even if analyzing
        # a radar system. This call determines the sky's contribution to the Rx
        # antenna noise temperature.
    else:
        antenna_noise_temp = args.antenna_noise_temp

    util.log_result("Antenna noise temp",
                    "{:.2f} K".format(antenna_noise_temp))

    T_syst = calc.rx_sys_noise_temp(antenna_noise_temp,
                                    effective_input_noise_temp)
    T_syst_db = util.abs_to_db(T_syst)  # in dBK (for T_syst in K)

    # -------- Received Power and Flux Density --------
    if (args.radar):
        rx_flux_dbw_m2 = calc.rx_flux_density(reflected_eirp_dbw,
                                              slant_range_m,
                                              (atmospheric_loss_db / 2))
        # TODO: use the radar-to-Rx distance instead of the Tx-to-radar slant
        # range when running in bistatic radar mode (when supported).
    else:
        rx_flux_dbw_m2 = calc.rx_flux_density(eirp_dbw, slant_range_m,
                                              atmospheric_loss_db)
    P_rx_dbw = calc.rx_power(eirp_dbw, path_loss_db, rx_dish.gain_db,
                             atmospheric_loss_db, args.mispointing_loss)

    # -------- Noise Power, G/T, and CNR --------
    N_dbw = calc.noise_power(
        T_syst_db,
        args.if_bw,
    )

    g_over_t_db = calc.g_over_t(rx_dish.gain_db, T_syst_db)

    cnr = calc.cnr(P_rx_dbw, N_dbw)

    # -------- Capacity --------
    capacity = calc.capacity(cnr, args.if_bw)

    # Results
    res = {
        'pointing': {
            'elevation': elevation,
            'azimuth': azimuth,
            'slant_range': slant_range_m
        },
        'eirp_db': eirp_dbw,
        'path_loss_db': path_loss_db,
        'atmospheric_loss_db': atmospheric_loss_db,
        'rx_dish_gain_db': rx_dish.gain_db,
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
            'noise': N_dbw
        },
        'rx_flux_dbw_m2': rx_flux_dbw_m2,
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
