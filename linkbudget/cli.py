"""Link budget command-line utility"""
import json
import logging
import argparse

from . import calc, constants, propagation, pointing, util
from .antenna import Antenna

__version__ = "0.1.5"


def get_parser():
    """Command-line arguments"""
    parser = argparse.ArgumentParser(
        description="Link Budget Calculator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    general_p = parser.add_argument_group('General Options')
    general_p.add_argument('--json',
                           action='store_true',
                           help='Print results in JSON format.')
    general_p.add_argument('-v',
                           '--version',
                           action='version',
                           version='%(prog)s {}'.format(__version__))

    margin_p = parser.add_argument_group('Margin Options')
    margin_p.add_argument(
        '--min-cnr',
        type=float,
        help='Target minimum carrier-to-noise ratio (CNR) in dB considering '
        'an ideal receiver (i.e., excluding the implementation margin).')
    margin_p.add_argument(
        '--impl-margin',
        type=float,
        default=0,
        help='Implementation margin in dB accounting for the non-ideal '
        'behavior of the receiver system.')

    freq_p = parser.add_argument_group('Frequency Options')
    freq_p.add_argument(
        '--freq',
        required=True,
        type=float,
        help='Downlink carrier frequency in Hz for satellite signals or '
        'simply the signal frequency in Hz for radar (passively reflected) '
        'signals.')
    freq_p.add_argument('--if-bw',
                        required=True,
                        type=float,
                        help='IF bandwidth in Hz.')

    pol_p = parser.add_argument_group('Polarization Options')
    pol_p.add_argument(
        '--polarization',
        choices=['circular', 'linear'],
        default='linear',
        help='Polarization of the transmitted electromagnetic wave.')

    tx_feed_p = parser.add_argument_group('Tx Feed Options')
    tx_pwr_group = tx_feed_p.add_mutually_exclusive_group(required=True)
    tx_pwr_group.add_argument(
        '--eirp',
        type=float,
        help="Carrier or transponder EIRP in dBW. If an output backoff is "
        "defined, this EIRP corresponds to the amplifier\'s saturated "
        "operation. Otherwise, it represents the actual operating EIRP.")
    tx_pwr_group.add_argument(
        '--tx-power',
        type=float,
        help='Power feeding the Tx antenna in dBW. If an output backoff is '
        'defined, this parameter represents the amplifier\'s saturated '
        'output power. Otherwise, it refers to the actual Tx power. ')
    tx_feed_p.add_argument('--obo',
                           type=float,
                           default=0,
                           help="Carrier or transponder output backoff in dB.")

    dish_p = parser.add_argument_group('Antenna Options')
    tx_dish_group = dish_p.add_mutually_exclusive_group()
    tx_dish_group.add_argument(
        '--tx-dish-size',
        type=float,
        help='Diameter in meters of the parabolic antenna used for '
        'transmission. Used when the power is specified through option '
        '--tx-power.')
    tx_dish_group.add_argument(
        '--tx-dish-gain',
        type=float,
        help='Gain in dBi of the parabolic antenna used for transmission. '
        'Used when the power is specified through option --tx-power.')
    dish_p.add_argument(
        '--tx-dish-efficiency',
        type=float,
        default=0.56,
        help='Aperture efficiency of the parabolic antenna used for '
        'transmission. Determines the antenna gain when the dish is specified '
        'by size (option ``--tx-dish-size``). Otherwise, when the gain is '
        'defined directly by option ``--tx-dish-gain``, this parameter is '
        'used to infer the diameter of an equivalent parabolic reflector.')
    rx_dish_group = dish_p.add_mutually_exclusive_group(required=True)
    rx_dish_group.add_argument('--rx-dish-size',
                               type=float,
                               help='Parabolic antenna (dish) diameter in m.')
    rx_dish_group.add_argument('--rx-dish-gain',
                               type=float,
                               help='Parabolic antenna (dish) gain in dBi.')
    dish_p.add_argument(
        '--rx-dish-efficiency',
        type=float,
        default=0.56,
        help='Aperture efficiency of the parabolic antenna used for '
        'reception. Determines the antenna gain when the dish is specified by '
        'size (option ``--rx-dish-size``). Otherwise, when the gain is '
        'defined directly by option ``--rx-dish-gain``, this parameter is '
        'used to infer the diameter of an equivalent parabolic reflector.')

    prop_p = parser.add_argument_group('Propagation Options')
    prop_p.add_argument(
        '--atmospheric-loss',
        type=float,
        help='Attenuation in dB experienced through the atmosphere. It should '
        'always include the clear air attenuation, and it could include other '
        'effects such as rain and cloud attenuation. When analyzing a radar '
        'system, note this option should determine the one-way attenuation, '
        'not the two-way. When omitted, the program assumes a reasonable '
        'atmospheric attenuation based on models from ITU-R recommendations.')
    prop_p.add_argument(
        '--availability',
        type=float,
        default=99,
        help='Target link availability in %% to consider on the atmospheric '
        'attenuation model. For instance, when targeting at a 99.9%% '
        'availability, the analysis is based on the atmospheric attenuation '
        'exceeded 0.1%% of the time.')

    interf_p = parser.add_argument_group('Interference Options')
    interf_p.add_argument(
        '--asi',
        action='store_true',
        default=False,
        help='Include adjacent satellite interference (ASI) in the link  '
        'budget considering neighbor satellites with overlapping coverage, '
        'frequency and polarization.')
    interf_p.add_argument(
        '--asi-eirp-ratio',
        type=float,
        default=1.0,
        help='Ratio between the aggregate downlink EIRP from adjacent '
        'satellites and the wanted signal\'s EIRP.')
    interf_p.add_argument(
        '--asi-long-separation',
        type=float,
        default=2.0,
        help='Longitudinal orbit separation in degrees between the wanted '
        'satellite and the adjacent satellite(s).')

    rx_p = parser.add_argument_group('Rx Options')
    rx_p.add_argument(
        '--antenna-noise-temp',
        type=float,
        help='Receive antenna\'s noise temperature in K. When omitted, '
        'this parameter is derived based on the atmospheric attenuation.')
    rx_p.add_argument('--mispointing-loss',
                      type=float,
                      default=0,
                      help='Loss in dB due to antenna mispointing.')
    lnb_noise_group = rx_p.add_mutually_exclusive_group(required=True)
    lnb_noise_group.add_argument('--lnb-noise-fig',
                                 type=float,
                                 help='LNB\'s noise figure in dB.')
    lnb_noise_group.add_argument('--lnb-noise-temp',
                                 type=float,
                                 help='LNB\'s noise temperature in K.')
    rx_p.add_argument('--lnb-gain',
                      required=True,
                      type=float,
                      help='LNB\'s gain.')
    rx_p.add_argument('--rx-noise-fig',
                      required=True,
                      type=float,
                      help='Receiver\'s noise figure in dB.')
    rx_p.add_argument(
        '--coax-length',
        required=True,
        type=float,
        help='Length of the coaxial transmission line between the LNB and the '
        'receiver in ft.')

    fdma_group = parser.add_argument_group(
        title='FDMA Carrier Power Options',
        description="Parameters to determine the power allocated to an FDMA "
        "carrier.")
    fdma_group.add_argument(
        '--carrier-peb',
        type=float,
        help="Power-equivalent bandwidth (PEB) in Hz assigned for the FDMA "
        "carrier. When provided, the EIRP computed from --eirp or --tx-power "
        "refers to the transponder, while the PEB determines the fraction of "
        "this transponder EIRP allocated to the carrier. In this case, note "
        "the output backoff must refer to the transponder, not the carrier. "
        "If the output backoff represents the carrier backoff, do not inform "
        "the carrier PEB.")
    fdma_group.add_argument(
        '--tp-bw',
        type=float,
        help="Transponder bandwidth in Hz, required if the PEB is provided.")

    pos_p = parser.add_argument_group(
        'Satellite and Earth Station Position Information')
    pos_p.add_argument(
        '--sat-long',
        type=float,
        help='Satellite\'s longitude. Negative to the West and positive to '
        'the East.')
    pos_p.add_argument(
        '--rx-long',
        type=float,
        help='Rx station\'s longitude. Negative to the West and positive '
        'to the East.')
    pos_p.add_argument(
        '--rx-lat',
        type=float,
        help='Rx station\'s latitude. Positive to the North and negative '
        'to the South.')
    pos_p.add_argument(
        '--slant-range',
        type=float,
        help='Slant path length in km between the Rx station and the '
        'satellite or reflector.')

    radar_p = parser.add_argument_group('Radar Options')
    radar_p.add_argument(
        '--radar',
        default=False,
        action='store_true',
        help='Activate radar mode so that the link budget considers the '
        'path loss to and back from object.')
    radar_p.add_argument('--radar-alt',
                         type=float,
                         help='Altitude of the radar object')
    radar_p.add_argument('--radar-cross-section',
                         type=float,
                         help='Radar cross-section of the radar object.')
    radar_p.add_argument(
        '--radar-bistatic',
        default=False,
        action='store_true',
        help='Bistatic radar scenario, i.e., radar transmitter and receiver '
        'are not collocated.')
    return parser


def validate(parser, args):
    """Validate command-line arguments"""

    # Validate the latitude and longitude coordinates
    lat = args.rx_lat
    lng1 = args.rx_long
    lng2 = args.sat_long
    if (lat is not None and (lat < -90 or lat > 90)):
        parser.error("--rx-lat must be within -90 (South) to +90 (North).")

    if (lng1 is not None and (lng1 < -180 or lng1 > 180)):
        parser.error("--rx-long must be within -180 (West) to +180 (East).")

    if (lng2 is not None and (lng2 < -180 or lng2 > 180)):
        parser.error("--sat-long must be within -180 (West) to +180 (East).")

    if (args.tx_power and args.tx_dish_size is None
            and args.tx_dish_gain is None):
        parser.error("Define either --tx-dish-size or --tx-dish-gain  "
                     "using option --tx-power")

    # Link availability must be within [95, 99.999] because that's the range
    # supported by the ITU-R models (see, e.g., Step 10 in Section 2.2.1.1 of
    # ITU-R P.618-13).
    if (args.availability > 99.999 or args.availability < 95):
        parser.error(
            "Target link availability must be between 95 and 99.999 %%")

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

    # The atmospheric loss models require the Rx/satellite coordinates.
    if len(missing_pos_args) > 0 and args.atmospheric_loss is None:
        logging.error(
            "Unable to model the atmospheric loss without the Rx and "
            "satellite coordinates")
        parser.error(
            "Define the --atmospheric-loss or substitute --slant-range by the "
            "Rx/satellite coordinates")

    # If the antenna noise temperature is given directly, it will take
    # precedence over the inferred temperature based on the atmospheric
    # loss. This approach makes more sense when the atmospheric loss is also
    # given directly. When the atmospheric loss is not given directly, it is
    # determined automatically from ITU-R models. In this case, it makes sense
    # to let the antenna noise temperature be determined by the same models
    # too. Warn the user to prevent unintenional usage of the option.
    if (args.atmospheric_loss is None and args.antenna_noise_temp is not None):
        logging.warning(
            "Atmospheric loss model used to determine the atmospheric "
            "attenuation but not the antenna noise temperature (overriden by "
            "option --antenna-noise-temp)")

    if (args.carrier_peb is not None and args.tp_bw is None):
        parser.error("Argument --tp-bw is required if option --carrier-peb "
                     "is defined")

    # Validate the implementation margin
    if (args.impl_margin < 0):
        parser.error("argument --impl-margin must be a non-negative number")


def configure_logging():
    """Configure the logging format"""
    logging_fmt = "%(levelname)s %(message)s"
    logging.basicConfig(level=logging.INFO, format=logging_fmt)
    # Keep the level name except if it's an INFO message
    logging.addLevelName(logging.INFO, '')


def analyze(args, verbose=False):
    """Main link budget analysis

    Args:
        args : Populated argparse namespace object.
        verbose : Verbose mode.

    Returns:
        Dictionary with the main link budget results.

    """
    if (verbose and not args.json):
        util.log_header()

    # -------- Look angles --------
    if (args.slant_range is None):
        # Satellite altitude
        sat_alt = constants.GEOSYNC_ORBIT if not args.radar else args.radar_alt

        # Look angles
        elevation, azimuth, slant_range_m = pointing.look_angles(
            args.sat_long, args.rx_long, args.rx_lat, sat_alt)

        # Polarization skew
        #
        # NOTE: this parameter is used in the atmospheric loss model. It is
        # assumed equal to 45 degrees for circular polarization, according to
        # Recommendation ITU-R P.838-3.
        if (args.polarization == 'linear'):
            pol_skew = pointing.polarization_angle(args.sat_long, args.rx_long,
                                                   args.rx_lat)
        else:
            pol_skew = 45
    else:
        elevation = azimuth = pol_skew = None
        slant_range_m = args.slant_range * 1e3  # km to m

    # -------- Tx dish gain --------
    # Skipped if the EIRP is given directly.
    if (args.eirp is None):
        if args.tx_dish_gain is None:
            tx_dish = Antenna(freq=args.freq,
                              diameter=args.tx_dish_size,
                              efficiency=args.tx_dish_efficiency,
                              label="Tx dish")
        else:
            tx_dish = Antenna(freq=args.freq,
                              gain=args.tx_dish_gain,
                              efficiency=args.tx_dish_efficiency,
                              label="Tx dish")

    # -------- Rx dish gain --------
    if (args.rx_dish_gain is None):
        rx_dish = Antenna(freq=args.freq,
                          diameter=args.rx_dish_size,
                          efficiency=args.rx_dish_efficiency,
                          label="Rx dish")
    else:
        rx_dish = Antenna(freq=args.freq,
                          gain=args.rx_dish_gain,
                          efficiency=args.rx_dish_efficiency,
                          label="Rx dish")

    # -------- EIRP --------
    if (args.eirp is None):
        eirp_dbw = calc.eirp(args.tx_power, tx_dish.gain_db)
        util.log_result(
            "Tx Power",
            "{}".format(util.format_power(util.db_to_lin(args.tx_power))))
    else:
        eirp_dbw = args.eirp

    # The EIRP computed above could refer to an entire
    # transponder. Furthermore, it could refer to the power level in a given
    # direction when the amplifier operates in saturation. The following
    # function converts the EIRP computed above to the actual carrier EIRP.
    eirp_dbw = calc.carrier_eirp(eirp_dbw,
                                 args.obo,
                                 peb=args.carrier_peb,
                                 tp_bw=args.tp_bw)

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
    #
    # When the atmospheric loss is defined on argument "--atmospheric-loss", it
    # is prioritized over the atmospheric loss models. In contrast, when it is
    # not defined (the default), as long as the Rx station's coordinates and
    # the elevation are known, the atmospheric loss is computed from ITU-R
    # models. If the Rx and satellite coordinates are not available, the
    # atmospheric loss must be provided, otherwise the parser throws an error.
    #
    # Note that, ultimately, the atmospheric loss has two effects. It defines
    # the antenna noise temperature on reception (i.e., adds noise), and it
    # increases the path loss (i.e., reduces the signal power). In other words,
    # it disturbs both the numerator and the denominator of the C/N.
    if (args.atmospheric_loss is not None):
        one_way_atmospheric_loss_db = args.atmospheric_loss
    elif (all(v is not None for v in [args.rx_lat, args.rx_lat, elevation])):
        # Model the atmospheric attenuation for a given link availability.
        one_way_atmospheric_loss_db = propagation.atmospheric_attenuation(
            args.rx_lat, args.rx_long, elevation, pol_skew, args.freq,
            args.availability, rx_dish.diameter, rx_dish.aperture_efficiency)

    # On radar systems, assume the atmospheric loss is experienced twice.
    if (args.radar):
        two_way_atmospheric_loss_db = 2 * one_way_atmospheric_loss_db
        util.log_result("Total one-way atmospheric loss",
                        "{:.2f} dB".format(one_way_atmospheric_loss_db))
        util.log_result("Total two-way atmospheric loss",
                        "{:.2f} dB".format(two_way_atmospheric_loss_db))
    else:
        util.log_result("Total atmospheric loss",
                        "{:.2f} dB".format(one_way_atmospheric_loss_db))

    total_atmospheric_loss_db = two_way_atmospheric_loss_db if (args.radar) \
        else one_way_atmospheric_loss_db

    # -------- Reflected EIRP (radar mode only) --------
    # Compute the equivalent EIRP reflected off the radar object
    if (args.radar):
        one_way_path_loss_db = calc._path_loss(slant_range_m, args.freq)
        reflected_eirp_dbw = eirp_dbw \
            - one_way_path_loss_db \
            - one_way_atmospheric_loss_db \
            + radar_obj_gain_db
        util.log_result("Reflected EIRP",
                        "{:.2f} dBW".format(reflected_eirp_dbw))

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

    util.log_result("Input-noise temperature",
                    "{:.2f} K".format(effective_input_noise_temp))

    if (args.antenna_noise_temp is None):
        antenna_noise_temp = calc.antenna_noise_temp(
            one_way_atmospheric_loss_db)
        # NOTE: consider the atmospheric loss only once here even if analyzing
        # a radar system. This call determines the sky's contribution to the Rx
        # antenna noise temperature.
    else:
        antenna_noise_temp = args.antenna_noise_temp

    util.log_result("Antenna noise temperature",
                    "{:.2f} K".format(antenna_noise_temp))

    T_syst = calc.rx_sys_noise_temp(antenna_noise_temp,
                                    effective_input_noise_temp)
    T_syst_db = util.lin_to_db(T_syst)  # in dBK (for T_syst in K)

    # -------- Received Power and Flux Density --------
    if (args.radar):
        rx_flux_dbw_m2 = calc.rx_flux_density(reflected_eirp_dbw,
                                              slant_range_m,
                                              one_way_atmospheric_loss_db)
        # TODO: use the radar-to-Rx distance instead of the Tx-to-radar slant
        # range when running in bistatic radar mode (when supported).
    else:
        assert (one_way_atmospheric_loss_db == total_atmospheric_loss_db)
        rx_flux_dbw_m2 = calc.rx_flux_density(eirp_dbw, slant_range_m,
                                              one_way_atmospheric_loss_db)
    P_rx_dbw = calc.rx_power(eirp_dbw, path_loss_db, rx_dish.gain_db,
                             total_atmospheric_loss_db, args.mispointing_loss)

    # -------- Intermediate frequency (IF) Power --------
    P_if_dbw = P_rx_dbw + args.lnb_gain
    util.log_result("IF Power (LNB Output)",
                    "{:.2f} dBm".format(util.dbw_to_dbm(P_if_dbw)))

    # -------- Noise Power --------
    N_dbw = calc.noise_power(
        T_syst_db,
        args.if_bw,
    )

    # -------- Signal and Noise Power Spectral Densities --------
    sig_psd_dbw_hz = calc.spectral_density(P_rx_dbw,
                                           args.if_bw,
                                           label="Rx signal")
    noise_psd_dbw_hz = calc.spectral_density(N_dbw, args.if_bw, label="Noise")

    # -------- G/T and CNR --------
    g_over_t_db = calc.g_over_t(rx_dish.gain_db, T_syst_db)
    cnr_db = calc.cnr(P_rx_dbw, N_dbw)

    # -------- ASI and C/N+I --------
    if (args.asi):
        ci_db = calc.carrier_to_asi_ratio(rx_dish, args.asi_long_separation,
                                          args.asi_eirp_ratio)
        cnir_db = calc.cnir(cnr_db, ci_db)
    else:
        cnir_db = cnr_db
        ci_db = None

    # -------- Capacity --------
    capacity = calc.capacity(cnir_db, args.if_bw)

    # Results
    res = {
        'pointing': {
            'elevation': elevation,
            'azimuth': azimuth,
            'polarization_skew': pol_skew,
            'slant_range': slant_range_m
        },
        'eirp_db': eirp_dbw,
        'path_loss_db': path_loss_db,
        'atmospheric_loss_db': total_atmospheric_loss_db,
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
            'noise': N_dbw,
            'if': P_if_dbw
        },
        'psd_dbw_hz': {
            'carrier': sig_psd_dbw_hz,
            'noise': noise_psd_dbw_hz
        },
        'rx_flux_dbw_m2': rx_flux_dbw_m2,
        'g_over_t_db': g_over_t_db,
        'cnr_db': cnr_db,
        'ci_db': ci_db,
        'cnir_db': cnir_db,
        'capacity_bps': capacity
    }

    # -------- Link Margin --------
    if (args.min_cnr is not None):
        effective_min_cnr = args.min_cnr + args.impl_margin
        margin_db = cnir_db - effective_min_cnr
        util.log_result("Link margin", "{:.2f} dB".format(margin_db))
        res['margin_db'] = margin_db

    if (verbose and args.json):
        print(json.dumps(res, indent=4, ensure_ascii=True))

    return res


def main():
    parser = get_parser()
    args = parser.parse_args()
    if (not args.json):
        configure_logging()
    validate(parser, args)
    analyze(args, verbose=True)
