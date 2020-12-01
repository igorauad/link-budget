"""Main link budget options and analysis

To reproduce Example SA8-1 from [1]:

./link_budget.py --eirp 52 \
  --freq 12.45e9 \
  --if-bw 24e6 \
  --rx-dish-size 0.46 \
  --antenna-noise-temp 20 \
  --lnb-noise-fig 0.6 \
  --lnb-gain 40 \
  --coax-length 110 \
  --rx-noise-fig 10 \
  --sat-long -101 \
  --rx-long -82.43 \
  --rx-lat 29.71

References:

[1] Couch, Leon W.. Digital & Analog Communication Systems.

"""
import logging
import argparse
from . import calc, pointing, util


__version__ = "0.1.0"


def parser():
    parser = argparse.ArgumentParser(
        description="Link Budget Calculator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
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


def analyze(args):
    sat_alt = 35786e3 if not args.radar else args.radar_alt

    elevation, azimuth, slant_range = pointing.look_angles(
        args.sat_long, args.rx_long, args.rx_lat, sat_alt)

    # Compute the EIRP
    if (args.eirp is None):
        if args.tx_dish_gain is None:
            tx_gain = calc.dish_gain(args.tx_dish_size, args.freq)
        else:
            tx_gain = args.tx_dish_gain
        eirp = calc.eirp(args.tx_power, tx_gain)
        logging.info("Tx Power:           {:6.2f} kW".format(
            util.db_to_abs(args.tx_power)/1e3))
    else:
        eirp = args.eirp

    logging.info("EIRP:               {:6.2f} dBW ({:6.2f} kW)".format(
        eirp, util.db_to_abs(eirp)/1e3))

    path_loss_db = calc.path_loss(slant_range, args.freq, args.radar,
                                  args.radar_cross_section,
                                  args.radar_bistatic)
    # TODO support bistatic radar. Add distance from radar object to rx
    # station.

    if (args.rx_dish_gain is None):
        dish_gain_db = calc.dish_gain(args.rx_dish_size, args.freq)
    else:
        dish_gain_db = args.rx_dish_gain

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

    effective_input_noise_temp = calc.noise_fig_to_noise_temp(noise_fig_db)

    logging.info("Antenna noise temp: {:6.2f} K".format(
        args.antenna_noise_temp))
    logging.info("Input-noise temp:   {:6.2f} K".format(
        effective_input_noise_temp))

    T_syst_db = calc.rx_sys_noise_temp(args.antenna_noise_temp,
                                       effective_input_noise_temp)

    cnr = calc.cnr(eirp, path_loss_db, dish_gain_db, T_syst_db, args.if_bw)

    calc.capacity(cnr, args.if_bw)


def main():
    logging.basicConfig(level=logging.INFO)
    args = parser()
    analyze(args)
