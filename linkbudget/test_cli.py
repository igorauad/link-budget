"""Link budget analysis unit tests

References:
 [1] Couch, Leon W.. Digital & Analog Communication Systems.
 [2] Lindgren, M. (2015). A 1296 MHz Earth–Moon–Earth Communication System
     (Master's thesis).
 [3] Timothy Pratt, Jeremy E. Allnutt, "Satellite Communications", 3rd Ed.

"""
import copy
import io
import json
import sys
import unittest
import unittest.mock
from itertools import combinations
from . import cli


class TestCli(unittest.TestCase):
    def setUp(self):
        # Arguments to configure Example SA8-1 from [1]:
        self.base_args = [
            '--eirp', '52', '--freq', '12.45e9', '--if-bw', '24e6',
            '--rx-dish-size', '0.46', '--rx-dish-efficiency', '0.557',
            '--atmospheric-loss', '0', '--antenna-noise-temp', '20',
            '--lnb-noise-fig', '0.6', '--lnb-gain', '40', '--coax-length',
            '110', '--rx-noise-fig', '10', '--sat-long', '-101', '--rx-long',
            '-82.43', '--rx-lat', '29.71'
        ]

    def test_ku_band_example(self):
        # Example SA8-1 from [1]:
        parser = cli.get_parser()
        args = parser.parse_args(self.base_args)
        cli.validate(parser, args)
        res = cli.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95, places=2)
        # The actual result in [1] is distinct due to various roundings.

    def test_c_band_rain_example(self):
        # Example 4.6 in [3] (clear air):
        parser = cli.get_parser()
        args = parser.parse_args([
            '--eirp',
            '34',  # 39 -2 dB backoff -3 dB beam edge contour
            '--freq',
            '4e9',
            '--if-bw',
            '30e6',
            '--rx-dish-gain',
            '49.7',
            '--mispointing-loss',
            '0.5',
            '--atmospheric-loss',
            '0.2',
            '--lnb-noise-temp',
            '45',
            '--lnb-gain',
            '60',  # not given - consider a high gain
            '--coax-length',
            '0',
            '--rx-noise-fig',
            '0',  # LNA noise only
            '--slant-range',
            '40e3'
        ])
        cli.validate(parser, args)
        res = cli.analyze(args, verbose=True)
        self.assertAlmostEqual(res['cnr_db'], 22.7, places=1)

        # Example 4.7 in [3] (same as Example 4.6, but with heavy rain):
        args = parser.parse_args([
            '--eirp',
            '34',  # 39 -2 dB backoff -3 dB beam edge contour
            '--freq',
            '4e9',
            '--if-bw',
            '30e6',
            '--rx-dish-gain',
            '49.7',
            '--mispointing-loss',
            '0.5',
            '--atmospheric-loss',
            '1.2',  # 1 dB rain attenuation
            '--lnb-noise-temp',
            '45',
            '--lnb-gain',
            '60',  # not given - consider a high gain
            '--coax-length',
            '0',
            '--rx-noise-fig',
            '0',  # LNA noise only
            '--slant-range',
            '40e3'
        ])
        cli.validate(parser, args)
        res = cli.analyze(args, verbose=True)
        # The book seems to have an error on the clear air noise power given in
        # Table 4.5b, which is given as −135.5 dBW instead of −136.2 dBW as
        # given in Table 4.5a. That is, it is considering a noise power 0.7 dBW
        # higher. In the end, the book considers a CNR of 18.2 dB in heavy
        # rain, but we assume the correct one is 18.9 dBW.
        self.assertAlmostEqual(res['cnr_db'], 18.9, places=1)

    def test_radar_example(self):
        # Based on [2]. Most of the calculation is in Chapters 10 and 11. Other
        # important sections are highlighted below.
        parser = cli.get_parser()
        radar_cross_section = 0.065 * 9.49e12  # See Section 10.2
        args = parser.parse_args([
            '--eirp',
            '55.6',
            '--freq',
            '1296e6',
            '--if-bw',
            '100',
            '--rx-dish-gain',
            '31.1',  # See Section 8.7.2
            '--atmospheric-loss',  # Neglected (See Section 3.5.6)
            '0',
            '--antenna-noise-temp',
            '51.8',  # See Section 8.10
            '--lnb-noise-fig',
            '0.52',  # Chosen to achieve a total NF of 0.54
            # dB, as in Table 8.4
            '--lnb-gain',
            '36.2',  # LNA gain in Table 8.4
            '--coax-length',
            '32.8',  # 32.8 ft ~= 10 m
            '--rx-noise-fig',
            '10',  # See Table 8.4
            '--slant-range',
            '364288',  # See Section 10.2
            '--radar',
            '--radar-cross-section',
            str(radar_cross_section)
        ])
        cli.validate(parser, args)
        res = cli.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 5.51, places=2)
        # The final result differs slightly from the one presented in [2] due
        # to various rounding of intermediate results.

    def test_ku_band_example_with_atmospheric_attn(self):
        """Test Ku band example while considering atmospheric attenuation"""
        # Same arguments based on Example SA8-1 from [1], but now with the
        # atmospheric attenuation omitted such that the tool computes it
        # automatically from models. The antenna noise temperature is also
        # omitted such that it can be derived from the atmospheric
        # attenuation. Assume also circular polarization to hit that condition
        # on code coverage.
        availability = 99.99
        args = [
            '--eirp', '52', '--freq', '12.45e9', '--if-bw', '24e6',
            '--polarization', 'circular', '--availability',
            str(availability), '--rx-dish-size', '0.46',
            '--rx-dish-efficiency', '0.557', '--lnb-noise-fig', '0.6',
            '--lnb-gain', '40', '--coax-length', '110', '--rx-noise-fig', '10',
            '--sat-long', '-101', '--rx-long', '-82.43', '--rx-lat', '29.71'
        ]
        parser = cli.get_parser()
        args = parser.parse_args(args)
        cli.validate(parser, args)
        res = cli.analyze(args)

        # With circular polarization, the LNB polarization skew should be
        # assumed equal to 45 degrees.
        self.assertEqual(res['pointing']['polarization_skew'], 45)

        # The atmospheric attenuation modeled for a high link availability
        # (99.99%) in Ku band is usually significant.
        expected_exceeded_attn_db = 8
        self.assertGreater(res['atmospheric_loss_db'],
                           expected_exceeded_attn_db)

        # The resulting C/N is expected to be lower than in
        # test_ku_band_example (which assumed zero atmospheric loss).
        self.assertLess(res['cnr_db'], 15.95 - expected_exceeded_attn_db)

        # TODO: replace this test by an example that can be referenced to a
        # publication. The current version is only guessing a reasonable
        # attenuation level so that the atmospheric modeling can be covered.

    def test_ku_band_example_with_margin(self):
        """Test Ku band example with the link margin for a target min CNR

        Introduce a hypothetical minimum CNR and implementation margin in the
        scenario of Example SA8-1 from [1].

        """
        min_cnr = 2
        impl_margin = 1
        parser = cli.get_parser()
        args = parser.parse_args(
            self.base_args +
            ['--min-cnr',
             str(min_cnr), '--impl-margin',
             str(impl_margin)])
        cli.validate(parser, args)
        res = cli.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95, places=2)
        self.assertAlmostEqual(res['margin_db'],
                               15.95 - (min_cnr + impl_margin),
                               places=2)

    def test_ku_band_example_with_asi(self):
        """Test Ku band scenario with adjacent satellite interference (ASI)

        According to ITU-R BO.1213-1, the difference between the on-axis and
        off-axis gain with a 45 cm dish of 0.65 aperture efficiency is roughly
        3.3 dB when operating at 12.2 GHz. Hence, the ASI (or C/I) will
        dominate the carrier to noise plus interference ratio (CNIR). In this
        case, the C/(N+I) has to be less than the C/I of 3.3 dB.

        """
        base_args = [
            '--eirp', '52', '--freq', '12.2e9', '--if-bw', '24e6',
            '--rx-dish-size', '0.45', '--rx-dish-efficiency', '0.65',
            '--atmospheric-loss', '0', '--antenna-noise-temp', '20',
            '--lnb-noise-fig', '0.6', '--lnb-gain', '40', '--coax-length',
            '110', '--rx-noise-fig', '10', '--sat-long', '-101', '--rx-long',
            '-82.43', '--rx-lat', '29.71', '--min-cnr', '1'
        ]
        parser = cli.get_parser()
        args = parser.parse_args(base_args)

        # Without ASI, the C/(N+I) (equal to the C/N) is more than 16.4 dB
        cli.validate(parser, args)
        res = cli.analyze(args)
        self.assertEqual(res['cnr_db'], res['cnir_db'])
        self.assertGreater(res['cnir_db'], 16.4)
        no_asi_cnir_db = res['cnir_db']
        no_asi_capacity_bps = res['capacity_bps']
        no_asi_margin_db = res['margin_db']

        # With ASI from an adjacent satellite 2° away, the C/(N+I) drops
        # drastically to a value close to the C/I
        parser = cli.get_parser()
        args = parser.parse_args(base_args +
                                 ['--asi', '--asi-long-separation',
                                  str(2)])
        cli.validate(parser, args)
        res = cli.analyze(args)
        self.assertAlmostEqual(res['ci_db'], 3.35, places=1)
        self.assertLess(res['cnir_db'], res['ci_db'])

        # The C/N should remain the same as the C/(N+I) from the evaluation
        # without ASI (where I=0):
        self.assertEqual(res['cnr_db'], no_asi_cnir_db)

        # The capacity and link margin should both consider the C/(N+I), not
        # the C/N, so they should have reduced:
        self.assertLess(res['capacity_bps'], no_asi_capacity_bps)
        self.assertLess(res['margin_db'], no_asi_margin_db)

    def test_opts(self):
        """Test program options"""
        parser = cli.get_parser()

        # Base parameters: same as the ones adopted in "test_ku_band_example"
        base_args = [
            '--freq', '12.45e9', '--if-bw', '24e6', '--rx-dish-size', '0.46',
            '--rx-dish-efficiency', '0.557', '--atmospheric-loss', '0',
            '--antenna-noise-temp', '20', '--lnb-noise-fig', '0.6',
            '--lnb-gain', '40', '--coax-length', '110', '--rx-noise-fig', '10',
            '--sat-long', '-101'
        ]

        # Test out-of-range latitude and longitude coordinates
        for rx_lat in [-91, 91]:
            with self.assertRaises(SystemExit):
                args = parser.parse_args(base_args + [
                    '--eirp', '52', '--sat-long', '-101', '--rx-long',
                    '-82.43', '--rx-lat',
                    str(rx_lat)
                ])
                cli.validate(parser, args)
        for rx_long in [-181, 181]:
            with self.assertRaises(SystemExit):
                args = parser.parse_args(base_args + [
                    '--eirp', '52', '--sat-long', '-101', '--rx-long',
                    str(rx_long), '--rx-lat', '29.71'
                ])
                cli.validate(parser, args)
        for sat_long in [-181, 181]:
            with self.assertRaises(SystemExit):
                args = parser.parse_args(base_args + [
                    '--eirp', '52', '--sat-long',
                    str(sat_long), '--rx-long', '-82.43', '--rx-lat', '29.71'
                ])
                cli.validate(parser, args)

        # Now with the valid coordinates adopted in "test_ku_band_example"
        base_args.extend(
            ['--sat-long', '-101', '--rx-long', '-82.43', '--rx-lat', '29.71'])

        # Set the EIRP of 52 dBW indirectly through the Tx power and dish gain
        args = parser.parse_args(base_args +
                                 ['--tx-power', '20', '--tx-dish-gain', '32'])
        cli.validate(parser, args)
        res = cli.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95, places=2)

        # Set the EIRP based on the Tx power and the Tx dish size. Consider a
        # Tx power of 5 dBW, and a 2.4m dish with 55.7% aperture efficiency,
        # which should lead to an EIRP of 52.37 dBW
        args = parser.parse_args(base_args + [
            '--tx-power', '5', '--tx-dish-size', '2.4', '--tx-dish-efficiency',
            '0.557'
        ])
        cli.validate(parser, args)
        res = cli.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95 + 0.37, places=2)

        # Tx power without the Tx dish gain or size should throw error:
        with self.assertRaises(SystemExit):
            args = parser.parse_args(base_args + ['--tx-power', '20'])
            cli.validate(parser, args)

        # Provide the carrier PEB but not the transponder BW
        with self.assertRaises(SystemExit):
            args = parser.parse_args(base_args +
                                     ['--eirp', '52', '--carrier-peb', '1e6'])
            cli.validate(parser, args)

        # Rx dish gain given directly instead of through the dish size
        base_args = [
            '--eirp', '52', '--freq', '12.45e9', '--if-bw', '24e6',
            '--atmospheric-loss', '0', '--antenna-noise-temp', '20',
            '--lnb-noise-fig', '0.6', '--lnb-gain', '40', '--coax-length',
            '110', '--rx-noise-fig', '10', '--sat-long', '-101', '--rx-long',
            '-82.43', '--rx-lat', '29.71'
        ]
        args = parser.parse_args(base_args + ['--rx-dish-gain', '33.024'])
        cli.validate(parser, args)
        res = cli.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95, places=2)

        # LNB noise temperature instead of the LNB noise figure
        base_args = [
            '--eirp', '52', '--freq', '12.45e9', '--if-bw', '24e6',
            '--rx-dish-size', '0.46', '--rx-dish-efficiency', '0.557',
            '--atmospheric-loss', '0', '--antenna-noise-temp', '20',
            '--lnb-noise-temp', '42.964', '--lnb-gain', '40', '--coax-length',
            '110', '--rx-noise-fig', '10', '--sat-long', '-101', '--rx-long',
            '-82.43', '--rx-lat', '29.71'
        ]
        args = parser.parse_args(base_args)
        res = cli.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95, places=2)

        # Radar mode requires the radar object's altitude and cross section
        with self.assertRaises(SystemExit):
            args = parser.parse_args(base_args +
                                     ['--radar', '--radar-alt', '0'])
            cli.validate(parser, args)

        with self.assertRaises(SystemExit):
            args = parser.parse_args(base_args +
                                     ['--radar', '--radar-cross-section', '0'])
            cli.validate(parser, args)

        # Slant range is mutual exclusive with satellite and Rx coordinates
        with self.assertRaises(SystemExit):
            args = parser.parse_args(base_args + ['--slant-range', 0])
            cli.validate(parser, args)

        # If the slant is not provided, the satellite and Rx station
        # coordinates are required
        base_args = [
            '--eirp', '52', '--freq', '12.45e9', '--if-bw', '24e6',
            '--rx-dish-size', '0.46', '--antenna-noise-temp', '20',
            '--lnb-noise-fig', '0.6', '--lnb-gain', '40', '--coax-length',
            '110', '--rx-noise-fig', '10'
        ]
        pos_args = [['--sat-long', '-101'], ['--rx-long', '-82.43'],
                    ['--rx-lat', '29.71']]
        # Test all combinations of one or two Rx/sat position args, all of
        # which should fail (all the three args are required)
        for n_comb in [1, 2]:
            for selected_pos_args in list(combinations(pos_args, n_comb)):
                with self.assertRaises(SystemExit):
                    test_args = copy.copy(base_args)
                    for pos_arg_pair in selected_pos_args:
                        test_args.extend(pos_arg_pair)
                    args = parser.parse_args(test_args)
                    cli.validate(parser, args)

        # When modeling the atmospheric loss, the satellite and earth station
        # coordinates are required.
        with self.assertLogs(level='ERROR'):  # An error message is logged.
            with self.assertRaises(SystemExit):
                args = parser.parse_args(base_args + ['--slant-range', 0])
                cli.validate(parser, args)

        # Also, the link availability must be within [95, 99.999].
        base_args += [
            '--sat-long', '-101', '--rx-long', '-82.43', '--rx-lat', '29.71'
        ]
        for availability in [50, 94.999, 99.9999, 100]:
            with self.assertRaises(SystemExit):
                args = parser.parse_args(
                    base_args +
                    ['--availability', str(availability)])
                cli.validate(parser, args)

        # A warning is printed if the --antenna-noise-temp option is given and
        # --atmospheric-loss is not.
        with self.assertLogs(level='WARNING'):
            args = parser.parse_args(base_args)
            cli.validate(parser, args)

        # Try a negative implementation margin
        with self.assertRaises(SystemExit):
            args = parser.parse_args(self.base_args + ['--impl-margin', '-1'])
            cli.validate(parser, args)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_json_verbose_output(self, mock_stdout):
        """Test the output when the analyzer runs in verbose and --json mode"""
        parser = cli.get_parser()
        args = parser.parse_args(self.base_args + ['--json'])
        # Run the analyzer in verbose mode, in which case it should print the
        # JSON results to stdout in addition to returning them.
        cli.validate(parser, args)
        json_res = cli.analyze(args, verbose=True)
        printed_json = json.loads(mock_stdout.getvalue())
        self.assertEqual(printed_json, json_res)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_main_entrypoint_json(self, mock_stdout):
        """Test the main program entrypoint with json output"""

        # Place arguments on argv
        old_sys_argv = sys.argv
        sys.argv = [old_sys_argv[0]] + self.base_args + ['--json']

        # Run the main entrypoint
        try:
            cli.main()
        finally:
            sys.argv = old_sys_argv

        # Check the JSON results printed to stdout
        printed_json = json.loads(mock_stdout.getvalue())
        self.assertAlmostEqual(printed_json['cnr_db'], 15.95, places=2)

    def test_main_entrypoint(self):
        """Test the main program entrypoint with normal (non-json) output"""

        # Place arguments on argv
        old_sys_argv = sys.argv
        sys.argv = [old_sys_argv[0]] + self.base_args

        # Run the main entrypoint and capture the info logs
        with self.assertLogs(level='INFO') as cm:
            try:
                cli.main()
            finally:
                sys.argv = old_sys_argv

        # Check the C/N logged to stdout
        for line in cm.output:
            if ('C/N' in line):
                cnr_db = float(line.split()[-2])
        self.assertAlmostEqual(cnr_db, 15.95, places=2)
