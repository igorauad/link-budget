"""Link budget analysis unit tests

References:
 [1] Couch, Leon W.. Digital & Analog Communication Systems.
 [2] Lindgren, M. (2015). A 1296 MHz Earth–Moon–Earth Communication System
     (Master's thesis).

"""
import unittest
from . import main


class TestBudgetAnalysis(unittest.TestCase):
    def test_ku_band_example(self):
        # Example SA8-1 from [1]:
        parser = main.get_parser()
        args = parser.parse_args(
            ['--eirp', '52',
             '--freq', '12.45e9',
             '--if-bw', '24e6',
             '--rx-dish-size', '0.46',
             '--antenna-noise-temp', '20',
             '--lnb-noise-fig', '0.6',
             '--lnb-gain', '40',
             '--coax-length', '110',
             '--rx-noise-fig', '10',
             '--sat-long', '-101',
             '--rx-long', '-82.43',
             '--rx-lat', '29.71']
        )
        res = main.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95, places=2)
        # The actual result in [1] is distinct due to various roundings.

    def test_radar_example(self):
        # Based on [2]. Most of the calculation is in Chapters 10 and 11. Other
        # important sections are highlighted below.
        parser = main.get_parser()
        args = parser.parse_args(
            ['--eirp', '55.6',
             '--freq', '1296e6',
             '--if-bw', '100',
             '--rx-dish-gain', '31.1',  # See Section 8.7.2
             '--antenna-noise-temp', '51.8',  # See Section 8.10
             '--lnb-noise-fig', '0.52',  # Chosen to achieve a total NF of 0.54
             # dB, as in Table 8.4
             '--lnb-gain', '36.2',  # LNA gain in Table 8.4
             '--coax-length', '32.8',  # 32.8 ft ~= 10 m
             '--rx-noise-fig', '10',  # See Table 8.4
             '--sat-long', '-172',  # Arbitrary moon longitude
             '--rx-long', '-46.6333',  # Arbitrary Rx position
             '--rx-lat', '-23.5505',  # Arbitrary Rx position
             '--radar',
             '--radar-alt', '355600e3',  # Earth surface to moon surface
             # distance
             '--radar-cross-section', '0.61685e12']
        )
        res = main.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 5.46, places=2)
        # In [2], the path loss is computed directly based on the distance
        # between the ground station and the moon, as given by the MoonSked
        # software. In contrast, the above example assumes arbitrary Rx and
        # moon positions. The moon is assumed to be above latitude 0 (like a
        # GEO satellite) and at an arbitrary longitude. The adopted values are
        # such that the transmission loss matches the 270.3 dB loss used in
        # [2]. Nevertheless, the final result still differs slightly from the
        # one presented in [2] due to rounding of intermediate results.

    def test_opts(self):
        """Test program options"""
        parser = main.get_parser()

        # Test EIRP determined indirectly
        base_args = [
            '--freq', '12.45e9',
            '--if-bw', '24e6',
            '--rx-dish-size', '0.46',
            '--antenna-noise-temp', '20',
            '--lnb-noise-fig', '0.6',
            '--lnb-gain', '40',
            '--coax-length', '110',
            '--rx-noise-fig', '10',
            '--sat-long', '-101',
            '--rx-long', '-82.43',
            '--rx-lat', '29.71'
        ]

        # EIRP based on the Tx power and the Tx dish gain
        args = parser.parse_args(base_args +
                                 ['--tx-power', '20',
                                  '--tx-dish-gain', '32'])
        main.validate(parser, args)
        res = main.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95, places=2)

        # EIRP based on the Tx power and the Tx dish size
        args = parser.parse_args(base_args +
                                 ['--tx-power', '5',
                                  '--tx-dish-size', '2.4'])  # 52.37 dBW EIRP
        main.validate(parser, args)
        res = main.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95 + 0.37, places=2)

        # Tx power without the Tx dish gain or size should throw error:
        with self.assertRaises(SystemExit):
            args = parser.parse_args(base_args + ['--tx-power', '20'])
            main.validate(parser, args)

        # Rx dish gain given directly instead of through the dish size
        base_args = [
            '--eirp', '52',
            '--freq', '12.45e9',
            '--if-bw', '24e6',
            '--antenna-noise-temp', '20',
            '--lnb-noise-fig', '0.6',
            '--lnb-gain', '40',
            '--coax-length', '110',
            '--rx-noise-fig', '10',
            '--sat-long', '-101',
            '--rx-long', '-82.43',
            '--rx-lat', '29.71'
        ]
        args = parser.parse_args(base_args + ['--rx-dish-size', '0.46'])
        main.validate(parser, args)
        res = main.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95, places=2)

        # LNB noise temperature instead of the LNB noise figure
        base_args = [
            '--eirp', '52',
            '--freq', '12.45e9',
            '--if-bw', '24e6',
            '--rx-dish-size', '0.46',
            '--antenna-noise-temp', '20',
            '--lnb-noise-temp', '42.964',
            '--lnb-gain', '40',
            '--coax-length', '110',
            '--rx-noise-fig', '10',
            '--sat-long', '-101',
            '--rx-long', '-82.43',
            '--rx-lat', '29.71'
        ]
        args = parser.parse_args(base_args)
        res = main.analyze(args)
        self.assertAlmostEqual(res['cnr_db'], 15.95, places=2)

        # Radar mode requires the radar object's altitude and cross section
        with self.assertRaises(SystemExit):
            args = parser.parse_args(base_args +
                                     ['--radar', '--radar-alt', '0'])
            main.validate(parser, args)

        with self.assertRaises(SystemExit):
            args = parser.parse_args(base_args +
                                     ['--radar', '--radar-cross-section', '0'])
            main.validate(parser, args)
