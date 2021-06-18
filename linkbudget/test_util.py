import unittest

from . import util


class TestUtil(unittest.TestCase):
    def test_format_rate(self):
        self.assertEqual(util.format_rate(1e-1), "0.10 bps")
        self.assertEqual(util.format_rate(1e0), "1.00 bps")
        self.assertEqual(util.format_rate(1e1), "10.00 bps")
        self.assertEqual(util.format_rate(1e2), "100.00 bps")
        self.assertEqual(util.format_rate(1e3), "1.00 kbps")
        self.assertEqual(util.format_rate(1e4), "10.00 kbps")
        self.assertEqual(util.format_rate(1e5), "100.00 kbps")
        self.assertEqual(util.format_rate(1e6), "1.00 Mbps")
        self.assertEqual(util.format_rate(1e7), "10.00 Mbps")
        self.assertEqual(util.format_rate(1e8), "100.00 Mbps")
        self.assertEqual(util.format_rate(1e9), "1.00 Gbps")
        self.assertEqual(util.format_rate(1e10), "10.00 Gbps")

    def test_format_area(self):
        self.assertEqual(util.format_area(1e1), "10.00 m2")
        self.assertEqual(util.format_area(1e0), "1.00 m2")
        self.assertEqual(util.format_area(1e-1), "10.00 dm2")
        self.assertEqual(util.format_area(1e-2), "1.00 dm2")
        self.assertEqual(util.format_area(1e-3), "10.00 cm2")
        self.assertEqual(util.format_area(1e-4), "1.00 cm2")
        self.assertEqual(util.format_area(1e-5), "0.10 cm2")

    def test_format_power(self):
        self.assertEqual(util.format_power(1e4), "10.00 kW")
        self.assertEqual(util.format_power(1e3), "1.00 kW")
        self.assertEqual(util.format_power(1e2), "100.00 W")
        self.assertEqual(util.format_power(1e1), "10.00 W")
        self.assertEqual(util.format_power(1e0), "1.00 W")
        self.assertEqual(util.format_power(1e-1), "100.00 mW")
        self.assertEqual(util.format_power(1e-2), "10.00 mW")
        self.assertEqual(util.format_power(1e-3), "1.00 mW")
        self.assertEqual(util.format_power(1e-4), "0.10 mW")
