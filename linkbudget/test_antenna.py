"""
References:

 [1] Couch, Leon W.. Digital & Analog Communication Systems.
 [2] www.everythingrf.com/rf-calculators/parabolic-reflector-antenna-gain

"""

import unittest
from . import antenna


class TestAntenna(unittest.TestCase):
    def test_dish_gain(self):
        """Test dish gain calculation"""
        # Example 8-5 in [1]:
        dish = antenna.Antenna(freq=4e9, diameter=3.05, efficiency=0.557)
        self.assertAlmostEqual(
            dish.gain_db,
            39.6,  # expected gain in dB
            places=1)

        # Using [2] with an aperture efficiency of 56%:
        dish = antenna.Antenna(freq=12.45e9, diameter=0.45, efficiency=0.56)
        self.assertAlmostEqual(
            dish.gain_db,
            32.84568544,  # expected gain in dB
            places=1)

    def test_dish_aperture_inference(self):
        """Test inference of the effective aperture area"""
        freq = 4e9  # antenna's operating frequency

        # Create two antenna objects. The first can calculate the effective
        # aperture area directly based on the given physical area and the
        # aperture efficiency. The second is created directly with the gain
        # from dish 1, in which case the effective aperture area is inferred
        # only. The two effective aperture computations should match.
        dish1 = antenna.Antenna(freq=freq, diameter=3.05, efficiency=0.557)
        dish2 = antenna.Antenna(freq=freq, gain=dish1.gain_db)
        self.assertAlmostEqual(dish1.effective_aperture,
                               dish2.effective_aperture)

    def test_required_params(self):
        """Test the required constructor parameters"""

        # The operating frequency is always required:
        with self.assertRaises(TypeError):
            antenna.Antenna(diameter=3.05, efficiency=0.557)

        with self.assertRaises(ValueError):
            antenna.Antenna(freq=None, diameter=3.05, efficiency=0.557)

        with self.assertRaises(ValueError):
            antenna.Antenna(freq=None, gain=30)

        # When the gain is not provided, both the diameter and aperture
        # efficiency are required.
        with self.assertRaises(ValueError):
            antenna.Antenna(freq=4e9, efficiency=0.557)

        with self.assertRaises(ValueError):
            antenna.Antenna(freq=4e9, diameter=3.05)
