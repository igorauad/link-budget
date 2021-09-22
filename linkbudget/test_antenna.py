"""
References:

 [1] Couch, Leon W.. Digital & Analog Communication Systems.
 [2] www.everythingrf.com/rf-calculators/parabolic-reflector-antenna-gain

"""

import unittest
from . import antenna
from .constants import SPEED_OF_LIGHT


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

    def test_dish_diameter_inference(self):
        """Test inference of the physical diameter"""
        freq = 12.45e9
        # Dish object with diameter informed by argument
        dish1 = antenna.Antenna(freq, diameter=0.45, efficiency=0.56)
        # Dish object informing the gain and the aperture efficiency only
        dish2 = antenna.Antenna(freq, gain=dish1.gain_db, efficiency=0.56)
        # Verify that the diameter can be inferred correctly
        self.assertEqual(dish1.diameter, dish2.diameter)
        # If the efficiency is not provided, the physical diameter cannot be
        # inferred. In this case, a warning is expected:
        with self.assertLogs(level='WARNING'):
            dish3 = antenna.Antenna(freq, gain=dish1.gain_db)
            self.assertIsNone(dish3.diameter)

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

    def test_antenna_off_axis_gain(self):
        """Test off-axis gain computation

        Test against the co-polar examples from ITU-R BO.1213-1.

        """

        # Co-polar example 1 (Figure 1)
        dish = antenna.Antenna(freq=11.7e9, diameter=0.6, efficiency=0.65)
        self.assertAlmostEqual(dish.gain_db, 35.5, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(0), dish.gain_db, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(2), 30, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(3), 23.1, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(4), 13.8, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(5), 11.5, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(10), 4, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(20), -3.5, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(30), -5, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(80), 0, places=1)

        # Co-polar example 2 (Figure 2)
        dish = antenna.Antenna(freq=12.2e9, diameter=0.45, efficiency=0.65)
        self.assertAlmostEqual(dish.gain_db, 33.3, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(0), dish.gain_db, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(2), 30, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(3), 25.8, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(4), 19.9, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(5), 12.4, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(10), 4, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(20), -3.5, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(30), -5, places=1)
        self.assertAlmostEqual(dish.off_axis_gain(80), 0, places=1)

        # The off-axis angle must be >=0 and < 180
        with self.assertRaises(ValueError):
            dish.off_axis_gain(-1)

        with self.assertRaises(ValueError):
            dish.off_axis_gain(180)

        with self.assertRaises(ValueError):
            dish.off_axis_gain(181)

        # The diameter/wavelength ratio must be >= 11
        diameter = 0.45
        wavelength = diameter * 10
        freq = SPEED_OF_LIGHT / wavelength
        dish = antenna.Antenna(freq=freq, diameter=diameter, efficiency=0.65)

        with self.assertRaises(ValueError):
            dish.off_axis_gain(0)
