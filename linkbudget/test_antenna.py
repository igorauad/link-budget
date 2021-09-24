"""
References:

 [1] Couch, Leon W.. Digital & Analog Communication Systems.
 [2] www.everythingrf.com/rf-calculators/parabolic-reflector-antenna-gain

"""

import unittest

import numpy as np

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

        Test a couple of arbitrary values, most of which can be computed
        directly as "32 - 25*log10(angle)".

        """
        # Antenna with diameter/wavelength = 23.4:
        dish1 = antenna.Antenna(freq=11.7e9, diameter=0.6, efficiency=0.65)
        self.assertAlmostEqual(dish1.gain_db, 35.5, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(0), dish1.gain_db, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(0.1), 35.4, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(1), 34.1, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(2), 30, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(3), 20.1, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(4), 16.9, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(5), 14.5, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(10), 7, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(20), -0.5, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(30), -4.92, places=1)
        self.assertAlmostEqual(dish1.off_axis_gain(80), -10, places=1)

        # Antenna with diameter/wavelength = 46.8 (phi_1 = 1.55):
        dish2 = antenna.Antenna(freq=11.7e9, diameter=1.2, efficiency=0.65)
        self.assertAlmostEqual(dish2.gain_db, 41.5, places=1)
        # The gains within the main lobe differ from the 0.6m antenna
        self.assertAlmostEqual(dish2.off_axis_gain(0), dish2.gain_db, places=1)
        self.assertAlmostEqual(dish2.off_axis_gain(0.1), 41.4, places=1)
        self.assertAlmostEqual(dish2.off_axis_gain(1), 36, places=1)
        self.assertAlmostEqual(dish2.off_axis_gain(2), 24.47, places=1)
        # In contrast, the sidelobe values follow the same expression:
        self.assertAlmostEqual(dish2.off_axis_gain(3), 20.1, places=1)
        self.assertAlmostEqual(dish2.off_axis_gain(4), 16.9, places=1)
        self.assertAlmostEqual(dish2.off_axis_gain(5), 14.5, places=1)
        self.assertAlmostEqual(dish2.off_axis_gain(10), 7, places=1)
        self.assertAlmostEqual(dish2.off_axis_gain(20), -0.5, places=1)
        self.assertAlmostEqual(dish2.off_axis_gain(30), -4.92, places=1)
        self.assertAlmostEqual(dish2.off_axis_gain(80), -10, places=1)

        # Ensure the antenna gain is monotonically decreasing
        for dish in [dish1, dish2]:
            last_gain = np.inf
            for angle in np.arange(0, 180, 1e-2):
                gain = dish.off_axis_gain(angle)
                self.assertGreaterEqual(last_gain, gain,
                                        "Failed for angle {}".format(angle))
                last_gain = gain

        # The off-axis angle must be >=0 and <= 180°:
        with self.assertRaises(ValueError):
            dish.off_axis_gain(-1)

        with self.assertRaises(ValueError):
            dish.off_axis_gain(180.01)

        # The model is valid for frequencies within [2, 31] GHz:
        freq = 1.9e9
        dish = antenna.Antenna(freq=freq, diameter=0.45, efficiency=0.65)
        with self.assertRaises(ValueError):
            dish.off_axis_gain(0)

        freq = 32.1e9
        dish = antenna.Antenna(freq=freq, diameter=0.45, efficiency=0.65)
        with self.assertRaises(ValueError):
            dish.off_axis_gain(0)
