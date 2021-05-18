"""
References:

 [1] Couch, Leon W.. Digital & Analog Communication Systems.
 [2] Lindgren, M. (2015). A 1296 MHz Earth–Moon–Earth Communication System
     (Master's thesis).
 [3] Timothy Pratt, Jeremy E. Allnutt, "Satellite Communications", 3rd Ed.
 [4] www.everythingrf.com/rf-calculators/parabolic-reflector-antenna-gain

"""

import unittest
from . import calc, util


class TestBudgetCalc(unittest.TestCase):
    def test_eirp(self):
        tx_power_dbw = 10
        tx_dish_gain_db = 17
        eirp_dbw = calc.eirp(tx_power_dbw, tx_dish_gain_db)
        self.assertEqual(eirp_dbw, 27)

    def test_path_loss(self):
        # Example 4.2 from [3]
        self.assertAlmostEqual(
            calc.path_loss(d=40e6, freq=11e9),
            205.3,  # expected path loss
            places=1)

        # Examples based on Section 10.2 from [2]:

        # One-way path loss
        self.assertAlmostEqual(
            calc.path_loss(d=364288e3, freq=1296e6),
            206,  # expected path loss
            delta=0.1)

        # Monostatic radar path loss
        rcs = 0.065 * 9.49e12
        self.assertAlmostEqual(
            calc.path_loss(d=364288e3, freq=1296e6, radar=True, rcs=rcs),
            270,  # expected transmission loss
            delta=0.25)

        # Equivalent bistatic radar path loss (assuming the object-to-Rx
        # distance is equal to the Tx-to-object distance)
        self.assertAlmostEqual(
            calc.path_loss(d=364288e3,
                           freq=1296e6,
                           radar=True,
                           rcs=rcs,
                           bistatic=True,
                           d_rx=364288e3),
            270,  # expected transmission loss
            delta=0.25)

        # The RCS must be provided in radar mode
        with self.assertRaises(ValueError):
            calc.path_loss(d=364288e3, freq=1296e6, radar=True)

        # The Rx distance must be provided in bistatic radar mode
        with self.assertRaises(ValueError):
            calc.path_loss(d=364288e3,
                           freq=1296e6,
                           radar=True,
                           rcs=rcs,
                           bistatic=True)

    def test_dish_gain(self):
        # Example 8-5 in [1]:
        self.assertAlmostEqual(
            calc.dish_gain(diameter=3.05, freq=4e9, efficiency=0.557),
            39.6,  # expected gain in dB
            places=1)

        # Using [4] with an aperture efficiency of 56%:
        self.assertAlmostEqual(
            calc.dish_gain(diameter=0.45, freq=12.45e9, efficiency=0.56),
            32.84568544,  # expected gain in dB
            places=1)

    def test_antenna_noise_temp(self):
        # Example 4.9 in [3]
        self.assertAlmostEqual(
            calc.antenna_noise_temp(attn_db=3.4),
            147,  # expected antenna noise temperature in K
            places=0)

    def test_coax_loss_nf(self):
        # Study aid example SA8-1 from [1]:
        loss_db, nf = calc.coax_loss_nf(length_ft=110, Tl=290)
        # At 290 K, the noise figure is equal to the line loss in dB:
        self.assertEqual(nf, 8.8)
        self.assertEqual(loss_db, 8.8)

    def test_total_nf(self):
        # Study aid example SA8-1 from [1]:
        nfs = [0.6, 8.8, 10]
        gains = [40, -8.8]
        total_nf = calc.total_noise_figure(nfs, gains)
        self.assertAlmostEqual(total_nf, 0.63, places=2)

    def test_noise_fig_temp_conv(self):
        # Table 4.4 from [3]:
        noise_temp = [0, 20, 40]
        noise_fig = [0, 0.29, 0.56]
        for i, nf in enumerate(noise_fig):
            # Noise figure to noise temperature:
            Te = calc.noise_fig_to_noise_temp(nf)
            self.assertAlmostEqual(Te, noise_temp[i], delta=0.1)
            # Noise temperature to noise figure:
            self.assertAlmostEqual(calc.noise_temp_to_noise_fig(Te), nf)

    def test_rx_sys_noise_temp(self):
        # Study aid example SA8-1 from [1]:
        Tsys = calc.rx_sys_noise_temp(Tar=20, Te=43.18)
        Tsys_dbk = util.abs_to_db(Tsys)
        self.assertAlmostEqual(Tsys_dbk, 18.01, places=2)

    def test_cnr(self):
        # Study aid example SA8-1 from [1]:
        C = calc.rx_power(eirp_db=52,
                          path_loss_db=205.73,
                          rx_ant_gain_db=32.96)
        N = calc.noise_power(T_sys_db=18.01, bw=24e6)
        cnr = calc.cnr(C, N)
        self.assertAlmostEqual(cnr, 16.03, places=1)

    def test_capacity(self):
        # At 0 dB SNR, the capacity in bits/complex-symbol becomes 1. The
        # corresponding capacity in bits per second (considering the symbol
        # rate) is equal to the passband (RF) bandwidth.
        self.assertEqual(
            calc.capacity(snr_db=0, bw=1e3),
            1e3  # expected capacity in bps
        )
