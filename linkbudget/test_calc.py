"""
References:

 [1] Couch, Leon W.. Digital & Analog Communication Systems.
 [2] Lindgren, M. (2015). A 1296 MHz Earth–Moon–Earth Communication System
     (Master's thesis).
 [3] Timothy Pratt, Jeremy E. Allnutt, "Satellite Communications", 3rd Ed.
 [4] www.everythingrf.com/rf-calculators/parabolic-reflector-antenna-gain

"""

import unittest
from math import pi

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

    def test_rx_flux_density(self):
        """Test flux density computation"""
        # Example 4.1 from [3]
        distance = 4e7  # distance in meters
        eirp_dbw = 27  # EIRP in dBW
        Ae = 10  # effective aperture in m2

        # Flux density
        F_dbw_m2 = calc.rx_flux_density(eirp_dbw, distance)
        self.assertAlmostEqual(F_dbw_m2, -136, places=1)

        # Corresponding Rx power
        Ae_db = util.abs_to_db(Ae)  # effective aperture in dBm^2
        Pr_dbw = F_dbw_m2 + Ae_db
        self.assertAlmostEqual(Pr_dbw, -126, places=1)

        # This Rx power computation based on the flux density should match with
        # the Rx power computed using the antenna gain (i.e., based on the
        # "link equation"). Consider Example 4.2 from [3].
        freq = 11e9
        wavelength = util.wavelength(freq)
        rx_gain = util.abs_to_db(4 * pi * Ae / wavelength**2)
        path_loss = calc.path_loss(d=distance, freq=freq)
        Pr_dbw_2 = calc.rx_power(eirp_dbw, path_loss, rx_gain)
        self.assertAlmostEqual(Pr_dbw, Pr_dbw_2)

        # Check the match between the two Rx power computations (based on flux
        # density and based on the link equation) when a non-zero atmospheric
        # loss is introduced in the analysis.
        atm_loss_db = 0.5
        F_dbw_m2 = calc.rx_flux_density(eirp_dbw,
                                        distance,
                                        atm_loss_db=atm_loss_db)
        Pr_dbw = F_dbw_m2 + Ae_db
        Pr_dbw_2 = calc.rx_power(eirp_dbw,
                                 path_loss,
                                 rx_gain,
                                 atm_loss_db=atm_loss_db)
        self.assertAlmostEqual(Pr_dbw, Pr_dbw_2)

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

        # If a single noise figure is provided, make sure the total noise
        # figure is equivalent to it.
        nf = 0.6
        total_nf = calc.total_noise_figure(nfs=[nf], gains=[])
        self.assertEqual(total_nf, nf)

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
