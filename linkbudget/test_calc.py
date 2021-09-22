"""
References:

 [1] Couch, Leon W.. Digital & Analog Communication Systems.
 [2] Lindgren, M. (2015). A 1296 MHz Earth–Moon–Earth Communication System
     (Master's thesis).
 [3] Timothy Pratt, Jeremy E. Allnutt, "Satellite Communications", 3rd Ed.
 [4] www.everythingrf.com/rf-calculators/parabolic-reflector-antenna-gain

"""

import unittest
from unittest.mock import Mock
from math import pi

from . import calc, util


class TestBudgetCalc(unittest.TestCase):
    def test_eirp(self):
        tx_power_dbw = 10
        tx_dish_gain_db = 17
        eirp_dbw = calc.eirp(tx_power_dbw, tx_dish_gain_db)
        self.assertEqual(eirp_dbw, 27)

    def test_carrier_eirp(self):
        # Example 6.2 from [3] for the case of Station A.
        tp_sat_power = 20
        carrier_power = calc.carrier_eirp(tp_sat_power,
                                          obo=3,
                                          peb=15e6,
                                          tp_bw=30e6)
        self.assertAlmostEqual(carrier_power, 14, places=1)
        # NOTE the carrier_eirp function takes the saturated EIRP, while the
        # above example considers the transponder output power only (before the
        # antenna). Nevertheless, the calculation is equivalent.

        # The transponder bandwidth must be provided for the PEB calculation:
        with self.assertRaises(ValueError):
            carrier_power = calc.carrier_eirp(tp_sat_power, obo=3, peb=15e6)

    def test_path_loss(self):
        # Example 4.2 from [3]
        self.assertAlmostEqual(
            calc.path_loss(d=40e6, freq=11e9),
            205.3,  # expected path loss
            places=1)

        # Examples based on Section 10.2 from [2]:
        distance = 364288e3
        freq = 1296e6
        rcs = 0.065 * 9.49e12

        # One-way path loss
        self.assertAlmostEqual(
            calc.path_loss(d=distance, freq=freq),
            206,  # expected one-way path loss
            delta=0.1)

        # Monostatic radar path loss
        radar_obj_gain_db = calc.radar_obj_gain(freq, rcs)
        self.assertAlmostEqual(
            calc.path_loss(d=distance,
                           freq=freq,
                           radar=True,
                           obj_gain=radar_obj_gain_db),
            270,  # expected two-way transmission loss
            delta=0.25)

        # Equivalent bistatic radar path loss (assuming the object-to-Rx
        # distance is equal to the Tx-to-object distance)
        self.assertAlmostEqual(
            calc.path_loss(d=distance,
                           freq=freq,
                           radar=True,
                           obj_gain=radar_obj_gain_db,
                           bistatic=True,
                           d_rx=distance),
            270,  # expected two-way transmission loss
            delta=0.25)

        # The gain of the radar object must be provided in radar mode
        with self.assertRaises(ValueError):
            calc.path_loss(d=distance, freq=freq, radar=True)

        # The Rx distance must be provided in bistatic radar mode
        with self.assertRaises(ValueError):
            calc.path_loss(d=distance,
                           freq=freq,
                           radar=True,
                           obj_gain=radar_obj_gain_db,
                           bistatic=True)

    def test_radar_object_gain(self):
        # Based on the example in Section 10.2 from [2]:
        freq = 1296e6
        rcs = 0.065 * 9.49e12
        Gobj = calc.radar_obj_gain(freq, rcs)
        self.assertAlmostEqual(Gobj, 141.7, delta=0.1)

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
        Ae_db = util.lin_to_db(Ae)  # effective aperture in dBm^2
        Pr_dbw = F_dbw_m2 + Ae_db
        self.assertAlmostEqual(Pr_dbw, -126, places=1)

        # This Rx power computation based on the flux density should match with
        # the Rx power computed using the antenna gain (i.e., based on the
        # "link equation"). Consider Example 4.2 from [3].
        freq = 11e9
        wavelength = util.wavelength(freq)
        rx_gain = util.lin_to_db(4 * pi * Ae / wavelength**2)
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
        Tsys_dbk = util.lin_to_db(Tsys)
        self.assertAlmostEqual(Tsys_dbk, 18.01, places=2)

    def test_spectral_density(self):
        # Double-sided power spectral density of the white noise produced by a
        # thermal noise source at 290 K:
        N0_over_2_dbm_hz = -174  # dBm/Hz
        N0_over_2_dbw_hz = util.dbm_to_dbw(N0_over_2_dbm_hz)  # dBW/Hz
        N0_over_2_w_hz = util.db_to_lin(N0_over_2_dbw_hz)  # Watts/Hz
        bw = 1e6  # Hz
        N = N0_over_2_w_hz * bw  # Watts
        N_dbw = util.lin_to_db(N)  # dBW
        self.assertEqual(N0_over_2_dbw_hz, calc.spectral_density(N_dbw, bw))

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

    def test_carrier_to_asi_ratio(self):
        # 45 cm dish example from ITU-R BO.1213-1, with an on-axis gain of 33.3
        # dB and an off-axis gain of roughly 30 dB at a 2° off-axis angle:
        rx_dish = Mock()
        rx_dish.gain_db = 33.3
        rx_dish.off_axis_gain.return_value = 30
        long_separation = 2
        # If the adjacent satellites transmit with the same EIRP as the wanted
        # satellite, the C/I becomes equal to the difference between the
        # on-axis and off-axis antenna gains, which is of 3.3 dB:
        c_to_i_db = calc.carrier_to_asi_ratio(rx_dish,
                                              long_separation,
                                              asi_eirp_ratio=1)
        self.assertAlmostEqual(c_to_i_db, 3.3, places=1)

        # If the adjacent satellites transmit with an aggregate EIRP twice as
        # stronger as the wanted EIRP, the C/I reduces by 3 dB:
        c_to_i_db = calc.carrier_to_asi_ratio(rx_dish,
                                              long_separation,
                                              asi_eirp_ratio=2.0)
        self.assertAlmostEqual(c_to_i_db, 0.3, places=1)

        # Lastly, in the opposite case, if the adjacent satellites collectively
        # transmit with half the wanted EIRP, the C/I increases by 3 dB:
        c_to_i_db = calc.carrier_to_asi_ratio(rx_dish,
                                              long_separation,
                                              asi_eirp_ratio=0.5)
        self.assertAlmostEqual(c_to_i_db, 6.3, places=1)

    def test_cnir(self):
        # Example 4.10.b from [3]:
        cnr_db = 17
        ci_db = 24
        self.assertAlmostEqual(calc.cnir(cnr_db, ci_db), 16.2, places=1)
