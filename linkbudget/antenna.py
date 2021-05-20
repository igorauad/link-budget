"""
Antenna Calculations

References:

 [1] Timothy Pratt, Jeremy E. Allnutt, "Satellite Communications", 3rd Ed.
 [2] Couch, Leon W.. Digital & Analog Communication Systems.

"""
from math import log10, pi
from . import util


class Antenna:
    """Parabolic antenna"""
    def __init__(self,
                 freq,
                 gain=None,
                 diameter=None,
                 efficiency=None,
                 label="Dish"):
        """Construct the antenna object

        Args:
            freq       : Operating frequency in Hz.
            gain       : Antenna gain in dB.
            diameter   : Diameter in m.
            efficiency : Aperture efficiency.
            label : Label used for logs.

        Note:
            According to [1], the aperture efficiency is typically in the range
            0.5–0.75 for parabolodial reflector antennas, lower for small
            antennas and higher for large Cassegrain and Gregorian
            antennas. For instance, Table 8-4 in [2] assumes an aperture
            efficiency of 0.56.

        """
        if (freq is None):
            raise ValueError("Operating frequency is required")

        able_to_calc_gain = diameter is not None and efficiency is not None
        if (gain is None and not able_to_calc_gain):
            raise ValueError("gain must be provided when the diameter, freq, "
                             "and efficiency parameters are not")

        self.diameter = diameter
        self.freq = freq
        self.aperture_efficiency = efficiency

        # When the gain is explicitly provided, infer the effective aperture
        # area from it. Otherwise, compute the effective aperture from the
        # given antenna diameter and aperture efficiency. Then, using the
        # computed effective aperture, compute the antenna gain.
        if (gain is None):
            self.effective_aperture = self._calc_eff_aperture(
                diameter, efficiency)
            self.gain_db = self._calc_gain(freq, self.effective_aperture)
        else:
            self.gain_db = gain
            self.effective_aperture = self._infer_eff_aperture(freq, gain)

        util.log_result("{} gain".format(label),
                        "{:.2f} dB".format(self.gain_db))
        util.log_result("{} eff. aperture".format(label),
                        "{:.2f} m2".format(self.effective_aperture))

    def _calc_eff_aperture(self, diameter, aperture_efficiency):
        """Compute the antenna's effective aperture area

        The effective aperture area is given by:

        Ae = 𝜂 * A,

        where A represents the antenna's physical aperture area.

        If the aperture is circular with a diameter D in meters (or radius r),
        the physical aperture area is given by:

        A = 𝜋 * r^2 = 𝜋 * D^2∕4 square meters.

        Returns:
            Effective aperture area in square meters (m^2).

        """
        radius = diameter / 2
        face_area = pi * (radius**2)  # assume circular area
        return face_area * aperture_efficiency

    def _calc_gain(self, freq, effective_aperture):
        """Calculate the parabolic dish gain

        The gain in linear units is given by:

        Gain = 4 * 𝜋 * Ae ∕ 𝜆**2,

        where Ae is the effective aperture and 𝜆 is the wavelength.

        Returns:
            Gain in dB.

        """
        wavelength = util.wavelength(freq)
        gain = effective_aperture * 4 * pi / (wavelength**2)
        return 10 * log10(gain)

    def _infer_eff_aperture(self, freq, gain_db):
        """Infer the effective aperture area from the antenna gain

        The effective aperture area can be inferred from the gain and
        wavelength, as follows:

        Ae = Gain * 𝜆**2 / 4 * 𝜋.

        Returns:
            Effective aperture area in square meters (m^2).

        """
        wavelength = util.wavelength(freq)
        gain = util.db_to_abs(gain_db)
        return gain * (wavelength**2) / (4 * pi)
