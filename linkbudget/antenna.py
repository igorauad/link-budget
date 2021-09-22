"""
Antenna Calculations

References:

- [1] Timothy Pratt, Jeremy E. Allnutt, "Satellite Communications", 3rd Ed.
- [2] Couch, Leon W.. Digital & Analog Communication Systems.
- [3] Recommendation ITU-R BO.1213-1.

"""
import logging
from math import log10, pi, sqrt

from . import util


class Antenna:
    """Parabolic antenna

    Args:
        freq (float): Operating frequency in Hz.
        gain (float): Antenna gain in dB.
        diameter (float): Diameter in m.
        efficiency (float): Aperture efficiency.
        label (str): Label used for logs.

    Note:
        According to [1], the aperture efficiency is typically in the range
        0.5–0.75 for parabolodial reflector antennas, lower for small antennas
        and higher for large Cassegrain and Gregorian antennas. For instance,
        Table 8-4 in [2] assumes an aperture efficiency of 0.56.

    Attributes:
        effective_aperture : Effective aperture area.
        gain_db : Antenna gain in dBi.

    """
    def __init__(self,
                 freq,
                 gain=None,
                 diameter=None,
                 efficiency=None,
                 label="Dish"):
        """Construct the antenna object"""
        if (freq is None):
            raise ValueError("Operating frequency is required")

        able_to_calc_gain = diameter is not None and efficiency is not None
        if (gain is None and not able_to_calc_gain):
            raise ValueError("gain must be provided when the diameter, freq, "
                             "and efficiency parameters are not")

        self.freq = freq
        self.diameter = diameter
        self.aperture_efficiency = efficiency

        # When the gain is explicitly provided, infer the effective aperture
        # area from it. Otherwise, compute the effective aperture from the
        # given antenna diameter and aperture efficiency. Then, using the
        # computed effective aperture, compute the antenna gain.
        if (gain is None):
            self.effective_aperture = self._calc_eff_aperture(
                diameter, efficiency)
            self.gain_db = self._calc_gain(freq, self.effective_aperture)
            util.log_result("{} aperture efficiency".format(label),
                            "{:.1f} %".format(efficiency * 100))
        else:
            self.gain_db = gain
            self.effective_aperture = self._infer_eff_aperture(freq, gain)
            # If the aperture efficiency is provided, infer the physical
            # diameter of an equivalent parabolic reflector.
            if (efficiency is not None):
                self.diameter = self._infer_diameter(self.effective_aperture,
                                                     self.aperture_efficiency)
                util.log_result("{} inferred diameter".format(label),
                                "{:.2f} m".format(self.diameter))
            else:
                logging.warning(
                    "Cannot infer the physical diameter of {} "
                    "(aperture efficiency required).".format(label))

        util.log_result("{} gain".format(label),
                        "{:.2f} dB".format(self.gain_db))
        util.log_result("{} effective aperture".format(label),
                        "{}".format(util.format_area(self.effective_aperture)))

    def _calc_eff_aperture(self, diameter, aperture_efficiency):
        """Compute the antenna's effective aperture area

        The effective aperture area is given by:

        .. math::

          A_e = \\eta A,

        where :math:`A` represents the antenna's physical aperture area and
        :math:`\\eta` is the aperture efficiency.

        If the aperture is circular with a diameter :math:`D` in meters (or
        radius :math:`r`), the physical aperture area is given by:

        .. math::

          A = 𝜋 r^2 = \\frac{𝜋 D^2}{4}.

        Hence, the effective aperture area becomes:

        .. math::

          A_e = \\frac{\\eta 𝜋 D^2}{4}.

        Args:
            diameter (float): Diameter in m.
            aperture_efficiency (float): Aperture efficiency.

        Returns:
            Effective aperture area in square meters (m^2).

        """
        radius = diameter / 2
        face_area = pi * (radius**2)  # assume circular area
        return face_area * aperture_efficiency

    def _calc_gain(self, freq, effective_aperture):
        """Calculate the parabolic dish gain

        The gain in linear units is given by:

        .. math::

          G = \\frac{4𝜋A_e}{𝜆^2},

        where :math:`A_e` is the effective aperture and 𝜆 is the wavelength.

        Args:
            freq (float): Operating frequency in Hz.
            effective_aperture (float): Effective aperture in square meters.

        Returns:
            Gain in dB.

        """
        wavelength = util.wavelength(freq)
        gain = effective_aperture * 4 * pi / (wavelength**2)
        return util.lin_to_db(gain)

    def _infer_eff_aperture(self, freq, gain_db):
        """Infer the effective aperture area from the antenna gain

        The effective aperture area can be inferred from the gain and
        wavelength, as follows:

        .. math::

          A_e = \\frac{G 𝜆^2}{4𝜋}.

        Args:
            freq (float): Operating frequency in Hz.
            gain_db (float): Antenna gain in dB.

        Returns:
            float: Effective aperture area in square meters (m^2).

        """
        wavelength = util.wavelength(freq)
        gain = util.db_to_lin(gain_db)
        return gain * (wavelength**2) / (4 * pi)

    def _infer_diameter(self, effective_aperture, aperture_efficiency):
        """Infer the diameter of an equivalent parabolic reflector

        The physical aperture area :math:`A` can be expressed in terms of the
        effective aperture :math:`A_e` and the aperture efficiency
        :math:`\\eta`, as follows:

        .. math::

          A = \\frac{A_e}{\\eta}.

        If the aperture is circular with a diameter D in meters, the physical
        aperture area is given by:

        .. math::

          A = \\frac{𝜋 D^2}{4}.

        Hence, it follows that:

        .. math::

          D = \\sqrt{\\frac{4 A_e}{\\eta \\pi}}.

        Note:
            This inference is useful when working with a non-parabolic antenna,
            such as a flat-panel antenna. In this case, you may know the
            antenna gain and aperture efficiency specifications, but not the
            physical diameter and aperture of the antenna. Meanwhile, the
            diameter may still be required for computations such as the
            tropospheric scintillation model from Recommendation ITU-R 618
            (see, e.g., `ITU-Rpy
            <https://itu-rpy.readthedocs.io/en/latest/apidoc/itu618.html>`_).

        Args:
            effective_aperture (float): Effective aperture in square meters.
            aperture_efficiency (float): Aperture efficiency.

        Returns:
            float: Diameter in meters (m) of an equivalent parabolic reflector
            with the same aperture efficiency.

        """
        return sqrt((4 * effective_aperture) / (aperture_efficiency * pi))

    def off_axis_gain(self, angle):
        """Compute the off-axis co-polar antenna gain

        Based on the co-polar antenna pattern formulae presented in Annex 1 of
        Recommendation ITU-R BO.1213-1.

        Args:
            angle (float): Off-axis angle in degrees relative to the boresight.
                Must be within the [0, 180°) range.

        Returns:
            (float) Off-axis antenna gain in dBi.

        """
        D_over_lambda = self.diameter / util.wavelength(self.freq)

        if (D_over_lambda < 11):
            raise ValueError(
                "Off-axis model from ITU-R BO.1213-1 requires D/lambda >= 11")

        Gmax = self.gain_db
        phi_r = 95 / D_over_lambda
        G1 = 29 - 25 * log10(phi_r)
        phi_m = sqrt((Gmax - G1) / 0.0025) / D_over_lambda
        phi_b = 10**(34 / 25)

        if (angle >= 0 and angle < phi_m):
            Gco = Gmax - 2.5e-3 * (D_over_lambda * angle)**2
        elif (angle >= phi_m and angle < phi_r):
            Gco = G1
        elif (angle >= phi_r and angle < phi_b):
            Gco = 29 - 25 * log10(angle)
        elif (angle >= phi_b and angle < 70):
            Gco = -5
        elif (angle >= 70 and angle < 180):
            Gco = 0
        else:
            raise ValueError(
                "Off-axis angle must be within the [0, 180°) range.")

        return Gco
