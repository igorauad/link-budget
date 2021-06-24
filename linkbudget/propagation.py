import itur

from . import util


def atmospheric_attenuation(lat, lng, elevation, pol_skew, freq, availability,
                            dish_size, dish_eff):
    """Total attenuation due to multiple sources of atmospheric attenuation

    A wrapper to obtain the total atmospheric attenuation through ITU-Rpy.

    Args:
        lat (float): Earth station's latitude in degrees.
        lng (float): Earth station's longitude in degrees.
        elevation (float): Elevation angle in degrees.
        freq (float): Carrier frequency in Hz.
        availability (float): Target link availability within [95 to 99.999%].
        dish_size (float): Earth station's dish diameter in m.
        dish_eff (float): Dish aperture efficiency within [0, 1].

    Returns:
        float: Total atmospheric attenuation in dB.

    """
    f_ghz = freq / 1e9
    unavailability = 100 - availability
    a_g, a_c, a_r, a_s, a_t = itur.atmospheric_attenuation_slant_path(
        lat,
        lng,
        f_ghz,
        elevation,
        unavailability,
        dish_size,
        eta=dish_eff,
        tau=pol_skew,
        return_contributions=True)
    util.log_result("Rain attenuation @ {:g}%".format(unavailability),
                    "{:.2f}".format(a_r))
    util.log_result("Cloud attenuation @ {:g}%".format(unavailability),
                    "{:.2f}".format(a_c))
    util.log_result("Gaseous attenuation @ {:g}%".format(unavailability),
                    "{:.2f}".format(a_g))
    util.log_result("Tropo. scintillation @ {:g}%".format(unavailability),
                    "{:.2f}".format(a_s))
    return a_t.value
