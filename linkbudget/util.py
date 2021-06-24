from math import log10
import logging

from .constants import SPEED_OF_LIGHT


def lin_to_db(val):
    """Linear value to decibels

    Args:
        val : Linear scale value.

    Returns:
        Corresponding value in decibels.

    """
    return 10 * log10(val)


def db_to_lin(val_db):
    """Decibels to linear value

    Args:
        val : Value in decibels.

    Returns:
        Corresponding value in linear scale.

    """
    return 10**(val_db / 10)


def dbw_to_dbm(val_dbw):
    """Convert dBW to dBm

    Args:
        val_dbw : Value in dBW (decibels above 1 Watt).

    Returns:
        Value in dBm (decibels above 1 milliwatt).

    """
    return val_dbw + 30


def dbm_to_dbw(val_dbm):
    """Convert dBm to dBW

    Args:
        val_dbm : Value in dBm (decibels above 1 milliwatt).

    Returns:
        Value in dBW (decibels above 1 Watt).

    """
    return val_dbm - 30


def wavelength(freq):
    """Compute the radio wavelength for a given frequency

    Args:
        freq : Radio frequency in Hz.

    Returns:
        Wavelength in meters.

    """
    return SPEED_OF_LIGHT / freq


def _format_unit(value, thresholds, units):
    """Generate "value unit" string with an adequately scaled unit

    Select the highest unit that keeps the scaled value above unit. If none of
    the scaled values are above unit, output the lowest value (based on the
    lowest threshold).

    For example, consider the scaling of a bit rate of 1e4 bps with thresholds
    of 1e6 for Mbps, 1e3 for kbps, and 1 for bps. The implementation scans the
    thresholds in descending order, i.e., starting from 1e6, then 1e3, then
    1. In the given example, the value of 1e4 is equivalent to 0.01 Mbps, which
    is below unit (i.e., < 1). Hence, this first threshold is not suitable, and
    the implementation proceeds to the next. The next threshold is of 1e3
    (kbps) and results in 10 kbps. This value is above unit. Hence, the loop
    stops there and output 10 kbps.

    Args:
        value : Value to output as a string accompanied by the unit.
        thresholds : List of thresholds used to select the appropriate unit.
        units : List of strings with the units corresponding to each threshold.

    Returns:
        "value unit" string.

    """
    original_value = value
    zipped = zip(thresholds, units)
    for t, u in sorted(zipped, key=lambda t: t[0], reverse=True):
        value = original_value / t
        unit = u
        if (original_value >= t):
            break

    return "{:.2f} {}".format(value, unit)


def format_rate(rate):
    """Format data rate given in bps

    Args:
        rate : Data rate in bits per second (bps).

    Returns:
        String with the data rate and an adjusted unit such as Gbps, Mbps, or
        kbps.

    """
    return _format_unit(rate,
                        thresholds=[1e9, 1e6, 1e3, 1],
                        units=["Gbps", "Mbps", "kbps", "bps"])


def format_area(area):
    """Format area given in m2

    Args:
        area : Area in square meters (m2).

    Returns:
        String with the area and an adjusted unit such as cm2, dm2, or m2.

    """
    return _format_unit(area,
                        thresholds=[1, 1e-2, 1e-4],
                        units=["m2", "dm2", "cm2"])


def format_power(power):
    """Format power given in Watts

    Args:
        power : Power in Watts (W).

    Returns:
        String with the power value and an adjusted unit such as kW or mW.

    """
    return _format_unit(power,
                        thresholds=[1e3, 1, 1e-3],
                        units=["kW", "W", "mW"])


def log_header():
    """Log header to form a table with result logs"""
    w = 31  # width
    logging.info(w * "-" + "|" + w * "-")
    log_result("Parameter", "Value")
    logging.info(w * "-" + "|" + w * "-")


def log_result(parameter, value):
    """Log parameter-value result

    Args:
        parameter : Parameter to log.
        value     : Corresponding value.

    """
    logging.info("{:30s} | {:30s}".format(parameter, value))
