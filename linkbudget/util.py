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


def format_rate(rate):
    """Format data rate given in bps

    Args:
        rate : Data rate in bits per second (bps).

    Returns:
        String with the data rate and an adjusted unit such as Gbps, Mbps, or
        kbps.

    """

    thresholds = [1e9, 1e6, 1e3]
    units = ["Gbps", "Mbps", "kbps"]

    # Default to bps
    value = rate
    unit = "bps"

    # Scale to a higher unit if possible
    for t, u in zip(thresholds, units):
        if (rate > t):
            value = rate / t
            unit = u
            break

    return "{:.2f} {}".format(value, unit)


def log_header():
    """Log header to form a table with result logs"""
    w = 25  # width
    logging.info(w * "-" + "|" + w * "-")
    log_result("Parameter", "Value")
    logging.info(w * "-" + "|" + w * "-")


def log_result(parameter, value):
    """Log parameter-value result

    Args:
        parameter : Parameter to log.
        value     : Corresponding value.

    """
    logging.info("{:24s} | {:24s}".format(parameter, value))
