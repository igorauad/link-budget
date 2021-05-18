from math import log10
import logging


def abs_to_db(val):
    return 10*log10(val)


def db_to_abs(val_db):
    return 10**(val_db/10)


def dbw_to_dbm(val_db):
    return val_db + 30


def format_rate(rate):
    """Format data rate given in bps"""

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
    logging.info("---------------------|--------------------")
    logging.info("Parameter            | Value")
    logging.info("---------------------|--------------------")


def log_result(parameter, value):
    """Log parameter-value result"""
    logging.info("{:20s} | {:20s}".format(parameter, value))
