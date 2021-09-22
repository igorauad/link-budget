"""
A collection of link budget and other RF calculations.

References:

- [1] Couch, Leon W.. Digital & Analog Communication Systems.
- [2] Lindgren, M. (2015). A 1296 MHz Earth–Moon–Earth Communication System.
- [3] Timothy Pratt, Jeremy E. Allnutt, "Satellite Communications", 3rd Ed.

"""
from math import log10, pi, log2

from . import util
from .constants import T0


def eirp(tx_power, tx_dish_gain):
    """Compute the effective isotropically radiated power (EIRP)

    EIRP (dB) = Tx Power (dB) + Tx Antenna Gain (dB).

    Args:
        tx_power     : Transmit power feeding the antenna (dBW).
        tx_dish_gain : Transmit antenna gain (dB).

    Returns:
        The EIRP in dBW.

    """

    eirp = tx_power + tx_dish_gain
    return eirp


def carrier_eirp(sat_eirp, obo, peb=None, tp_bw=None):
    """Compute the carrier EIRP using transponder/amplifier information

    Converts the saturated EIRP obtained by an amplifier or transponder into
    the corresponding carrier EIRP, after output backoff and power allocation.

    There are two main use cases for this function:

    1. To compute the carrier EIRP from the transponder's saturated EIRP and
       carrier output backoff (OBO).

    2. To compute the carrier EIRP based on the transponder's saturated EIRP,
       the transponder's OBO, the carrier's power-equivalent bandwidth (PEB),
       and the transponder's bandwidth.

    In scenario 1, we start with the transponder's EIRP in dBW obtained when
    its amplifier operates in saturation. This metric is often given as the
    transponder's downlink EIRP at beam peak (BP) or beam center (BC). Then,
    assuming the transponder is shared by multiple FDMA carriers, we subtract
    the carrier OBO to obtain the EIRP assigned to the carrier of interest. In
    this case, note the **carrier OBO** is often a much larger value than the
    **transponder OBO**. The carrier OBO is easily in excess of 10 dB,
    depending on how much of the transponder it occupies. In contrast, as
    discussed in [3], the transponder OBO is typically within the range from 1
    to 7 dB in FDMA systems.

    The OBO interpretation changes in scenario 2 when considering FDMA
    systems. In this case, the OBO is interpreted as the **transponder OBO**,
    not the **carrier OBO**. Hence, the first step in the computation is to
    obtain the transponder's operating EIRP, which is given by the saturated
    EIRP minus the transponder OBO. The second step is to allocate a fraction
    of the available transponder EIRP to the carrier of interest. This fraction
    is determined by the ratio between the PEB and the transponder bandwidth.

    To compute scenario 1, call this function with parameter `obo` as the
    carrier OBO, while leaving the PEB undefined. To compute scenario 2, call
    this function with `obo` as the transponder OBO, while informing the PEB
    and transponder bandwidth. In any case, however, the `sat_eirp` parameter
    represents the transponder's saturated EIRP in the direction of interest.

    Note that, for convenience, the saturated EIRP is **in the direction of
    interest**, not necessarily at BP or BC. For example, it could be at beam
    edge (BE) or any arbitrary coverage contour. Furthermore, not this function
    takes the saturated EIRP (i.e., including the on-axis antenna gain),
    instead of the amplifier's saturated output power. This is intentional
    because it is often easier to find the EIRP information for a given
    satellite contour than to find the transponder and antenna gain separately.

    Lastly, note that both use cases are equivalent in TDMA systems. Since
    there is a single carrier in TDMA, the transponder OBO is the same as the
    carrier OBO. Furthermore, the PEB is the full transponder bandwidth in
    TDMA, so the `peb` parameter does not need to be informed.

    Args:
        sat_eirp  : Saturated amplifier/transponder EIRP in dBW in the
                    direction if interest.
        obo       : Output backoff in dB.
        peb       : Power equivalent bandwidth if the carrier shares the
                    amplifier/transponder with other FDMA carriers.
        tp_bw     : Transponder bandwidth used when the PEB is given.

    Returns:
        Carrier EIRP in dBW.

    """
    if (peb is None):
        carrier_eirp = sat_eirp - obo
    else:
        if (tp_bw is None):
            raise ValueError("The transponder bandwidth is required when "
                             "computing the carrier EIRP with the PEB")
        tp_eirp = sat_eirp - obo
        carrier_eirp = tp_eirp + util.lin_to_db(peb / tp_bw)

    return carrier_eirp


def _path_loss(d, freq):
    """One-way free-space path loss

    Args:
        d    : Distance in meters.
        freq : Carrier frequency in Hz.

    Returns:
        Path loss in dB.

    """
    wavelength = util.wavelength(freq)
    return 20 * log10(4 * pi * d / wavelength)


def path_loss(d, freq, radar=False, obj_gain=None, bistatic=False, d_rx=None):
    """Calculate the free-space path loss (or transmission loss)

    This function supports radar mode, in which case it computes the
    transmission loss considering the path loss in the forward and reverse
    paths (to and from the radar object) as well as the scattering at the radar
    object based on its radar cross section.

    Args:
        d        : Distance in meters between transmitter and receiver (e.g.,
                   satellite to ground station) or between the transmitter and
                   the radar object.
        freq     : Carrier frequency in Hz.
        radar    : Radar mode.
        obj_gain : Gain of the radar object in dB.
        bistatic : Bistatic radar mode.
        d_rx     : Bistatic radar mode only: distance between the radar object
                   and the Rx that is not collocated with the Tx.

    Returns:
        Path loss in dB.

    """
    # Eq. 8-11 from [1], or Eq. 3.16 from [2]:
    Lfs_one_way_db = _path_loss(d, freq)

    if (radar):
        if (obj_gain is None):
            raise ValueError("Radar object's gain is required in radar mode")

        if (bistatic):
            if (d_rx is None):
                raise ValueError("Rx distance required in bistatic radar mode")

            Lfs_tx_db = Lfs_one_way_db
            Lfs_rx_db = _path_loss(d_rx, freq)
            util.log_result("Uplink path loss", "{:.2f} dB".format(Lfs_tx_db))
            util.log_result("Downlink path loss",
                            "{:.2f} dB".format(Lfs_rx_db))
            # Bistatic radar transmission loss in dB, equation 3.24 in [2]:
            Lfs_db = Lfs_tx_db + Lfs_rx_db - obj_gain
        else:
            util.log_result("One-way path loss",
                            "{:.2f} dB".format(Lfs_one_way_db))
            # Monostatic radar transmission loss in dB, equation 3.26 in [2]:
            Lfs_db = 2 * Lfs_one_way_db - obj_gain

        util.log_result("Total path loss", "{:.2f} dB".format(Lfs_db))
    else:
        Lfs_db = Lfs_one_way_db
        util.log_result("Path loss", "{:.2f} dB".format(Lfs_db))
    return Lfs_db


def radar_obj_gain(freq, rcs):
    """Compute the gain of the radar object

    Args:
        freq : Carrier frequency in Hz.
        rcs  : Radar cross section (RCS).

    Returns:
        Gain of the radar object in dB.

    Note:
        The RCS definition repeated in [2] is the following: "the RCS of a
        radar object is the hypothetical area intercepting that amount of power
        which, when scattered isotropically, produces a power density at the
        receiver equal to that from the actual object."

    """
    wavelength = util.wavelength(freq)

    # Radar object gain in dB from equation 3.23 in [2]:
    G_obj_db = 10 * log10(4 * pi * rcs / (wavelength**2))
    util.log_result("Radar object's gain", "{:.2f} dB".format(G_obj_db))
    return G_obj_db


def antenna_noise_temp(attn_db, T_medium=270, coupling_eff=1.0):
    """Compute the antenna noise temperature based on the atmosphere attenuation

    Follow the theory in Section 4.5.2 of [3].

    Args:

        attn_db      : Total path attenuation (in dB) experienced through the
                       atmosphere, including clear air and rain attenuation.
        T_medium     : Medium temperature (in K) assumed for the rain.
        coupling_eff : Sky noise coupling coefficient determining the fraction
                       of incident sky noise energy output by the antenna.

    Note:
        In [3], a medium temperature of 270 K is considered when computing the
        sky noise temperature for rain. In contrast, a temperature of 290 K is
        considered when analyzing clear sky. The rationale is to be confirmed.

    """

    # Sky noise temperature (Eq. 4.30)
    T_sky = T_medium * (1 - (10**(-attn_db / 10)))

    # Antenna noise temperature (Eq. 4.31)
    T_a = coupling_eff * T_sky

    return T_a


def coax_loss_nf(length_ft, Tl=T0):
    """Compute the loss and noise figure of a coaxial RG6 transmission line

    Args:
        length_ft : Line length in feet.
        Tl        : temperature of the line in Kelvin.

    Returns:
        Tuple with line loss (dB) and noise figure (dB).

    """
    loss_db_per_ft = 8 / 100
    loss_db = length_ft * loss_db_per_ft
    loss = util.db_to_lin(loss_db)

    # The noise figure (dB) of a coaxial line is equal to the loss in dB if the
    # physical temperature of the line is equal to T0=290 K. See Equation 8.32a
    # on Example 8-2 in [1]. More generally, any passive two-port element (or
    # attenuator) at room temperature will have this property (noise figure =
    # attenuation in dB), see Equation 4.22 in [2].
    noise_factor = 1 + (Tl / T0) * (loss - 1)
    noise_fig = 10 * log10(noise_factor)

    util.log_result("Coax loss", "{:.2f} dB".format(loss_db))
    util.log_result("Coax noise figure", "{:.2f} dB".format(noise_fig))

    return loss_db, noise_fig


def total_noise_figure(nfs, gains):
    """Calculate the overall noise figure of the receiver system

    Args:
        nfs   : List with the noise figures (in dB) corresponding to the
                cascaded linear devices.
        gains : List with the gains (also in dB) of the cascaded linear
                devices, in the same order as given for "nfs".

    Note: The list of gains should not include the gain of the last device in
    the chain, as it is irrelevant for the overall noise figure computation.

    Returns:
        The overall noise figure in dB.

    """

    assert (len(nfs) > 0)
    assert (len(gains) == len(nfs) - 1)

    if (len(nfs) == 1):
        return nfs[0]

    # Implement Equation 8-34 from [1]:
    F = util.db_to_lin(nfs[0])
    G_prod = 1
    for i, nf in enumerate(nfs[1:]):
        nf_abs = util.db_to_lin(nf)
        G_prod *= util.db_to_lin(gains[i])
        F += (nf_abs - 1) / G_prod

    F_db = 10 * log10(F)
    util.log_result("Rx noise figure", "{:.2f} dB".format(F_db))
    return F_db


def noise_fig_to_noise_temp(nf):
    """Convert noise figure to the effective input-noise temperature

    Note the noise figure is always referenced to a noise source at the
    standard noise temperature of T0 = 290 K. In contrast, the noise
    temperature is independent of the temperature of the noise source.

    Args:
        nf : Noise figure in dB.

    Returns:
        Noise temperature in K.

    """
    nf_abs = util.db_to_lin(nf)

    # Using Equation 8-30b in [1]:
    Te = T0 * (nf_abs - 1)

    return Te


def noise_temp_to_noise_fig(Te):
    """Convert an effective input-noise temperature to a noise figure in dB

    Args:
        Te : Noise temperature in K.

    Returns:
        Noise figure in dB.

    """
    # Noise factor
    nf_abs = 1 + Te / T0
    # Return the noise figure
    return 10 * log10(nf_abs)


def rx_sys_noise_temp(Tar, Te):
    """Compute the receiver system noise temperature

    The receiver noise temperature is the sum of the effective input-noise
    temperature (Te) of the entire receiver seen as a blackbox and the antenna
    noise temperature (Tar). The Te term represents the noise introduced by the
    cascaded linear components (e.g., the LNB, the coax line, and the radio
    interface) of the receiver. The Tar component, in turn, is the noise
    captured by the antenna due to the received cosmic noise and Earth
    blackbody radiation. The simplified model is as follows:

    .. code-block::

         Rx Antenna (Tar) -----> Sum ----> Noise-free Gain Stage ---> Detector
                                  ^
                                  |
                                  |
                           Receiver Noise (Te)

    Note that this is peculiar because the antenna is not treated as another
    cascaded device within the receiver. Instead, the antenna adds to the
    cascaded devices. See Figure 8-24 in [1].

    As explained in [2], around equation 4.39, this is just a convenient choice
    in terms of where the effective input-noise temperature is observed. Note
    that an equivalent (or effective) input noise temperature represents the
    thermodynamic temperature of a noisy resistor, connected to the input of a
    noiseless two-port element, which gives the same output noise power as the
    noisy but otherwise equivalent two-port element, with an ideal noiseless
    source at its input [2]. Hence, if we group the entire receiver into a
    single equivalent two-port element, the Te term represents the noise
    generated by the entire receiver, which has power k*Te*B. When combined to
    the noise collected by the antenna, one obtaines the total system noise
    temperature.

    Args:
        Tar : Antenna noise temperature in K.
        Te  : Effective input-noise temperature in K.

    Returns:
        Receiver system noise temperature in K.

    """

    # Equation 8-41 from [1], or 4.39 from [2]:
    Tsyst = Tar + Te
    util.log_result("System noise temperature", "{:.2f} K".format(Tsyst))
    return Tsyst


def g_over_t(rx_ant_gain_db, T_sys_db):
    """Compute the G/T ratio in dB/K

    Compute the ratio between the Rx antenna gain and the receiver system noise
    temperature, known as G/T. This ratio determines the quality of the
    satellite receiving system. The CNR is proportional to it, and the G/T
    ratio represents the part of the CNR that can be improved based on the
    receiver hardware alone. The other terms of the CNR are EIRP, slant range,
    noise bandwidth, and downlink frequency, which are usually fixed for a
    specific location and service.

    Args:
        rx_ant_gain_db : Receiver antenna gain in dB.
        T_sys_db       : Receiver system noise temperature in dBK.

    Returns:
        G/T in decibels with units of dBK^−1 (or dB/K).

    """
    g_over_t_db = rx_ant_gain_db - T_sys_db
    util.log_result("G/T", "{:.2f} dB/K".format(g_over_t_db))
    return g_over_t_db


def rx_flux_density(eirp_db, distance, atm_loss_db=0):
    """Compute the received flux density in in dBW/m^2

    Note the Rx flux density refers the incident power density (power per area)
    arriving at the receiving antenna, i.e., at the antenna's input. In
    contrast, the value computed by function "rx_power" refers to the power
    collected by the antenna, i.e., at the antenna's output.

    Args:
        eirp_db     : EIRP in dBW.
        distance    : Distance between the Tx and Rx in meters.
        atm_loss_db : Total atmospheric loss in dB.

    Returns:
        Received flux density in dBW/m^2.

    Note:

        The flux density unit of :math:`\\text{dBW}/m^2` should be interpreted
        as :math:`10\\log10(\\frac{W}{m^2})`, and not
        :math:`\\frac{10\\log10(W)}{m^2}`. With a flux density F in
        :math:`\\text{dBW}/m^2` and a given area A in :math:`m^2`, the Rx power
        can be obtained by :math:`F + 10\\log10(A)`. Equivalently, if the area
        is given as :math:`A_{db}` in :math:`\\text{dBm}^2` (decibels greater
        than 1 square meter), the Rx power can be obtained :math:`F + A_{db}`.

    """
    Pt = eirp_db - atm_loss_db  # Tx power minus atmospheric losses
    flux = util.db_to_lin(Pt) / (4 * pi * distance**2)  # in W / m^2
    flux_dbw_m2 = util.lin_to_db(flux)
    util.log_result("Rx flux density", "{:.2f} dBW/m2".format(flux_dbw_m2))
    return flux_dbw_m2


def rx_power(eirp_db,
             path_loss_db,
             rx_ant_gain_db,
             atm_loss_db=0,
             mispointing_db=0):
    """Compute the received carrier power in dBW

    Args:
        eirp_db        : EIRP in dBW.
        path_loss_db   : Free-space path loss in dB.
        rx_ant_gain_db : Receiver antenna gain in dB.
        atm_loss_db    : Total atmospheric loss in dB.
        mispointing_db : Antenna mispointing loss in dB.

    Returns:
        Received power in dBW.

    """
    P_rx_dbw = eirp_db - path_loss_db - atm_loss_db + rx_ant_gain_db \
        - mispointing_db

    util.log_result("Rx Power (LNB Input)",
                    "{:.2f} dBm".format(util.dbw_to_dbm(P_rx_dbw)))

    return P_rx_dbw


def noise_power(T_sys_db, bw):
    """Compute the receiver noise power in dBW

    According to Equation 8-40 in [1], the noise power is given by `N =
    k*Tsyst*bw`, where k is the Boltzmann constant, Tsyst is the receiver
    system noise temperature (in absolute units) and bw is the IF equivalent
    bandwidth in Hz.

    Args:
        T_sys_db : Receiver system noise temperature in dBK.
        bw       : Nominal signal bandwidth (also known as noise bandwidth).

    Returns:
        Receiver noise power in dBW.

    """

    k_db = -228.6  # Boltzmann’s constant (of 1.38e-23) in dB
    bw_db = 10 * log10(bw)

    # Noise power:
    N_dbw = k_db + T_sys_db + bw_db
    util.log_result("Noise Power", "{:.2f} dBm".format(util.dbw_to_dbm(N_dbw)))

    return N_dbw


def spectral_density(power_dbw, bw, label=""):
    """Compute the power spectral density (PSD) of a given signal

    Assume the signal has a constant power spectral density over the given
    bandwidth. Equivalently, assume the signal behaves as white noise, i.e.,
    with uncorrelated zero-mean samples in time.

    Args:
        power_dbw : Signal power in dBW.
        bw        : Bandwidth in Hz.
        label     : Optional label used for logging.

    Returns:
        PSD in dBW/Hz.

    """
    bw_db = util.lin_to_db(bw)
    psd_dbw_hz = power_dbw - bw_db
    util.log_result("{} PSD".format(label),
                    "{:.2f} dBm/Hz".format(util.dbw_to_dbm(psd_dbw_hz)))
    return psd_dbw_hz


def cnr(P_rx_dbw, N_dbw):
    """Compute the carrier-to-noise ratio (CNR) in dB

    Args:
        P_rx_dbw : Received (carrier) power in dBW
        N_dbw    : Receiver noise power in dBW

    Returns:
        CNR in dB.

    """
    cnr_db = P_rx_dbw - N_dbw
    util.log_result("C/N", "{:.2f} dB".format(cnr_db))

    return cnr_db


def capacity(snr_db, bw):
    """Compute the channel capacity in bps

    Args:
        snr_db : signal-to-noise ratio in dB.
        bw     : nominal bandwidth.

    Returns:
        Capacity in bits per second (bps).

    """
    snr = util.db_to_lin(snr_db)
    c = bw * log2(1 + snr)
    util.log_result("Capacity", util.format_rate(c))
    return c


def carrier_to_asi_ratio(rx_dish, long_separation, asi_eirp_ratio):
    """Compute the carrier to adjacent satellite interference (ASI) ratio

    Args:
        rx_dish (Antenna): Rx antenna object.
        long_separation (float): Longitudinal separation in degrees between
            the wanted satellite and the adjacent interferer(s).
        asi_eirp_ratio (float): Ratio between the aggregate adjacent downlink
            EIRP and the wanted signal's EIRP.

    Returns:
        (float): Carrier-to-ASI ratio in dB.

    """
    # Compute the difference between the on-axis (boresight) antenna gain and
    # the off-axis gain at the specified longitudinal separation.
    on_axis_gain = rx_dish.gain_db
    off_axis_gain = rx_dish.off_axis_gain(long_separation)
    delta_gain_db = on_axis_gain - off_axis_gain
    # The adjacent signal is received with the off-axis gain. Hence, the
    # difference between the wanted carrier and the adjacent signal power is
    # delta_gain_db if the adjacent signal has the same EIRP as the wanted
    # signal. Otherwise, if the aggregate adjacent downlink EIRP is greater
    # than the EIRP from the wanted carrier (asi_eirp_ratio > 1), the
    # carrier-to-asi ratio reduces. Similarly, in the opposite case, when
    # asi_eirp_ratio < 1, the carrier-to-asi increases.
    c_to_i_db = delta_gain_db - util.lin_to_db(asi_eirp_ratio)
    util.log_result("ASI C/I ({}° separation)".format(long_separation),
                    "{:.2f} dB".format(c_to_i_db))
    return c_to_i_db


def cnir(cnr, ci):
    """Compute the carrier to noise plus interference ratio

    Based on the reciprocal CNR formula in Eq. (4.42) from [3].

    Args:
        cnr (float): Carrier-to-noise ratio in dB.
        ci (float): Carrier-to-interference ratio in dB.

    Returns:
        (float): Carrier to noise plus interference ratio C/(N+I) in dB.

    """
    cnr_lin = util.db_to_lin(cnr)
    ci_lin = util.db_to_lin(ci)
    cnir_lin = 1 / ((1 / cnr_lin) + (1 / ci_lin))
    cnir_db = util.lin_to_db(cnir_lin)
    util.log_result("C/(N+I)", "{:.2f} dB".format(cnir_db))
    return cnir_db
