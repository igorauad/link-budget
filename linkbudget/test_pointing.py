import os
import shutil
import unittest

from skyfield.api import load, wgs84

from . import pointing
from .constants import GEOSYNC_ORBIT


class TestPointing(unittest.TestCase):

    def test_gso_look_angles(self):
        """Test look angle calculations for GSO satellites

        Test Rx stations at all four quadrants (SW, NW, NE, SE). The look
        angles are compared to those obtained with the dish pointer tool
        (https://www.dishpointer.com).

        SW Station:
        1) Sao Paulo (Latitude: -23.5505°, Longitude: -46.6333°)
           Satellite: Eutelsat 113 (113° West)

        NW Stations:
        1) Washington DC (Latitude: 38.9072°, Longitude: -77.0369°)
           Satellite: Eutelsat 113 (113° West)
        2) Los Angeles (Latitude: 34.0522°, Longitude: -118.2437°)
           Satellite: Galaxy 18 (123° West)

        NE Station:
        1) Berlin (Latitude: 52.5200°, Longitude: 13.4050°)
           Satellite: Telstar 11N (37.5° West)

        SE Stations:
        1) Sydney (Latitude: -33.8688°, Longitude: 151.2093°)
           Satellite: Telstar 18V (138° East)
        2) Cape Town (Latitude: -33.9249°, Longitude: 18.4241°)
           Satellite: Telstar 11N (37.5° West)

        """
        setup_info = [
            {
                'sat_long': -113,
                'rx_lat': -23.5505,
                'rx_long': -46.6333,
                'elevation': 13.1,
                'azimuth': 279.9,
                'distance': 40260,
                'lnb_skew': -64.6
            },
            {
                'sat_long': -113,
                'rx_lat': 38.9072,
                'rx_long': -77.0369,
                'elevation': 31.6,
                'azimuth': 229.1,
                'distance': 38469,
                'lnb_skew': 36
            },
            {
                'sat_long': -123,
                'rx_lat': 34.0522,
                'rx_long': -118.2437,
                'elevation': 50.1,
                'azimuth': 188.4,
                'distance': 37077,
                'lnb_skew': 7
            },
            {
                'sat_long': -37.5,
                'rx_lat': 52.5200,
                'rx_long': 13.4050,
                'elevation': 14.1,
                'azimuth': 237.2,
                'distance': 40151,
                'lnb_skew': 30.8
            },
            {
                'sat_long': 138,
                'rx_lat': -33.8688,
                'rx_long': 151.2093,
                'elevation': 48.1,
                'azimuth': 337.2,
                'distance': 37202,
                'lnb_skew': -18.8
            },
            {
                'sat_long': -37.5,
                'rx_lat': -33.9249,
                'rx_long': 18.4241,
                'elevation': 19.5,
                'azimuth': 290.7,
                'distance': 39604,
                'lnb_skew': -50.9
            },
        ]

        # Geostationary latitude and altitude
        sat_lat = 0
        sat_alt = GEOSYNC_ORBIT

        for info in setup_info:
            elevation, azimuth, slant_range = pointing.look_angles(
                info['sat_long'], sat_lat, sat_alt, info['rx_long'],
                info['rx_lat'])
            slant_range_km = slant_range / 1e3
            pol_angle = pointing.polarization_angle(info['sat_long'], sat_lat,
                                                    info['rx_long'],
                                                    info['rx_lat'])

            self.assertAlmostEqual(elevation, info['elevation'], delta=0.1)
            self.assertAlmostEqual(azimuth, info['azimuth'], delta=0.1)
            # The distance sometimes differs by a few km. Maybe dishpointer
            # considers the actual longitude of the satellite based on
            # ephemeris data instead of the nominal.
            self.assertAlmostEqual(slant_range_km, info['distance'], delta=12)
            self.assertAlmostEqual(pol_angle, info['lnb_skew'], delta=0.05)

    def test_ngso_look_angles(self):
        """Test look angle calculations for NGSO satellites

        Tested by comparison to the results obtained with Skyfield.

        """
        # Rx in Washington DC
        rx_lat = 38.9072
        rx_long = -77.0369
        rx_height = 0
        observer_pos = wgs84.latlon(rx_lat, rx_long, rx_height)

        # NGSO Satellite: Starlink-3699 covering DC at 2023-03-17T18:31:45
        ts = load.timescale()
        sat_long, sat_lat, sat_alt = pointing.get_sat_pos_by_tle(
            'STARLINK-3699')
        t = ts.utc(2023, 3, 17, 18, 31, 45)
        sat_pos = wgs84.latlon(sat_lat, sat_long, sat_alt)

        difference = sat_pos - observer_pos
        topocentric = difference.at(t)
        skyfield_el, skyfield_az, _ = topocentric.altaz()

        el, az, _ = pointing.look_angles(sat_long, sat_lat, sat_alt, rx_long,
                                         rx_lat, rx_height)

        self.assertAlmostEqual(el, skyfield_el.degrees, delta=0.25)
        self.assertAlmostEqual(az, skyfield_az.degrees, delta=0.25)
        # TODO: the slant range differs significantly. Understand why.

    def test_get_sat_pos_by_tle(self):
        """Test satellite position through TLE-based prediction"""

        # Basic query using the default download directory
        lon, lat, alt = pointing.get_sat_pos_by_tle('DIRECTV 9S')
        self.assertAlmostEqual(lon, -101, delta=0.3)
        self.assertAlmostEqual(lat, 0, delta=0.1)
        self.assertAlmostEqual(alt, GEOSYNC_ORBIT, delta=20e3)

        # The dataset should be saved
        expected_file = os.path.join(pointing.get_default_tle_dataset_dir(),
                                     'tle-name-DIRECTV-9S.txt')
        self.assertTrue(os.path.exists(expected_file))
        os.remove(expected_file)

        # Same thing but specifying another download directory
        tle_download_dir = '/tmp/link-budget-test/'
        lon, lat, alt = pointing.get_sat_pos_by_tle('DIRECTV 9S',
                                                    save_dir=tle_download_dir)
        expected_file = os.path.join(tle_download_dir,
                                     'tle-name-DIRECTV-9S.txt')
        self.assertTrue(os.path.exists(expected_file))

        # If the no_save option is specified but the dataset is already
        # available, the dataset must be kept (not removed in the end)
        lon, lat, alt = pointing.get_sat_pos_by_tle('DIRECTV 9S',
                                                    save_dir=tle_download_dir,
                                                    no_save=True)
        self.assertTrue(os.path.exists(expected_file))

        # On the other hand, if no_save is True and the dataset does not exist
        # yet, it should be removed in the end
        os.remove(expected_file)  # pretend the following is the first download
        lon, lat, alt = pointing.get_sat_pos_by_tle('DIRECTV 9S',
                                                    save_dir=tle_download_dir,
                                                    no_save=True)
        self.assertFalse(os.path.exists(expected_file))

        # Query under a specific group
        lon, lat, alt = pointing.get_sat_pos_by_tle('DIRECTV 9S',
                                                    group='active',
                                                    save_dir=tle_download_dir)
        self.assertAlmostEqual(lon, -101, delta=0.3)
        self.assertAlmostEqual(lat, 0, delta=0.1)
        self.assertAlmostEqual(alt, GEOSYNC_ORBIT, delta=20e3)

        # The group query downloads the entire group's dataset (prefixed with
        # 'tle-group'), unlike the name-only query, which downloads the TLE
        # entries matching the given name (saved with 'tle-name' prefix).
        expected_file = os.path.join(tle_download_dir, 'tle-group-active.txt')
        self.assertTrue(os.path.exists(expected_file))

        # Invalid satellite name
        with self.assertRaises(ValueError):
            pointing.get_sat_pos_by_tle('EUTELSAT 113',
                                        group='active',
                                        save_dir=tle_download_dir)

        # The same saving behavior must apply even when the satellite is not
        # found. In the above command, the active group's TLE dataset will be
        # downloaded even though the satellite name is invalid. If we add the
        # no_save option, the dataset should be removed in the end, except if
        # the dataset already existed before due to a previous query.
        with self.assertRaises(ValueError):
            pointing.get_sat_pos_by_tle('EUTELSAT 113',
                                        group='active',
                                        save_dir=tle_download_dir,
                                        no_save=True)
        self.assertTrue(
            os.path.exists(expected_file))  # still there (existed before)

        os.remove(expected_file)  # pretend the following is the first download
        with self.assertRaises(ValueError):
            pointing.get_sat_pos_by_tle('EUTELSAT 113',
                                        group='active',
                                        save_dir=tle_download_dir,
                                        no_save=True)
        self.assertFalse(os.path.exists(expected_file))  # removed in the end

        shutil.rmtree(tle_download_dir)
