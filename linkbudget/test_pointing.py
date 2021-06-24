import unittest
from . import pointing


class TestPointing(unittest.TestCase):
    def test_look_angles(self):
        """Test Rx stations at all four quadrants (SW, NW, NE, SE)

        Check the elevation angle, azimumth angle (true north), and slant range
        according to https://www.dishpointer.com.

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

        for info in setup_info:
            elevation, azimuth, slant_range = pointing.look_angles(
                info['sat_long'], info['rx_long'], info['rx_lat'])
            slant_range_km = slant_range / 1e3
            pol_angle = pointing.polarization_angle(info['sat_long'],
                                                    info['rx_long'],
                                                    info['rx_lat'])

            self.assertAlmostEqual(elevation, info['elevation'], delta=0.1)
            self.assertAlmostEqual(azimuth, info['azimuth'], delta=0.1)
            # The distance sometimes differs by a few km. Maybe dishpointer
            # considers the actual longitude of the satellite based on
            # ephemeris data instead of the nominal.
            self.assertAlmostEqual(slant_range_km, info['distance'], delta=12)
            self.assertAlmostEqual(pol_angle, info['lnb_skew'], delta=0.05)
