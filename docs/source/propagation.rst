Propagation Models
=======================================

By default, the link budget tool assumes the location of interest experiences
the atmospheric attenuation given by models from ITU-R P. Recommendations. In
particular, it uses the model implementation imported from the `ITU-Rpy package
<https://itu-rpy.readthedocs.io/>`_. Nevertheless, it is possible to override
the modeled values and specify the atmospheric attenuation directly through
option ``--atmospheric-attenuation`` (see :ref:`Propagation Options`).

The atmospheric attenuation models include the effects of rain, clouds,
atmospheric gases, and tropospheric scintillation. Each of these effects
involves several ITU-R P. recommendations. For instance, the following
recommendations are used to model the rain attenuation:

1. **ITU-R P.838-3:** Specific attenuation model for rain for use in prediction
   methods (1992-1999-2003-2005).
2. **ITU-R P.618-13:** Propagation data and prediction methods required for the
   design of Earth-space telecommunication systems (12/2017).
3. **ITU-R P.1511-2:** Topography for Earth-space propagation modelling
   (08/2019).
4. **ITU-R P.839-4:** Rain height model for prediction methods (09/2013).
5. **ITU-R P.837-7:** Characteristics of precipitation for propagation modelling
   (06/2017).

The sequence used to predict the rain attenuation exceeded over 0.01% of an
average year is as follows:

1. Find the rainfall rate in mm/h exceeded for 0.01% of an average year on the
   specific location of interest using the model and the database provided by
   Recommendation ITU-R P.837-7.
2. Using the carrier frequency, find the coefficients to obtain the rain
   attenuation per kilometer following Recommendation ITU-R P.838-3.
3. Find the mean annual rain height at the given location, i.e., how high is the
   melting layer, where the frozen precipitation turns into liquid
   precipitation. Use the model and the database from Recommendation ITU-R
   P.839-4.
4. Find the earth station's height at the given latitude/longitude by
   interpolating on the database from ITU-R P.1511-2.
5. Using the parameters from steps 2, 3, and 4, as well as the elevation,
   latitude, and carrier frequency, compute the effective slant path subject to
   rain attenuation using the procedure from Section 2.2.1.1 of Recommendation
   ITU-R P.618-13.
6. Finally, multiply the attenuation in dB/km (from step 2) by the effective
   slant path subject to rain attenuation (from step 5) to obtain the final rain
   attenuation in dB.

Distinct procedures apply to the other forms of atmospheric attenuation. Section
2.5 from Recommendation ITU-R P.618-13 outlines the references used for each
computation. Refer to `the ITU-Rpy's docs <https://itu-rpy.readthedocs.io/>`_
for further information.

.. automodule:: linkbudget.propagation
   :members:
