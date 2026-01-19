"""CalculateUTCI."""

import numpy as np
import math

from ..datetime_utils import as_datetime
from ..logs import logger
from ..toolbox import FileManager
from .base import Task


class Utci():
    """Calculate UTCI form the operational procedure."""

    def __init__(self, values):
        self.values = values

    def calc(self):
        """ Process self.values and calculate UTCI."""
        # Self.values are in order [t2, mrt, r2, u10, v10]
        data_all = [self.values[0], self.values[1], self.values[2]/100, self.values[3], self.values[4]]

        # Calculate wind speed
        u = data_all[3]
        v = data_all[4]
        va = np.hypot(u,v)

        # Set up physical constants
        RKBOL = 1.380658E-23
        RNAVO = 6.0221367E+23
        R     = RNAVO*RKBOL
        RMV   = 18.0153
        RV    = 1000.*R/ RMV
        RCPV  = 4.*RV
        RCW   = 4218.
        RLVTT = 2.5008E+06
        RESTT = 611.14
        RTT   = 273.16
        RGAMW = (RCW-RCPV)/RV
        RBETW = RLVTT/RV+RGAMW*RTT
        RALPW = math.log(RESTT) + RBETW/ RTT + RGAMW * math.log(RTT)

        # Calculate saturation pressure over water [Pa], consistently with model
        zesat = np.exp( RALPW - RBETW/ data_all[0] - RGAMW * np.log(data_all[0]) )

        # Change of units
        Ta = data_all[0] - 273.15           # [deg C]
        D_Tmrt = data_all[1] - data_all[0]  # [deg C]
        Pa = data_all[2] * zesat/ 1000.     # [kPa]
        #va is in [m/s]

        # Clip inputs to valid fitting ranges - this is important since
        # 6th order polynomial fit is used!
        Ta     = np.clip(Ta, -50, 50)        # -50 to 50 deg C
        D_Tmrt = np.clip(D_Tmrt, -30, 70)   # -30 to 70 deg C
        Pa     = np.clip(Pa, 0.0, 5.)       #   0 to  5 kPa
        va     = np.clip(va, 0.5, 17.)      # 0.5 to 17 m/s
                
        # Inlined UTCI temperature function [deg C]
        #!!! DO NOT TOUCH !!!
                
        zutci = Ta + \
          (  6.07562052E-01 )   + \
          ( -2.27712343E-02 ) * Ta + \
          (  8.06470249E-04 ) * np.prod([Ta,Ta], 0) + \
          ( -1.54271372E-04 ) * np.prod([Ta,Ta,Ta], 0) + \
          ( -3.24651735E-06 ) * np.prod([Ta,Ta,Ta,Ta], 0) + \
          (  7.32602852E-08 ) * np.prod([Ta,Ta,Ta,Ta,Ta], 0) + \
          (  1.35959073E-09 ) * np.prod([Ta,Ta,Ta,Ta,Ta,Ta], 0) + \
          ( -2.25836520E+00 ) * va + \
          (  8.80326035E-02 ) * np.prod([Ta,va], 0) + \
          (  2.16844454E-03 ) * np.prod([Ta,Ta,va], 0) + \
          ( -1.53347087E-05 ) * np.prod([Ta,Ta,Ta,va], 0) + \
          ( -5.72983704E-07 ) * np.prod([Ta,Ta,Ta,Ta,va], 0) + \
          ( -2.55090145E-09 ) * np.prod([Ta,Ta,Ta,Ta,Ta,va], 0) + \
          ( -7.51269505E-01 ) * np.prod([va,va], 0) + \
          ( -4.08350271E-03 ) * np.prod([Ta,va,va], 0) + \
          ( -5.21670675E-05 ) * np.prod([Ta,Ta,va,va], 0) + \
          (  1.94544667E-06 ) * np.prod([Ta,Ta,Ta,va,va], 0) + \
          (  1.14099531E-08 ) * np.prod([Ta,Ta,Ta,Ta,va,va], 0) + \
          (  1.58137256E-01 ) * np.prod([va,va,va], 0) + \
          ( -6.57263143E-05 ) * np.prod([Ta,va,va,va], 0) + \
          (  2.22697524E-07 ) * np.prod([Ta,Ta,va,va,va], 0) + \
          ( -4.16117031E-08 ) * np.prod([Ta,Ta,Ta,va,va,va], 0) + \
          ( -1.27762753E-02 ) * np.prod([va,va,va,va], 0) + \
          (  9.66891875E-06 ) * np.prod([Ta,va,va,va,va], 0) + \
          (  2.52785852E-09 ) * np.prod([Ta,Ta,va,va,va,va], 0) + \
          (  4.56306672E-04 ) * np.prod([va,va,va,va,va], 0) + \
          ( -1.74202546E-07 ) * np.prod([Ta,va,va,va,va,va], 0) + \
          ( -5.91491269E-06 ) * np.prod([va,va,va,va,va,va], 0) + \
          (  3.98374029E-01 ) * D_Tmrt + \
          (  1.83945314E-04 ) * np.prod([Ta,D_Tmrt], 0) + \
          ( -1.73754510E-04 ) * np.prod([Ta,Ta,D_Tmrt], 0) + \
          ( -7.60781159E-07 ) * np.prod([Ta,Ta,Ta,D_Tmrt], 0) + \
          (  3.77830287E-08 ) * np.prod([Ta,Ta,Ta,Ta,D_Tmrt], 0) + \
          (  5.43079673E-10 ) * np.prod([Ta,Ta,Ta,Ta,Ta,D_Tmrt], 0) + \
          ( -2.00518269E-02 ) * np.prod([va,D_Tmrt], 0) + \
          (  8.92859837E-04 ) * np.prod([Ta,va,D_Tmrt], 0) + \
          (  3.45433048E-06 ) * np.prod([Ta,Ta,va,D_Tmrt], 0) + \
          ( -3.77925774E-07 ) * np.prod([Ta,Ta,Ta,va,D_Tmrt], 0) + \
          ( -1.69699377E-09 ) * np.prod([Ta,Ta,Ta,Ta,va,D_Tmrt], 0) + \
          (  1.69992415E-04 ) * np.prod([va,va,D_Tmrt], 0) + \
          ( -4.99204314E-05 ) * np.prod([Ta,va,va,D_Tmrt], 0) + \
          (  2.47417178E-07 ) * np.prod([Ta,Ta,va,va,D_Tmrt], 0) + \
          (  1.07596466E-08 ) * np.prod([Ta,Ta,Ta,va,va,D_Tmrt], 0) + \
          (  8.49242932E-05 ) * np.prod([va,va,va,D_Tmrt], 0) + \
          (  1.35191328E-06 ) * np.prod([Ta,va,va,va,D_Tmrt], 0) + \
          ( -6.21531254E-09 ) * np.prod([Ta,Ta,va,va,va,D_Tmrt], 0) + \
          ( -4.99410301E-06 ) * np.prod([va,va,va,va,D_Tmrt], 0) + \
          ( -1.89489258E-08 ) * np.prod([Ta,va,va,va,va,D_Tmrt], 0) + \
          (  8.15300114E-08 ) * np.prod([va,va,va,va,va,D_Tmrt], 0) + \
          (  7.55043090E-04 ) * np.prod([D_Tmrt,D_Tmrt], 0) + \
          ( -5.65095215E-05 ) * np.prod([Ta,D_Tmrt,D_Tmrt], 0) + \
          ( -4.52166564E-07 ) * np.prod([Ta,Ta,D_Tmrt,D_Tmrt], 0) + \
          (  2.46688878E-08 ) * np.prod([Ta,Ta,Ta,D_Tmrt,D_Tmrt], 0) + \
          (  2.42674348E-10 ) * np.prod([Ta,Ta,Ta,Ta,D_Tmrt,D_Tmrt], 0) + \
          (  1.54547250E-04 ) * np.prod([va,D_Tmrt,D_Tmrt], 0) + \
          (  5.24110970E-06 ) * np.prod([Ta,va,D_Tmrt,D_Tmrt], 0) + \
          ( -8.75874982E-08 ) * np.prod([Ta,Ta,va,D_Tmrt,D_Tmrt], 0) + \
          ( -1.50743064E-09 ) * np.prod([Ta,Ta,Ta,va,D_Tmrt,D_Tmrt], 0) + \
          ( -1.56236307E-05 ) * np.prod([va,va,D_Tmrt,D_Tmrt], 0) + \
          ( -1.33895614E-07 ) * np.prod([Ta,va,va,D_Tmrt,D_Tmrt], 0) + \
          (  2.49709824E-09 ) * np.prod([Ta,Ta,va,va,D_Tmrt,D_Tmrt], 0) + \
          (  6.51711721E-07 ) * np.prod([va,va,va,D_Tmrt,D_Tmrt], 0) + \
          (  1.94960053E-09 ) * np.prod([Ta,va,va,va,D_Tmrt,D_Tmrt], 0) + \
          ( -1.00361113E-08 ) * np.prod([va,va,va,va,D_Tmrt,D_Tmrt], 0) + \
          ( -1.21206673E-05 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          ( -2.18203660E-07 ) * np.prod([Ta,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  7.51269482E-09 ) * np.prod([Ta,Ta,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  9.79063848E-11 ) * np.prod([Ta,Ta,Ta,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  1.25006734E-06 ) * np.prod([va,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          ( -1.81584736E-09 ) * np.prod([Ta,va,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          ( -3.52197671E-10 ) * np.prod([Ta,Ta,va,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          ( -3.36514630E-08 ) * np.prod([va,va,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  1.35908359E-10 ) * np.prod([Ta,va,va,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  4.17032620E-10 ) * np.prod([va,va,va,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          ( -1.30369025E-09 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  4.13908461E-10 ) * np.prod([Ta,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  9.22652254E-12 ) * np.prod([Ta,Ta,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          ( -5.08220384E-09 ) * np.prod([va,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          ( -2.24730961E-11 ) * np.prod([Ta,va,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  1.17139133E-10 ) * np.prod([va,va,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  6.62154879E-10 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  4.03863260E-13 ) * np.prod([Ta,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  1.95087203E-12 ) * np.prod([va,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          ( -4.73602469E-12 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt], 0) + \
          (  5.12733497E+00 ) * Pa + \
          ( -3.12788561E-01 ) * np.prod([Ta,Pa], 0) + \
          ( -1.96701861E-02 ) * np.prod([Ta,Ta,Pa], 0) + \
          (  9.99690870E-04 ) * np.prod([Ta,Ta,Ta,Pa], 0) + \
          (  9.51738512E-06 ) * np.prod([Ta,Ta,Ta,Ta,Pa], 0) + \
          ( -4.66426341E-07 ) * np.prod([Ta,Ta,Ta,Ta,Ta,Pa], 0) + \
          (  5.48050612E-01 ) * np.prod([va,Pa], 0) + \
          ( -3.30552823E-03 ) * np.prod([Ta,va,Pa], 0) + \
          ( -1.64119440E-03 ) * np.prod([Ta,Ta,va,Pa], 0) + \
          ( -5.16670694E-06 ) * np.prod([Ta,Ta,Ta,va,Pa], 0) + \
          (  9.52692432E-07 ) * np.prod([Ta,Ta,Ta,Ta,va,Pa], 0) + \
          ( -4.29223622E-02 ) * np.prod([va,va,Pa], 0) + \
          (  5.00845667E-03 ) * np.prod([Ta,va,va,Pa], 0) + \
          (  1.00601257E-06 ) * np.prod([Ta,Ta,va,va,Pa], 0) + \
          ( -1.81748644E-06 ) * np.prod([Ta,Ta,Ta,va,va,Pa], 0) + \
          ( -1.25813502E-03 ) * np.prod([va,va,va,Pa], 0) + \
          ( -1.79330391E-04 ) * np.prod([Ta,va,va,va,Pa], 0) + \
          (  2.34994441E-06 ) * np.prod([Ta,Ta,va,va,va,Pa], 0) + \
          (  1.29735808E-04 ) * np.prod([va,va,va,va,Pa], 0) + \
          (  1.29064870E-06 ) * np.prod([Ta,va,va,va,va,Pa], 0) + \
          ( -2.28558686E-06 ) * np.prod([va,va,va,va,va,Pa], 0) + \
          ( -3.69476348E-02 ) * np.prod([D_Tmrt,Pa], 0) + \
          (  1.62325322E-03 ) * np.prod([Ta,D_Tmrt,Pa], 0) + \
          ( -3.14279680E-05 ) * np.prod([Ta,Ta,D_Tmrt,Pa], 0) + \
          (  2.59835559E-06 ) * np.prod([Ta,Ta,Ta,D_Tmrt,Pa], 0) + \
          ( -4.77136523E-08 ) * np.prod([Ta,Ta,Ta,Ta,D_Tmrt,Pa], 0) + \
          (  8.64203390E-03 ) * np.prod([va,D_Tmrt,Pa], 0) + \
          ( -6.87405181E-04 ) * np.prod([Ta,va,D_Tmrt,Pa], 0) + \
          ( -9.13863872E-06 ) * np.prod([Ta,Ta,va,D_Tmrt,Pa], 0) + \
          (  5.15916806E-07 ) * np.prod([Ta,Ta,Ta,va,D_Tmrt,Pa], 0) + \
          ( -3.59217476E-05 ) * np.prod([va,va,D_Tmrt,Pa], 0) + \
          (  3.28696511E-05 ) * np.prod([Ta,va,va,D_Tmrt,Pa], 0) + \
          ( -7.10542454E-07 ) * np.prod([Ta,Ta,va,va,D_Tmrt,Pa], 0) + \
          ( -1.24382300E-05 ) * np.prod([va,va,va,D_Tmrt,Pa], 0) + \
          ( -7.38584400E-09 ) * np.prod([Ta,va,va,va,D_Tmrt,Pa], 0) + \
          (  2.20609296E-07 ) * np.prod([va,va,va,va,D_Tmrt,Pa], 0) + \
          ( -7.32469180E-04 ) * np.prod([D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -1.87381964E-05 ) * np.prod([Ta,D_Tmrt,D_Tmrt,Pa], 0) + \
          (  4.80925239E-06 ) * np.prod([Ta,Ta,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -8.75492040E-08 ) * np.prod([Ta,Ta,Ta,D_Tmrt,D_Tmrt,Pa], 0) + \
          (  2.77862930E-05 ) * np.prod([va,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -5.06004592E-06 ) * np.prod([Ta,va,D_Tmrt,D_Tmrt,Pa], 0) + \
          (  1.14325367E-07 ) * np.prod([Ta,Ta,va,D_Tmrt,D_Tmrt,Pa], 0) + \
          (  2.53016723E-06 ) * np.prod([va,va,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -1.72857035E-08 ) * np.prod([Ta,va,va,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -3.95079398E-08 ) * np.prod([va,va,va,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -3.59413173E-07 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          (  7.04388046E-07 ) * np.prod([Ta,D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -1.89309167E-08 ) * np.prod([Ta,Ta,D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -4.79768731E-07 ) * np.prod([va,D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          (  7.96079978E-09 ) * np.prod([Ta,va,D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          (  1.62897058E-09 ) * np.prod([va,va,D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          (  3.94367674E-08 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -1.18566247E-09 ) * np.prod([Ta,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          (  3.34678041E-10 ) * np.prod([va,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -1.15606447E-10 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,Pa], 0) + \
          ( -2.80626406E+00 ) * np.prod([Pa,Pa], 0) + \
          (  5.48712484E-01 ) * np.prod([Ta,Pa,Pa], 0) + \
          ( -3.99428410E-03 ) * np.prod([Ta,Ta,Pa,Pa], 0) + \
          ( -9.54009191E-04 ) * np.prod([Ta,Ta,Ta,Pa,Pa], 0) + \
          (  1.93090978E-05 ) * np.prod([Ta,Ta,Ta,Ta,Pa,Pa], 0) + \
          ( -3.08806365E-01 ) * np.prod([va,Pa,Pa], 0) + \
          (  1.16952364E-02 ) * np.prod([Ta,va,Pa,Pa], 0) + \
          (  4.95271903E-04 ) * np.prod([Ta,Ta,va,Pa,Pa], 0) + \
          ( -1.90710882E-05 ) * np.prod([Ta,Ta,Ta,va,Pa,Pa], 0) + \
          (  2.10787756E-03 ) * np.prod([va,va,Pa,Pa], 0) + \
          ( -6.98445738E-04 ) * np.prod([Ta,va,va,Pa,Pa], 0) + \
          (  2.30109073E-05 ) * np.prod([Ta,Ta,va,va,Pa,Pa], 0) + \
          (  4.17856590E-04 ) * np.prod([va,va,va,Pa,Pa], 0) + \
          ( -1.27043871E-05 ) * np.prod([Ta,va,va,va,Pa,Pa], 0) + \
          ( -3.04620472E-06 ) * np.prod([va,va,va,va,Pa,Pa], 0) + \
          (  5.14507424E-02 ) * np.prod([D_Tmrt,Pa,Pa], 0) + \
          ( -4.32510997E-03 ) * np.prod([Ta,D_Tmrt,Pa,Pa], 0) + \
          (  8.99281156E-05 ) * np.prod([Ta,Ta,D_Tmrt,Pa,Pa], 0) + \
          ( -7.14663943E-07 ) * np.prod([Ta,Ta,Ta,D_Tmrt,Pa,Pa], 0) + \
          ( -2.66016305E-04 ) * np.prod([va,D_Tmrt,Pa,Pa], 0) + \
          (  2.63789586E-04 ) * np.prod([Ta,va,D_Tmrt,Pa,Pa], 0) + \
          ( -7.01199003E-06 ) * np.prod([Ta,Ta,va,D_Tmrt,Pa,Pa], 0) + \
          ( -1.06823306E-04 ) * np.prod([va,va,D_Tmrt,Pa,Pa], 0) + \
          (  3.61341136E-06 ) * np.prod([Ta,va,va,D_Tmrt,Pa,Pa], 0) + \
          (  2.29748967E-07 ) * np.prod([va,va,va,D_Tmrt,Pa,Pa], 0) + \
          (  3.04788893E-04 ) * np.prod([D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          ( -6.42070836E-05 ) * np.prod([Ta,D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          (  1.16257971E-06 ) * np.prod([Ta,Ta,D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          (  7.68023384E-06 ) * np.prod([va,D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          ( -5.47446896E-07 ) * np.prod([Ta,va,D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          ( -3.59937910E-08 ) * np.prod([va,va,D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          ( -4.36497725E-06 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          (  1.68737969E-07 ) * np.prod([Ta,D_Tmrt,D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          (  2.67489271E-08 ) * np.prod([va,D_Tmrt,D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          (  3.23926897E-09 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt,D_Tmrt,Pa,Pa], 0) + \
          ( -3.53874123E-02 ) * np.prod([Pa,Pa,Pa], 0) + \
          ( -2.21201190E-01 ) * np.prod([Ta,Pa,Pa,Pa], 0) + \
          (  1.55126038E-02 ) * np.prod([Ta,Ta,Pa,Pa,Pa], 0) + \
          ( -2.63917279E-04 ) * np.prod([Ta,Ta,Ta,Pa,Pa,Pa], 0) + \
          (  4.53433455E-02 ) * np.prod([va,Pa,Pa,Pa], 0) + \
          ( -4.32943862E-03 ) * np.prod([Ta,va,Pa,Pa,Pa], 0) + \
          (  1.45389826E-04 ) * np.prod([Ta,Ta,va,Pa,Pa,Pa], 0) + \
          (  2.17508610E-04 ) * np.prod([va,va,Pa,Pa,Pa], 0) + \
          ( -6.66724702E-05 ) * np.prod([Ta,va,va,Pa,Pa,Pa], 0) + \
          (  3.33217140E-05 ) * np.prod([va,va,va,Pa,Pa,Pa], 0) + \
          ( -2.26921615E-03 ) * np.prod([D_Tmrt,Pa,Pa,Pa], 0) + \
          (  3.80261982E-04 ) * np.prod([Ta,D_Tmrt,Pa,Pa,Pa], 0) + \
          ( -5.45314314E-09 ) * np.prod([Ta,Ta,D_Tmrt,Pa,Pa,Pa], 0) + \
          ( -7.96355448E-04 ) * np.prod([va,D_Tmrt,Pa,Pa,Pa], 0) + \
          (  2.53458034E-05 ) * np.prod([Ta,va,D_Tmrt,Pa,Pa,Pa], 0) + \
          ( -6.31223658E-06 ) * np.prod([va,va,D_Tmrt,Pa,Pa,Pa], 0) + \
          (  3.02122035E-04 ) * np.prod([D_Tmrt,D_Tmrt,Pa,Pa,Pa], 0) + \
          ( -4.77403547E-06 ) * np.prod([Ta,D_Tmrt,D_Tmrt,Pa,Pa,Pa], 0) + \
          (  1.73825715E-06 ) * np.prod([va,D_Tmrt,D_Tmrt,Pa,Pa,Pa], 0) + \
          ( -4.09087898E-07 ) * np.prod([D_Tmrt,D_Tmrt,D_Tmrt,Pa,Pa,Pa], 0) + \
          (  6.14155345E-01 ) * np.prod([Pa,Pa,Pa,Pa], 0) + \
          (  1.33374846E-03 ) * np.prod([Ta,Ta,Pa,Pa,Pa,Pa], 0) + \
          (  3.55375387E-03 ) * np.prod([va,Pa,Pa,Pa,Pa], 0) + \
          ( -5.13027851E-04 ) * np.prod([Ta,va,Pa,Pa,Pa,Pa], 0) + \
          (  1.02449757E-04 ) * np.prod([va,va,Pa,Pa,Pa,Pa], 0) + \
          ( -1.48526421E-03 ) * np.prod([D_Tmrt,Pa,Pa,Pa,Pa], 0) + \
          ( -4.11469183E-05 ) * np.prod([Ta,D_Tmrt,Pa,Pa,Pa,Pa], 0) + \
          ( -6.16755931E-02 ) * np.prod([Ta,Pa,Pa,Pa,Pa], 0) + \
          ( -6.80434415E-06 ) * np.prod([va,D_Tmrt,Pa,Pa,Pa,Pa], 0) + \
          ( -9.77675906E-06 ) * np.prod([D_Tmrt,D_Tmrt,Pa,Pa,Pa,Pa], 0) + \
          (  8.82773108E-02 ) * np.prod([Pa,Pa,Pa,Pa,Pa], 0) + \
          ( -3.01859306E-03 ) * np.prod([Ta,Pa,Pa,Pa,Pa,Pa], 0) + \
          (  1.04452989E-03 ) * np.prod([va,Pa,Pa,Pa,Pa,Pa], 0) + \
          (  2.47090539E-04 ) * np.prod([D_Tmrt,Pa,Pa,Pa,Pa,Pa], 0) + \
          (  1.48348065E-03 ) * np.prod([Pa,Pa,Pa,Pa,Pa,Pa], 0)

        calculated_utci_values = zutci + 273.15 # in degrees °K

        return calculated_utci_values

