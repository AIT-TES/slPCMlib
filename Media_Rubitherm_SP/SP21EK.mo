within slPCMlib.Media_Rubitherm_SP;
package SP21EK "Rubitherm SP21EK; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP21EK";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.731500000000000e+02+1.500000000000000e+01, 2.731500000000000e+02+2.600000000000000e+01}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.731500000000000e+02+1.400000000000000e+01, 2.731500000000000e+02+2.300000000000000e+01}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.700000000000000e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+1.500000000000000e+01
             "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy  href = 0.0
             "reference enthalpy at Tref";

  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  protected
    constant Integer len_x =    data_H.len_x;
    constant Real data_x[:] =   data_H.data_x;
    constant Real data_y[:] =   data_H.data_y;
    constant Real m_k[:] =      data_H.m_k;
    constant Real iy_start[:] = data_H.iy_start;
    constant Real iy_scaler =   data_H.iy_scaler;
  algorithm
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15,
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
  protected
    constant Integer len_x =    data_C.len_x;
    constant Real data_x[:] =   data_C.data_x;
    constant Real data_y[:] =   data_C.data_y;
    constant Real m_k[:] =      data_C.m_k;
    constant Real iy_start[:] = data_C.iy_start;
    constant Real iy_scaler =   data_C.iy_scaler;
  algorithm
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15,
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complSolidification;
  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01};
    constant Real[12] data_y =   {0.000000000000000e+00, 6.369426751592357e-03, 1.910828025477707e-02, 2.547770700636943e-02, 4.458598726114649e-02, 1.082802547770701e-01, 2.484076433121019e-01, 2.993630573248408e-01, 1.656050955414013e-01, 7.643312101910828e-02, 6.369426751592357e-03, 0.000000000000000e+00};
    constant Real[12] m_k =      {6.369426751592357e-03, 9.554140127388535e-03, 9.554140127388536e-03, 1.273885350318471e-02, 4.140127388535032e-02, 1.019108280254777e-01, 9.554140127388534e-02, 0.000000000000000e+00, -1.114649681528662e-01, -7.961783439490445e-02, -1.884829153819383e-02, -3.141381923032305e-03};
    constant Real[12] iy_start = {0.000000000000000e+00, 2.917008668534110e-03, 1.564577376759205e-02, 3.765593008471306e-02, 7.027339065104901e-02, 1.416075117270195e-01, 3.203405883262914e-01, 6.019645161429481e-01, 8.435458704188185e-01, 9.618173127975651e-01, 9.981256812826937e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.992080602760478e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01};
    constant Real[10] data_y =   {0.000000000000000e+00, 7.430329047995139e-03, 2.240140672350838e-02, 2.273430636969565e-02, 5.279346931568574e-02, 9.815998930085812e-02, 2.407949632524718e-01, 5.044348515524332e-01, 5.125068443735692e-02, 0.000000000000000e+00};
    constant Real[10] m_k =      {7.430329047995139e-03, 1.120070336175419e-02, 4.491644770021025e-04, 8.919926246802399e-04, 3.771284146558124e-02, 9.400074696839306e-02, 2.031374311257875e-01, 0.000000000000000e+00, -1.506728456950410e-01, -3.061678408710208e-02};
    constant Real[10] iy_start = {0.000000000000000e+00, 3.390217664930116e-03, 1.915207270092999e-02, 4.161181620172871e-02, 7.619764240414463e-02, 1.467599884569030e-01, 3.066358398366995e-01, 6.949476877654787e-01, 9.844286998335743e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.968394281250451e-01;
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.500000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.400000000000000e+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 5.000000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 5.000000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------

annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP21EK</strong>.<br><br>
  Information taken from: data sheet - last access 02.12.2019.<br><br>
  It also contains the phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  These functions are modelled by piece-wise splines using <strong>pchip</strong> method,
  see also 
  <blockquote>
  <p>
  Barz, T., Krämer, J., & Emhofer, J. (2020). Identification of Phase
  Fraction–Temperature Curves from Heat Capacity Data for Numerical
  Modeling of Heat Transfer in Commercial Paraffin Waxes.
  Energies, 13(19), 5149.
  <a href>doi.org/10.3390/en13195149</a>.
  </p>
  </blockquote>
  </p></html>",
  revisions="<html>
  <ul>
  <li>file creation date: 18-Jul-2022  </ul>
  </html>"));
end SP21EK;
