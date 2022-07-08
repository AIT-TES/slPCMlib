within slPCMlib.Media_Rubitherm_SP;
package SP90 "Rubitherm SP90; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP90";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.571500000000000e+02, 3.651500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.551500000000000e+02, 3.631500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.500000000000000e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+2.500000000000000e+01
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
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {8.400000000000000e+01, 8.500000000000000e+01, 8.600000000000000e+01, 8.700000000000000e+01, 8.800000000000000e+01, 8.900000000000000e+01, 9.000000000000000e+01, 9.100000000000000e+01, 9.200000000000000e+01};
    constant Real[9] data_y =   {0.000000000000000e+00, 1.006506966623132e-02, 3.060043150823838e-02, 5.196255274048978e-02, 1.144449975913506e-01, 3.002238390783484e-01, 4.733520491611969e-01, 1.935106025416807e-02, 0.000000000000000e+00};
    constant Real[9] m_k =      {1.006506966623132e-02, 1.530021575411919e-02, 2.094874153712922e-02, 4.192228304155613e-02, 1.241306431689293e-01, 1.794535257849231e-01, 0.000000000000000e+00, -5.786010568208064e-02, -4.730746992013904e-03};
    constant Real[9] iy_start = {0.000000000000000e+00, 4.590612504044992e-03, 2.442819317548207e-02, 6.391320575723397e-02, 1.401722579995789e-01, 3.426467881142957e-01, 7.438944604665799e-01, 9.947583792686996e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.987685336576415e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {8.200000000000000e+01, 8.300000000000000e+01, 8.400000000000000e+01, 8.500000000000000e+01, 8.600000000000000e+01, 8.700000000000000e+01, 8.800000000000000e+01, 8.900000000000000e+01, 9.000000000000000e+01};
    constant Real[9] data_y =   {0.000000000000000e+00, 9.009009009009009e-03, 9.009009009009009e-03, 2.702702702702703e-02, 9.009009009009009e-02, 1.801801801801802e-01, 5.495495495495496e-01, 1.351351351351351e-01, 0.000000000000000e+00};
    constant Real[9] m_k =      {9.009009009009009e-03, 0.000000000000000e+00, 0.000000000000000e+00, 4.054054054054054e-02, 7.657657657657657e-02, 2.297297297297298e-01, 0.000000000000000e+00, -2.747747747747748e-01, -1.351351351351351e-01};
    constant Real[9] iy_start = {0.000000000000000e+00, 5.192878338278932e-03, 1.409495548961424e-02, 2.856083086053413e-02, 8.345697329376855e-02, 2.043768545994065e-01, 5.838278931750742e-01, 9.447329376854601e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.881305637982196e-01;
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.700000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.650000000000000e+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 6.000000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 6.000000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------

annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP90</strong>.<br><br>
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
  <li>file creation date: 07-Jul-2022  </ul>
  </html>"));
end SP90;
