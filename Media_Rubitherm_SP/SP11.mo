within slPCMlib.Media_Rubitherm_SP;
package SP11 "Rubitherm SP11; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP11";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.731500000000000e+02+8.000000000000000e+00, 2.731500000000000e+02+1.500000000000000e+01}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.731500000000000e+02+6.000000000000000e+00, 2.731500000000000e+02+1.200000000000000e+01}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.750000000000000e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+8.000000000000000e+00
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
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01};
    constant Real[8] data_y =   {0.000000000000000e+00, 1.587301587301587e-02, 7.936507936507936e-02, 1.428571428571428e-01, 5.238095238095238e-01, 2.063492063492063e-01, 3.174603174603174e-02, 0.000000000000000e+00};
    constant Real[8] m_k =      {1.587301587301587e-02, 3.968253968253968e-02, 5.232783388521486e-02, 1.831474185982520e-01, 0.000000000000000e+00, -2.460317460317460e-01, -9.102657225923333e-02, -2.800817607976410e-02};
    constant Real[8] iy_start = {0.000000000000000e+00, 5.930693792985688e-03, 5.232630898713522e-02, 1.521706805360480e-01, 4.994962105647947e-01, 8.836733751548675e-01, 9.894172159069550e-01, 9.999999999999999e-01};
    constant Real    iy_scaler = 9.963565572215956e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    7;
    constant Real[7] data_x =   {6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01};
    constant Real[7] data_y =   {0.000000000000000e+00, 3.125000000000000e-02, 1.562500000000000e-02, 1.093750000000000e-01, 2.187500000000000e-01, 6.250000000000000e-01, 0.000000000000000e+00};
    constant Real[7] m_k =      {3.125000000000000e-02, 0.000000000000000e+00, 0.000000000000000e+00, 1.015625000000000e-01, 2.578125000000000e-01, 0.000000000000000e+00, -6.250000000000000e-01};
    constant Real[7] iy_start = {0.000000000000000e+00, 1.728395061728395e-02, 3.950617283950618e-02, 9.074074074074076e-02, 2.339506172839507e-01, 6.543209876543210e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.481481481481482e-01;
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.330000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.320000000000000e+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP11</strong>.<br><br>
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
end SP11;
