within slPCMlib.Media_Rubitherm_SP;
package SP25E2 "Rubitherm SP25E2; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP25E2";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.731500000000000e+02+1.800000000000000e+01, 2.731500000000000e+02+2.800000000000000e+01}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.731500000000000e+02+1.600000000000000e+01, 2.731500000000000e+02+2.500000000000000e+01}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.800000000000000e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+1.800000000000000e+01
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
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01};
    constant Real[11] data_y =   {0.000000000000000e+00, 7.793844889164175e-03, 1.552694576017200e-02, 2.319977604247893e-02, 6.978202992152556e-02, 1.471766084959079e-01, 1.849987612837984e-01, 3.628153420844410e-01, 1.651414933557215e-01, 2.356519816678688e-02, 0.000000000000000e+00};
    constant Real[11] m_k =      {7.793844889164175e-03, 7.763472880085999e-03, 6.287614531851921e-03, 2.214309879778846e-02, 6.198841622671449e-02, 5.347162810890326e-02, 1.000770810974119e-01, 0.000000000000000e+00, -1.696250719588271e-01, -6.798125505402933e-02, -1.940144435444008e-02};
    constant Real[11] iy_start = {0.000000000000000e+00, 3.890636197278066e-03, 1.564737572330138e-02, 3.364865044639485e-02, 7.672149543614998e-02, 1.856636541948293e-01, 3.475007840392678e-01, 6.291093901900263e-01, 9.065943735414277e-01, 9.922832068433995e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.977388502814215e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01};
    constant Real[10] data_y =   {0.000000000000000e+00, 6.993006993006993e-03, 1.398601398601399e-02, 2.797202797202797e-02, 3.496503496503497e-02, 4.195804195804196e-02, 1.328671328671329e-01, 4.755244755244755e-01, 2.657342657342657e-01, 0.000000000000000e+00};
    constant Real[10] m_k =      {6.993006993006993e-03, 6.993006993006993e-03, 1.048951048951049e-02, 1.048951048951049e-02, 2.966881599384115e-03, 2.076817119568880e-02, 2.167832167832168e-01, 0.000000000000000e+00, -2.377622377622378e-01, -2.657342657342657e-01};
    constant Real[10] iy_start = {0.000000000000000e+00, 3.418803418803420e-03, 1.339031339031339e-02, 3.390313390313390e-02, 6.528531961896757e-02, 1.014416817772174e-01, 1.709401709401709e-01, 4.860398860398860e-01, 8.678062678062677e-01, 9.999999999999999e-01};
    constant Real    iy_scaler = 9.777777777777777e-01;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP25E2</strong>.<br><br>
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
end SP25E2;
