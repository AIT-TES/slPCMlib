within slPCMlib.Media_Rubitherm_SP;
package SP15 "Rubitherm SP15; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP15";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.831500000000000e+02, 2.941500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.801500000000000e+02, 2.901500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.800000000000000e+05
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
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01};
    constant Real[12] data_y =   {0.000000000000000e+00, 1.010101010101010e-02, 4.040404040404041e-02, 3.030303030303030e-02, 6.060606060606061e-02, 3.131313131313131e-01, 3.434343434343434e-01, 8.080808080808081e-02, 4.040404040404041e-02, 3.030303030303030e-02, 5.050505050505050e-02, 0.000000000000000e+00};
    constant Real[12] m_k =      {1.010101010101010e-02, 2.020202020202020e-02, 0.000000000000000e+00, 0.000000000000000e+00, 9.090909090909091e-02, 9.090909090909083e-02, 0.000000000000000e+00, -1.195628998584417e-01, -1.992714997640695e-02, 0.000000000000000e+00, 0.000000000000000e+00, -5.050505050505050e-02};
    constant Real[12] iy_start = {0.000000000000000e+00, 4.187604690117254e-03, 3.098827470686767e-02, 6.616415410385260e-02, 1.038525963149079e-01, 2.897822445561139e-01, 6.239530988274706e-01, 8.449218819815625e-01, 8.969621338757656e-01, 9.304857621440534e-01, 9.706867671691790e-01, 9.999999999999998e-01};
    constant Real    iy_scaler = 9.949748743718592e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01};
    constant Real[11] data_y =   {0.000000000000000e+00, 1.208699079732104e-02, 1.179480010425200e-02, 2.359666362094040e-02, 3.511322910619391e-02, 3.426440255205239e-02, 1.301320218939759e-01, 3.687260287382580e-01, 3.718994433033949e-01, 1.238641988361096e-02, 0.000000000000000e+00};
    constant Real[11] m_k =      {1.208699079732104e-02, 0.000000000000000e+00, 0.000000000000000e+00, 1.165921450097095e-02, 0.000000000000000e+00, 0.000000000000000e+00, 1.672308130931028e-01, 9.520243695410602e-03, 0.000000000000000e+00, -3.707709306186725e-02, -2.469766766245731e-03};
    constant Real[11] iy_start = {0.000000000000000e+00, 7.042201996172946e-03, 1.896862993694837e-02, 3.567249774495908e-02, 6.596230187232965e-02, 1.006090889614518e-01, 1.687886939283505e-01, 4.310421367120408e-01, 8.016985966885477e-01, 9.966947434010793e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.987884066089748e-01;
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.400000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.350000000000000e+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP15</strong>.<br><br>
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
end SP15;
