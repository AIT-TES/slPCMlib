within slPCMlib.Media_Rubitherm_SP;
package SP_minus_24 "Rubitherm SP-24; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP_minus_24";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.491500000000000e+02, 2.561500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.451500000000000e+02, 2.511500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.500000000000000e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 2.491500000000000e+02
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
    constant Real[8] data_x =   {-2.400000000000000e+01, -2.300000000000000e+01, -2.200000000000000e+01, -2.100000000000000e+01, -2.000000000000000e+01, -1.900000000000000e+01, -1.800000000000000e+01, -1.700000000000000e+01};
    constant Real[8] data_y =   {0.000000000000000e+00, 5.662100456621004e-01, 3.059360730593607e-01, 7.305936073059360e-02, 3.196347031963470e-02, 1.826484018264840e-02, 4.566210045662100e-03, 0.000000000000000e+00};
    constant Real[8] m_k =      {5.662100456621004e-01, 0.000000000000000e+00, -2.465753424657534e-01, -1.208935079618943e-01, -2.417870159237885e-02, -1.369863013698630e-02, -9.132420091324200e-03, -4.566210045662100e-03};
    constant Real[8] iy_start = {0.000000000000000e+00, 3.152924082818743e-01, 7.511805303305485e-01, 9.220761635465511e-01, 9.645096751357541e-01, 9.876498365419543e-01, 9.981837994914639e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.545949872865965e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    7;
    constant Real[7] data_x =   {-2.800000000000000e+01, -2.700000000000000e+01, -2.600000000000000e+01, -2.500000000000000e+01, -2.400000000000000e+01, -2.300000000000000e+01, -2.200000000000000e+01};
    constant Real[7] data_y =   {0.000000000000000e+00, 4.469274935036961e-03, 3.126495012681377e-02, 1.786112913250865e-01, 7.811651445876423e-01, 4.489339025418144e-03, 0.000000000000000e+00};
    constant Real[7] m_k =      {3.685588235567026e-03, 1.177338496003903e-02, 7.952019421148038e-02, 3.749500972304143e-01, 0.000000000000000e+00, -1.346712752757611e-02, -1.547905755619573e-04};
    constant Real[7] iy_start = {0.000000000000000e+00, 1.560155108260831e-03, 1.377779017266417e-02, 9.407105593399451e-02, 6.050415886256562e-01, 9.988650551181484e-01, 9.999999999999999e-01};
    constant Real    iy_scaler = 9.996800708198453e-01;
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.200000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.300000000000000e+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP-24</strong>.<br><br>
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
end SP_minus_24;
