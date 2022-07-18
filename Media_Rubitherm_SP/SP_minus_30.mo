within slPCMlib.Media_Rubitherm_SP;
package SP_minus_30 "Rubitherm SP-30; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP_minus_30";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.731500000000000e+02+ERROR
                                                                                      -3.400000000000000e+01, 2.731500000000000e+02+ERROR
                                                                                                                                    -2.600000000000000e+01}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.731500000000000e+02+ERROR
                                                                                            -3.700000000000000e+01, 2.731500000000000e+02+ERROR
                                                                                                                                          -2.700000000000000e+01}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.500000000000000e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+ERROR
                                                                    -3.400000000000000e+01
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
    constant Real[9] data_x =   {-3.400000000000000e+01, -3.300000000000000e+01, -3.200000000000000e+01, -3.100000000000000e+01, -3.000000000000000e+01, -2.900000000000000e+01, -2.800000000000000e+01, -2.700000000000000e+01, -2.600000000000000e+01};
    constant Real[9] data_y =   {0.000000000000000e+00, 4.345578815741169e-03, 4.326694760497350e-03, 2.169020803036863e-02, 6.939731849499568e-02, 6.861679388131205e-01, 1.660622628809743e-01, 4.800999820429320e-02, 0.000000000000000e+00};
    constant Real[9] m_k =      {4.345578815741169e-03, 0.000000000000000e+00, 0.000000000000000e+00, 1.394878665417161e-02, 1.424399763087370e-01, 0.000000000000000e+00, -3.190789703044137e-01, -8.303113144048717e-02, -4.800999820429320e-02};
    constant Real[9] iy_start = {0.000000000000000e+00, 2.523909248760284e-03, 6.841209806722829e-03, 1.863580292029350e-02, 5.332063833078210e-02, 4.412806093320346e-01, 8.920190683621503e-01, 9.790050289860321e-01, 9.999999999999998e-01};
    constant Real    iy_scaler = 9.956559880194912e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-3.700000000000000e+01, -3.600000000000000e+01, -3.500000000000000e+01, -3.400000000000000e+01, -3.300000000000000e+01, -3.200000000000000e+01, -3.100000000000000e+01, -3.000000000000000e+01, -2.900000000000000e+01, -2.800000000000000e+01, -2.700000000000000e+01};
    constant Real[11] data_y =   {0.000000000000000e+00, 8.940404834332715e-03, 2.231104666653091e-02, 2.668151418895903e-02, 7.573446900856781e-02, 1.960913858656345e-01, 4.276653433712030e-01, 2.290646863642071e-01, 4.530599991216530e-03, 8.980549709350936e-03, 0.000000000000000e+00};
    constant Real[11] m_k =      {8.940404834332715e-03, 1.115552333326545e-02, 4.132205318356469e-03, 1.244322130673314e-02, 8.470493583833776e-02, 1.759654371813176e-01, 0.000000000000000e+00, -2.115673716899933e-01, 0.000000000000000e+00, 0.000000000000000e+00, -8.980549709350936e-03};
    constant Real[11] iy_start = {0.000000000000000e+00, 4.279218568800788e-03, 2.046604720212962e-02, 4.423424722815550e-02, 8.935304828405273e-02, 2.174696030297477e-01, 5.435248190212325e-01, 8.890045043850052e-01, 9.880236567778018e-01, 9.947691578099712e-01, 9.999999999999999e-01};
    constant Real    iy_scaler = 9.985088140776179e-01;
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.100000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.200000000000000e+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP-30</strong>.<br><br>
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
end SP_minus_30;
