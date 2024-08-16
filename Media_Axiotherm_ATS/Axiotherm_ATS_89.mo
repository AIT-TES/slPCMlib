within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_89 "Axiotherm GmbH, ATS 89; data taken from: Axiotherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS 89";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.5714999999999998E+02, 3.6614999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.5614999999999998E+02, 3.6414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.1400000000000000E+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature  Tref = rangeTmelting[1]
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
    constant Real[12] data_x =   {8.4000000000000000E+01, 8.5625000000000000E+01, 8.7875000000000000E+01, 8.8875000000000000E+01, 8.9375000000000000E+01, 8.9625000000000000E+01, 9.0125000000000000E+01, 9.0375000000000000E+01, 9.0625000000000000E+01, 9.0875000000000000E+01, 9.1625000000000000E+01, 9.3000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.7761791886000002E-02, 4.3835038743000002E-02, 1.5396922247100001E-01, 4.0963935141800001E-01, 6.9714621823599998E-01, 6.4999376467299996E-01, 3.4528733011899998E-01, 7.7555628227000004E-02, 1.8497225491999999E-02, 6.7639216580000000E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -5.0333048639999997E-03, 3.6830533034000003E-02, 2.4763934603099999E-01, 7.9035457218299998E-01, 9.0355967067900000E-01, -1.1386084390900000E+00, -1.1555808540970001E+00, -6.9983246703299995E-01, -6.0572303564999999E-02, 9.6138621389999998E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.5677027973999998E-02, 6.7770662874000007E-02, 1.4982762451000001E-01, 2.8057397716299998E-01, 4.1955582536899999E-01, 8.0225434729800005E-01, 9.2785839807799997E-01, 9.7878834308499996E-01, 9.8754252050299995E-01, 9.9378037725299995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0088797505763683E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {8.3000000000000000E+01, 8.4875000000000000E+01, 8.6625000000000000E+01, 8.7375000000000000E+01, 8.7625000000000000E+01, 8.8375000000000000E+01, 8.8625000000000000E+01, 8.9125000000000000E+01, 8.9375000000000000E+01, 8.9625000000000000E+01, 8.9875000000000000E+01, 9.1000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.4592691288000000E-02, 1.2073377367400000E-01, 1.5361009521100000E-01, 1.4106009376699999E-01, 2.9323223768200002E-01, 4.6912429180800003E-01, 5.2434306052299995E-01, 3.7831043767099998E-01, 1.8902275304200000E-01, 9.2596013800999993E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.1331706214000000E-02, 1.0468450372599999E-01, -1.7270447836999998E-02, -2.1068881176000000E-02, 4.9256922601199998E-01, 5.6947099069300000E-01, -4.7185452496600000E-01, -6.7806254967400004E-01, -5.7990138082800002E-01, -2.3541849489000000E-01, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.0370047788999999E-02, 1.0504060588100000E-01, 2.1373299595500000E-01, 2.5061940393200000E-01, 3.8952595141300000E-01, 4.8450457505299999E-01, 7.5480642415999999E-01, 8.6881364581300002E-01, 9.3928179714000004E-01, 9.7271974096000002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0008913887653359E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.6810000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.5800000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 5.9999999999999998E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS 89</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
  <br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-05-18.<br><br>
   See:<br>
    Barz, T., Bres, A., & Emhofer, J. (2022).
    slPCMlib: A Modelica Library for the Prediction of Effective 
    Thermal Material Properties of Solid/Liquid Phase Change  
    Materials (PCM). 
    In Proceedings of Asian Modelica Conference 2022 (pp. 63-74). 
    Linkoping University Electronic Press. 
    <a href>https://doi.org/10.3384/ecp19363</a>.
    </p>
    </blockquote>
    </p></html>",
    revisions="<html>
    <ul>
    <li>file creation date: 2023-05-18 </ul>
    </p></html>"));
end Axiotherm_ATS_89;
