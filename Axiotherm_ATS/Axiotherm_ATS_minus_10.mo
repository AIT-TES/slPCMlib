
// within slPCMlib.Axiotherm_ATS;
package Axiotherm_ATS_minus_10 "Axiotherm GmbH, ATS -10; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -10";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6114999999999998E+02, 2.6614999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.5414999999999998E+02, 2.6414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.0213362218512589E+05
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
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-1.2000000000000000E+01, -1.0875000000000000E+01, -1.0625000000000000E+01, -1.0375000000000000E+01, -1.0125000000000000E+01, -9.6250000000000000E+00, -9.3750000000000000E+00, -9.1250000000000000E+00, -8.6250000000000000E+00, -7.6250000000000000E+00, -7.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.1136724809999999E-03, 1.2958214806400001E-01, 6.0868565311199996E-01, 1.0940322627120000E+00, 6.9205124887199998E-01, 1.6304189567499999E-01, 4.0836700155000002E-02, 9.2080324279999997E-03, 1.7538231012000000E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 0.0000000000000000E+00, 1.5549857767630000E+00, 1.9275420780940000E+00, 1.9195403621520000E+00, -1.9710162713990000E+00, -1.3689634076990000E+00, -1.3053986575099999E-01, 2.8512594533999999E-02, -3.1679009442000000E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.7587218980000001E-03, 1.0282102097000000E-02, 1.0100075603800000E-01, 3.1476716634200003E-01, 8.4453455006500000E-01, 9.4871681683200004E-01, 9.6783065017100001E-01, 9.7706647681500003E-01, 9.9553202309199995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0041572216250580E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    17;
    constant Real[17] data_x =   {-1.9000000000000000E+01, -1.6375000000000000E+01, -1.4375000000000000E+01, -1.3625000000000000E+01, -1.2875000000000000E+01, -1.2375000000000000E+01, -1.2125000000000000E+01, -1.1875000000000000E+01, -1.1625000000000000E+01, -1.1375000000000000E+01, -1.1125000000000000E+01, -1.0875000000000000E+01, -1.0625000000000000E+01, -1.0375000000000000E+01, -1.0125000000000000E+01, -9.3750000000000000E+00, -9.0000000000000000E+00};
    constant Real[17] data_y =   {0.0000000000000000E+00, 7.9037835930000008E-03, 6.1333239990000003E-03, 1.8810918619000000E-02, 5.4504294607999999E-02, 1.0882394708000001E-02, 1.0613575836000000E-02, 3.2311089278999998E-02, 1.4958277068600001E-01, 7.5908216222299996E-01, 1.2378312016130000E+00, 1.0099818406340000E+00, 3.9649462603800001E-01, 1.2310629791000000E-02, 1.0676225920000001E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[17] m_k =      {0.0000000000000000E+00, -1.1122313058000001E-02, -8.5512173179999994E-03, 6.1562144572000002E-02, -7.0971273325000001E-02, -8.1757479566000005E-02, 1.6076508089000002E-02, 1.0312173995600001E-01, 1.3005344860570001E+00, 2.2232576772780002E+00, 1.1993485990510000E+00, -2.0760116941220002E+00, -1.9876664050480000E+00, -1.0889649222400000E-01, -4.0156548319999997E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[17] iy_start = {0.0000000000000000E+00, 1.6935845921999999E-02, 3.0253923536000000E-02, 3.6384980700000000E-02, 7.0443605128999998E-02, 8.7188500714000000E-02, 8.9388743936000006E-02, 9.4352399901000000E-02, 1.1102537304399999E-01, 2.2094159168499999E-01, 4.7855804391500001E-01, 7.7971442515099998E-01, 9.5689984731500000E-01, 9.9864783704799998E-01, 9.9978565429800004E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0104704964216680E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.2890000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1600000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.1600000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS -10</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    macroencapsulation<br>  Data taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-04-20.<br><br>
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
    <li>file creation date: 2023-04-20 </ul>
    </p></html>"));
end Axiotherm_ATS_minus_10;