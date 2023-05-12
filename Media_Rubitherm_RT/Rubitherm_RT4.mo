
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT4 "Rubitherm GmbH, RT4; data taken from: Rubitherm datasheet; last access: 2020-09-30."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT4";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7214999999999998E+02, 2.8014999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7014999999999998E+02, 2.7814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4266713184559127E+05
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
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {-1.0000000000000000E+00, 1.3750000000000000E+00, 2.6250000000000000E+00, 3.6250000000000000E+00, 3.8750000000000000E+00, 4.6250000000000000E+00, 4.8750000000000000E+00, 6.1250000000000000E+00, 7.0000000000000000E+00};
    constant Real[9] data_y =   {0.0000000000000000E+00, 8.6501082485999997E-02, 2.4955492585200001E-01, 4.5908906808700001E-01, 4.4159320288400000E-01, 1.1634679336900000E-01, 5.8016228581000001E-02, 4.6085590690000001E-03, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 7.0066406558000005E-02, 2.0560288021699999E-01, 2.0318809134700000E-01, -2.3392564016299999E-01, -4.7427583832699999E-01, -1.3316646494100001E-01, -1.0254120919000000E-02, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 7.0025480052999994E-02, 2.6307499089399999E-01, 6.1881901961499997E-01, 7.3407646406799998E-01, 9.5532964812800003E-01, 9.7541734909500000E-01, 9.9863329960699998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0034434960680041E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-3.0000000000000000E+00, -1.1250000000000000E+00, 1.1250000000000000E+00, 2.3750000000000000E+00, 3.1250000000000000E+00, 3.3750000000000000E+00, 3.6250000000000000E+00, 3.8750000000000000E+00, 4.6250000000000000E+00, 4.8750000000000000E+00, 5.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.9227460447000001E-02, 8.0443734727000005E-02, 2.2517526822799999E-01, 3.0501122051300000E-01, 3.5317832410500000E-01, 4.1641008541199998E-01, 4.1134756577199999E-01, 1.0702967868200000E-01, 2.3848316865000000E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 5.5577031819999999E-03, 7.8354098817000006E-02, 8.7892913857999994E-02, 1.6342770594200001E-01, 2.0176129842400001E-01, 2.0253121319199999E-01, -1.8437807092700001E-01, -4.5873054623100001E-01, -2.7238269603699999E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.6556520798000001E-02, 9.8765173607999995E-02, 2.9037525384800000E-01, 4.8754815608099999E-01, 5.7041808360599999E-01, 6.6754544276899996E-01, 7.7405320972600000E-01, 9.8331471075700005E-01, 9.9885313055000002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0096972025007047E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.7000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.7000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT4</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-09-30.<br><br>
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
end Rubitherm_RT4;