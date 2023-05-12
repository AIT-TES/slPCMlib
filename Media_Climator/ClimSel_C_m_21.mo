
// within slPCMlib.Media_Climator;
package ClimSel_C_m_21 "Climator Sweden AB, ClimSel C-21; data taken from: Climator Sweden AB datasheet; last access: 2022-10-14."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C-21";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.4814999999999998E+02, 2.5514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4714999999999998E+02, 2.5314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.1000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {4.2000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.1717441478481572E+05
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
    constant Real[12] data_x =   {-2.5000000000000000E+01, -2.3375000000000000E+01, -2.1875000000000000E+01, -2.1625000000000000E+01, -2.1375000000000000E+01, -2.1125000000000000E+01, -2.0875000000000000E+01, -2.0625000000000000E+01, -2.0375000000000000E+01, -2.0125000000000000E+01, -1.9125000000000000E+01, -1.8000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.3902826918999999E-02, 2.5208319073399998E-01, 4.1240451083899998E-01, 6.6706373588800005E-01, 7.7937863748100000E-01, 6.9086417015500001E-01, 4.4761342705399998E-01, 1.7614280548799999E-01, 7.3086718622000002E-02, 1.6255542416999998E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 6.4465426880000003E-03, 4.7110914995100001E-01, 7.5249831732000005E-01, 8.3665920219800005E-01, 6.4355028652999996E-02, -8.9674620323700005E-01, -1.0270537592810001E+00, -8.8056786663999997E-01, -1.8839144141100000E-01, -1.7028952375999998E-02, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 9.9029286979999999E-03, 1.2255779945500000E-01, 2.0436348746200000E-01, 3.3920531165800000E-01, 5.2450937802099995E-01, 7.1378201849900003E-01, 8.5713892393900004E-01, 9.3454448835600001E-01, 9.6216409369400002E-01, 9.9263334260500002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0025772764804808E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {-2.6000000000000000E+01, -2.4625000000000000E+01, -2.3875000000000000E+01, -2.3625000000000000E+01, -2.3375000000000000E+01, -2.3125000000000000E+01, -2.2625000000000000E+01, -2.2375000000000000E+01, -2.1625000000000000E+01, -2.0875000000000000E+01, -2.0375000000000000E+01, -2.0000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 5.9959160960000001E-03, 1.3931952106000000E-02, 7.3729877119000004E-02, 4.1171652678600001E-01, 8.1450764771399997E-01, 8.2991841849000003E-01, 4.3223631131599999E-01, 4.2436761859000000E-02, 1.8300426229999999E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -1.3368645238000000E-02, 4.8254479527000002E-02, 7.2012946944800005E-01, 1.4966607759170001E+00, 1.4735927963030000E+00, -1.2735378310059999E+00, -1.0769370920460000E+00, -1.5289399433299999E-01, -8.0429455099999998E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 6.2398222710000001E-03, 1.0832559114999999E-02, 1.8304556542999999E-02, 7.5044331116999999E-02, 2.2872257654200001E-01, 6.9791608413799999E-01, 8.5494765248100002E-01, 9.8988145093400004E-01, 9.9970952131299995E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0018258228825794E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.3000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.3000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 1.4500000000000000E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.3000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C-21</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    macroencapsulation<br>  Data taken from: Climator Sweden AB datasheet - last access 2022-10-14.<br><br>
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
end ClimSel_C_m_21;