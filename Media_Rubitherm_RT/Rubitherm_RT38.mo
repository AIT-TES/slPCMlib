
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT38 "Rubitherm GmbH, RT38; data taken from: Rubitherm datasheet; last access: 2022-11-30."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT38";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.0114999999999998E+02, 3.1514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.0114999999999998E+02, 3.1414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.5368058990392770E+05
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
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {2.8000000000000000E+01, 3.2125000000000000E+01, 3.3875000000000000E+01, 3.5625000000000000E+01, 3.6625000000000000E+01, 3.7625000000000000E+01, 3.8625000000000000E+01, 4.0125000000000000E+01, 4.1875000000000000E+01, 4.2000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 3.2561741666000003E-02, 1.0313295107600000E-01, 1.3797391233600001E-01, 1.6568272596299999E-01, 1.4673950111799999E-01, 1.3461921393300000E-01, 3.7474841170000002E-02, 8.6074934599999998E-04, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.8598466620999999E-02, 5.8075885180000002E-02, -1.2303759231000000E-02, 4.6795677740999998E-02, -5.1259418303000001E-02, 2.8237008189999999E-03, -3.2454768490999999E-02, -9.5967657360000007E-03, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 4.1228966367999997E-02, 1.5106548701100000E-01, 3.8247870680800000E-01, 5.3097555877699998E-01, 6.9714101397899997E-01, 8.3479053052499996E-01, 9.7194758597399999E-01, 9.9995825095500002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0108471707817359E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {2.8000000000000000E+01, 3.0875000000000000E+01, 3.3625000000000000E+01, 3.5375000000000000E+01, 3.6375000000000000E+01, 3.7875000000000000E+01, 3.8875000000000000E+01, 3.9875000000000000E+01, 4.0875000000000000E+01, 4.1000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.7121576264000000E-02, 9.1694231252999997E-02, 1.4276373847500001E-01, 1.4537606995399999E-01, 1.6510183772199999E-01, 1.1395677922100000E-01, 7.0539152859999997E-03, 5.0794831000000002E-05, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 7.4749554020000003E-03, 1.5643912688000000E-02, -1.4953075340000000E-02, 4.1982772893000002E-02, -2.5263647159000001E-02, -1.2970492361399999E-01, -2.1105858088000001E-02, -5.1539326499999997E-04, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.4134069097000001E-02, 1.9373969000100000E-01, 4.0855914886599998E-01, 5.4910135886000000E-01, 7.9671256728999995E-01, 9.4624007930099996E-01, 9.9814495048100005E-01, 9.9999747453999999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0087346679547855E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.5000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.5000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT38</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2022-11-30.<br><br>
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
end Rubitherm_RT38;