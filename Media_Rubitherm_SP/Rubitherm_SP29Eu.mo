
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP29Eu "Rubitherm GmbH, SP29Eu; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP29Eu";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9514999999999998E+02, 3.0914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9314999999999998E+02, 3.0614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7044676647367701E+05
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
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {2.2000000000000000E+01, 2.4625000000000000E+01, 2.5625000000000000E+01, 2.6625000000000000E+01, 2.7375000000000000E+01, 2.7625000000000000E+01, 2.7875000000000000E+01, 2.8375000000000000E+01, 2.8625000000000000E+01, 2.8875000000000000E+01, 3.0125000000000000E+01, 3.2625000000000000E+01, 3.5375000000000000E+01, 3.6000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 1.0296367912999999E-02, 5.7531814664000001E-02, 2.0339131052100001E-01, 3.6690736227600002E-01, 4.4602531490000003E-01, 4.5389330396999999E-01, 2.6995884283499999E-01, 1.3922927843100000E-01, 7.4674454176999999E-02, 1.6437341974000000E-02, 2.0013965057000000E-02, 9.8229417059999997E-03, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, -3.6054595150000002E-03, 1.2913313845100000E-01, 1.5136017278300001E-01, 2.6872238216700001E-01, 2.7146495032300000E-01, -1.3940634458000001E-01, -4.9425494339199999E-01, -4.4889597391600000E-01, -1.5295973966900001E-01, -5.1967802949999999E-03, -9.2610637570000008E-03, -2.1366020529000000E-02, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 1.5701787524000001E-02, 3.8726602712999998E-02, 1.6830543286600000E-01, 3.7823680695900003E-01, 4.8060503566000001E-01, 5.9609895352499997E-01, 7.8587458767399998E-01, 8.3717066000999996E-01, 8.6255723667499995E-01, 9.0054637876700006E-01, 9.4858676936700004E-01, 9.9760794160199995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0075384933182903E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    18;
    constant Real[18] data_x =   {2.0000000000000000E+01, 2.1625000000000000E+01, 2.3375000000000000E+01, 2.4125000000000000E+01, 2.4625000000000000E+01, 2.5125000000000000E+01, 2.5375000000000000E+01, 2.5625000000000000E+01, 2.5875000000000000E+01, 2.6125000000000000E+01, 2.6375000000000000E+01, 2.6625000000000000E+01, 2.6875000000000000E+01, 2.7125000000000000E+01, 2.8125000000000000E+01, 2.9375000000000000E+01, 3.1875000000000000E+01, 3.3000000000000000E+01};
    constant Real[18] data_y =   {0.0000000000000000E+00, 1.3884147693000001E-02, 4.2267591983000001E-02, 7.4730916539999998E-02, 3.8769206395999997E-02, 8.4673845846000004E-02, 2.1425341332100001E-01, 6.0755445956800003E-01, 8.5416447983099997E-01, 7.4185254723600003E-01, 3.6736683008400001E-01, 5.1151477066999997E-02, 7.7185135089999999E-03, 1.9136560479999999E-03, 3.7425456598000000E-02, 3.2536328800000000E-02, 7.8283674849999995E-03, 0.0000000000000000E+00};
    constant Real[18] m_k =      {0.0000000000000000E+00, -4.2200617150000000E-03, 6.4156671082999994E-02, -3.6799852394000000E-02, -6.7432955192000002E-02, 2.0319223305900000E-01, 9.1720242241799999E-01, 1.3246794446300001E+00, 4.7134924265299999E-01, -1.3514688203530001E+00, -1.3837132353429999E+00, -6.0186317549699997E-01, -2.6633414740000000E-02, -2.8095335719999998E-03, 4.3774219295000000E-02, -7.7842818250000003E-03, -1.3165335294000000E-02, 0.0000000000000000E+00};
    constant Real[18] iy_start = {0.0000000000000000E+00, 1.2306144006000001E-02, 4.4239377141999998E-02, 9.3230887131000001E-02, 1.2247375288200001E-01, 1.4789613414700001E-01, 1.8180956195600001E-01, 2.8320956814699999E-01, 4.7185026653200002E-01, 6.8250048343299996E-01, 8.2241963586699995E-01, 8.7104413866599995E-01, 8.7544143126200002E-01, 8.7652991766900001E-01, 8.9244245637899999E-01, 9.4328115706899995E-01, 9.9696121101699997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0079151976946745E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.6000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.5000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.0000000000000000E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.5000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP29Eu</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
  <br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-05-15.<br><br>
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
    <li>file creation date: 2023-05-15 </ul>
    </p></html>"));
end Rubitherm_SP29Eu;