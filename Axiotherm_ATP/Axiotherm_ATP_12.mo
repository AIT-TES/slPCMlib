
// within slPCMlib.Axiotherm_ATP;
package Axiotherm_ATP_12 "Axiotherm GmbH, ATP 12; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 12";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7714999999999998E+02, 2.9014999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7714999999999998E+02, 2.8814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8882806734418208E+05
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
    constant Real[12] data_x =   {4.0000000000000000E+00, 6.3750000000000000E+00, 9.1250000000000000E+00, 1.1625000000000000E+01, 1.2125000000000000E+01, 1.2375000000000000E+01, 1.3125000000000000E+01, 1.4125000000000000E+01, 1.4625000000000000E+01, 1.4875000000000000E+01, 1.6375000000000000E+01, 1.7000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.4354897162000000E-02, 2.2099131748000000E-02, 7.3134060853000005E-02, 1.1839401506100000E-01, 1.6934053657299999E-01, 3.8760568914799998E-01, 2.7819409373499998E-01, 1.0060142580200000E-01, 5.5880609892999998E-02, 4.2372237339999997E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -9.1682476509999996E-03, 1.0426521169000001E-02, 2.2877247844000002E-02, 1.4307375381000001E-01, 3.0328988962999998E-01, 1.7782996338000001E-01, -3.5650469291499998E-01, -3.2080424032299998E-01, -1.1028686809900000E-01, -1.1916963538000000E-02, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 2.1695281883000000E-02, 6.0070933144999998E-02, 1.7441588520300000E-01, 2.2051473984400000E-01, 2.5620525564300001E-01, 4.7435258089799998E-01, 8.5777664233999995E-01, 9.5322444941200002E-01, 9.7198159802100004E-01, 9.9904891578300004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0158872300920994E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {4.0000000000000000E+00, 6.3750000000000000E+00, 9.6250000000000000E+00, 1.0875000000000000E+01, 1.1375000000000000E+01, 1.1625000000000000E+01, 1.1875000000000000E+01, 1.2125000000000000E+01, 1.2375000000000000E+01, 1.2625000000000000E+01, 1.3375000000000000E+01, 1.4875000000000000E+01, 1.5000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.6595713945999999E-02, 4.7217743219999997E-02, 1.8685541050400001E-01, 4.8051055625900002E-01, 8.0009508415999997E-01, 8.0415740025500004E-01, 4.8922300868399998E-01, 1.2391581124900000E-01, 0.0000000000000000E+00, 0.0000000000000000E+00, 9.4763725900000003E-04, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -1.2515136810000000E-03, 1.4034423426000000E-02, 2.9518990442999998E-01, 8.4965578704199995E-01, 9.0558793793900005E-01, -1.0840723240650001E+00, -1.3744349297920000E+00, -1.1821615360100000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00, -1.0605200972000000E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.0354789966999999E-02, 1.1085964241900000E-01, 2.2086599608700000E-01, 3.7660833537700000E-01, 5.3685803506300001E-01, 7.4836655471299995E-01, 9.1202662147699998E-01, 9.8788782389200003E-01, 9.9724738579100003E-01, 9.9724738579100003E-01, 9.9995444926400001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0029120986997861E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8900000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.0000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.0000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 12</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATP_12;