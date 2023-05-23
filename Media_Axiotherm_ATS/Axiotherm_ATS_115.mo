
// within slPCMlib.Axiotherm_ATS;
package Axiotherm_ATS_115 "Axiotherm GmbH, ATS 115; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS 115";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.8614999999999998E+02, 3.9314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.8414999999999998E+02, 3.9014999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.3030245799302441E+05
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
    constant Real[12] data_x =   {1.1300000000000000E+02, 1.1462500000000000E+02, 1.1512500000000000E+02, 1.1537500000000000E+02, 1.1562500000000000E+02, 1.1587500000000000E+02, 1.1637500000000000E+02, 1.1662500000000000E+02, 1.1687500000000000E+02, 1.1737500000000000E+02, 1.1837500000000000E+02, 1.2000000000000000E+02};
    constant Real[12] data_y =   {0.0000000000000000E+00, 5.9856433420000001E-03, 3.5488691461999997E-02, 1.3236155928500001E-01, 5.1742966562000003E-01, 8.8538106661799998E-01, 6.8832999846999998E-01, 2.7079901791700001E-01, 1.1527036662800000E-01, 4.2495174784999999E-02, 3.7335809850999997E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -5.2737457160000004E-03, 1.1463703369700000E-01, 1.3205386962250001E+00, 1.5062646704600000E+00, 1.4171827043790000E+00, -1.4221554969790000E+00, -1.0484121193139999E+00, -2.7822817724299997E-01, 1.5698105011000000E-02, -4.2145767367000000E-02, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 6.0325569100000004E-03, 1.3914396363999999E-02, 2.8636227370999998E-02, 1.0900902524800000E-01, 2.8507892910400001E-01, 7.3831494528399999E-01, 8.5643028478799998E-01, 9.0074165656399996E-01, 9.3410782455200003E-01, 9.7890842038299997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0014480749870300E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {1.1100000000000000E+02, 1.1337500000000000E+02, 1.1387500000000000E+02, 1.1412500000000000E+02, 1.1437500000000000E+02, 1.1462500000000000E+02, 1.1512500000000000E+02, 1.1537500000000000E+02, 1.1562500000000000E+02, 1.1637500000000000E+02, 1.1700000000000000E+02};
    constant Real[11] data_y =   {0.0000000000000000E+00, 6.4913183021999998E-02, 8.0115822188999999E-02, 1.3040736850900000E-01, 2.4551244376000000E-01, 4.9981279268500001E-01, 7.1989070072500005E-01, 5.7688066170100005E-01, 3.3079693393499998E-01, 4.5775984669000000E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 4.5162157769999996E-03, 1.0319535381300000E-01, 2.5002680049300002E-01, 7.3843876274200004E-01, 8.9567573144500001E-01, -2.9818612692500002E-01, -8.3447294400399996E-01, -7.3294127488299998E-01, -1.4017984431100000E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 7.4789832346999999E-02, 1.0891292238899999E-01, 1.3440504474000001E-01, 1.7874939528799999E-01, 2.7088456691300000E-01, 5.9992708771799996E-01, 7.6443895279899998E-01, 8.7711114852299998E-01, 9.9028046660600000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9770928080877441E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.5960000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.5000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.5000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS 115</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATS_115;