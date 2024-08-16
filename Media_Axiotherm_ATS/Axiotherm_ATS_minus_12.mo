within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_minus_12 "Axiotherm GmbH, ATS -12; data taken from: Axiotherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -12";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.5614999999999998E+02, 2.6814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.5314999999999998E+02, 2.6414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.1200000000000000E+05
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
    constant Real[11] data_x =   {-1.7000000000000000E+01, -1.5375000000000000E+01, -1.4375000000000000E+01, -1.3625000000000000E+01, -1.2625000000000000E+01, -1.1375000000000000E+01, -1.1125000000000000E+01, -1.0375000000000000E+01, -1.0125000000000000E+01, -8.8750000000000000E+00, -5.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 5.8636669539999997E-03, 6.7113989122000006E-02, 1.6246874891300001E-01, 2.2027087332500001E-01, 3.9992920848300001E-01, 3.7542039765000002E-01, 7.0540360997000007E-02, 3.0467027989000001E-02, 1.3083834969000001E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 4.5978395699999997E-03, 1.5018701146399999E-01, 4.6872097006999999E-02, 1.0733076973100000E-01, 1.5601216946900001E-01, -2.6168933995100002E-01, -3.5880608534399999E-01, -7.2611566867000002E-02, -5.1626871289999999E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 3.7300583460000000E-03, 2.7941017448999999E-02, 1.1833442748600000E-01, 3.0355338484900002E-01, 6.8256295907900000E-01, 7.8106547290399997E-01, 9.5182731929300002E-01, 9.6289615473199996E-01, 9.8122294677099997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9402875788160550E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-2.0000000000000000E+01, -1.8375000000000000E+01, -1.6625000000000000E+01, -1.5875000000000000E+01, -1.5125000000000000E+01, -1.4375000000000000E+01, -1.3375000000000000E+01, -1.2625000000000000E+01, -1.2375000000000000E+01, -1.1625000000000000E+01, -1.0375000000000000E+01, -9.8750000000000000E+00, -9.0000000000000000E+00};
    constant Real[13] data_y =   {0.0000000000000000E+00, 7.6983146980000002E-03, 7.6938462055000004E-02, 1.0896730320700000E-01, 3.8217903790999999E-02, 6.2849129771000004E-02, 1.0894095654700001E-01, 2.1688976352900000E-01, 3.1591525869400000E-01, 3.7050524819600000E-01, 1.9876496260000001E-02, 1.0538020990000000E-03, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -4.9675925420000003E-03, 9.2233810816999995E-02, -7.4721369297999995E-02, -2.6415299273000001E-02, 7.8217320426000006E-02, 2.4213045147000001E-02, 3.0940228866399999E-01, 3.4149323316800001E-01, -3.0103865559400000E-01, -1.0615140651299999E-01, -1.9482759819999999E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 7.3830835120000000E-03, 5.6868739266000001E-02, 1.3477954473900000E-01, 1.8796229783200000E-01, 2.2111527636000000E-01, 3.1194215371200001E-01, 4.2127984246599998E-01, 4.8803043214000003E-01, 7.7692923365699995E-01, 9.9658536822099997E-01, 9.9966165833800003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0047732287527227E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.0640000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.0000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATS -12</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATS_minus_12;
