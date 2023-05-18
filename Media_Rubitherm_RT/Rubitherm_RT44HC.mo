
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT44HC "Rubitherm GmbH, RT44HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT44HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1014999999999998E+02, 3.2014999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.0914999999999998E+02, 3.1814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.2067121693222999E+05
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
    constant Real[10] data_x =   {3.7000000000000000E+01, 3.8625000000000000E+01, 4.0625000000000000E+01, 4.1125000000000000E+01, 4.2125000000000000E+01, 4.3375000000000000E+01, 4.3875000000000000E+01, 4.5375000000000000E+01, 4.6625000000000000E+01, 4.7000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 5.1433333740000002E-03, 2.4473561275999998E-02, 5.7800956095000001E-02, 3.4691875685300000E-01, 3.1480674383399998E-01, 1.7586016332400001E-01, 3.5007127348000003E-02, 8.1181533919999996E-03, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, -7.5669214459999998E-03, 1.9518132938000001E-02, 1.0729226608400000E-01, 2.8895623298500001E-01, -2.8683946957400003E-01, -2.1195750374200001E-01, 1.2916414000000001E-04, -3.5959016333999998E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 5.8841486350000001E-03, 2.6613869636000001E-02, 4.5482373372999998E-02, 2.3398736140900001E-01, 7.2588924960000001E-01, 8.4782637675799999E-01, 9.6702236760600002E-01, 9.9889169294100000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0068571055088800E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {3.6000000000000000E+01, 3.7875000000000000E+01, 3.9375000000000000E+01, 3.9875000000000000E+01, 4.0375000000000000E+01, 4.0625000000000000E+01, 4.1125000000000000E+01, 4.2125000000000000E+01, 4.2625000000000000E+01, 4.2875000000000000E+01, 4.3375000000000000E+01, 4.3625000000000000E+01, 4.3875000000000000E+01, 4.4625000000000000E+01, 4.5000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 1.2142522260000000E-02, 1.6747605916300001E-01, 2.4987438420700001E-01, 1.9708290335500001E-01, 1.3294032977600001E-01, 9.9085142710000004E-02, 2.6785827917299998E-01, 4.0561848717999999E-01, 3.7180069253199999E-01, 1.0678561918399999E-01, 1.4383574317999999E-02, 2.6790338899999998E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, 1.7917526622999998E-02, 1.9239983471899999E-01, 4.0748446085999999E-02, -2.1364151901100001E-01, -1.8960626685800000E-01, 1.7874125216999999E-02, 2.4515922437200000E-01, 2.5658595058000000E-01, -4.2817910115899999E-01, -5.3023445592100005E-01, -9.8485155017000003E-02, -8.7137442119999999E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 6.1477506370000001E-03, 1.0836925145400000E-01, 2.1610128504199999E-01, 3.3339627840699998E-01, 3.7461391574800001E-01, 4.2841514377200002E-01, 5.9330614161399997E-01, 7.6180485852199997E-01, 8.6276899519100003E-01, 9.8480795590500003E-01, 9.9773360829400004E-01, 9.9940251563100002E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0021862822495136E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.0000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.0000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT44HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
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
end Rubitherm_RT44HC;