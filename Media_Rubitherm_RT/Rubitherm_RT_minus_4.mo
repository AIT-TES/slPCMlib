
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT_minus_4 "Rubitherm GmbH, RT-4; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT-4";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6414999999999998E+02, 2.7114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6314999999999998E+02, 2.6914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.5137953669636909E+05
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
    constant Real[11] data_x =   {-9.0000000000000000E+00, -7.6250000000000000E+00, -5.1250000000000000E+00, -4.6250000000000000E+00, -4.3750000000000000E+00, -4.1250000000000000E+00, -3.8750000000000000E+00, -3.3750000000000000E+00, -3.1250000000000000E+00, -2.1250000000000000E+00, -2.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 6.7022212578000007E-02, 2.3131779341100001E-01, 3.8312700340599998E-01, 4.9828669064300002E-01, 5.1867378890500004E-01, 4.3604318913200002E-01, 1.2351638351200001E-01, 5.3076308592000000E-02, 5.0832343599999997E-04, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.0197200094999999E-02, 2.1927449103600000E-01, 3.6771171922899998E-01, 3.7848260107400000E-01, -1.5606337484999999E-01, -5.8857292509799997E-01, -5.8010472813699998E-01, -1.4935616728199999E-01, -5.2256756129999997E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 4.4695191059000001E-02, 3.1005574984200002E-01, 4.6133269918999997E-01, 5.7200801165600001E-01, 7.0256651783599999E-01, 8.2477127805799999E-01, 9.6518851451700005E-01, 9.8511900957300003E-01, 9.9997490829199998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0050371954511657E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {-1.0000000000000000E+01, -8.6250000000000000E+00, -6.3750000000000000E+00, -5.6250000000000000E+00, -5.3750000000000000E+00, -5.1250000000000000E+00, -4.3750000000000000E+00, -4.1250000000000000E+00, -4.0000000000000000E+00};
    constant Real[9] data_y =   {0.0000000000000000E+00, 5.4658891773999997E-02, 1.9032031609500000E-01, 3.8103139262100000E-01, 4.9544594176500001E-01, 5.2194000772600002E-01, 1.4911120747500001E-01, 3.3570987658999997E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 2.2299091564000000E-02, 1.4300144347900001E-01, 3.6124657273700000E-01, 3.7467652384700001E-01, -1.1116790874099999E-01, -5.9143485335799995E-01, -3.8608021979399998E-01, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 3.4082340582999998E-02, 2.5887886617500000E-01, 4.6301105434899997E-01, 5.7255740975600000E-01, 7.0232818563900001E-01, 9.7662671872399998E-01, 9.9840369658999994E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0005172747776336E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.6000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.6000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT-4</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-09-29.<br><br>
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
end Rubitherm_RT_minus_4;