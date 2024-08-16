within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP21EK "Rubitherm GmbH, SP21EK; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP21EK";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8514999999999998E+02, 2.9814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8414999999999998E+02, 2.9614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.5100000000000000E+05
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
    constant Real[11] data_x =   {1.2000000000000000E+01, 1.4625000000000000E+01, 1.9125000000000000E+01, 2.0125000000000000E+01, 2.1125000000000000E+01, 2.1625000000000000E+01, 2.1875000000000000E+01, 2.2625000000000000E+01, 2.2875000000000000E+01, 2.4125000000000000E+01, 2.5000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 5.7677651760000004E-03, 5.5645492483999998E-02, 9.6785489710999997E-02, 3.6500167964000002E-01, 4.7790106836700003E-01, 4.5418125802900000E-01, 1.1281877695000000E-01, 5.4787377838000001E-02, 4.4488327209999999E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -9.4480404199999998E-04, 1.6276620727999999E-02, 1.1646816527700000E-01, 2.9110414493999998E-01, 1.3701748896400001E-01, -2.6522571512300003E-01, -4.9896734037000001E-01, -1.2872647227600001E-01, -9.7537527029999994E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 8.1172711040000001E-03, 1.1729721432800000E-01, 1.8520151545499999E-01, 4.0166357240700001E-01, 6.1571951843600004E-01, 7.3439142025000004E-01, 9.5809860833500005E-01, 9.7713172051499997E-01, 9.9867520199899995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0005614784728911E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {1.1000000000000000E+01, 1.3875000000000000E+01, 1.5875000000000000E+01, 1.7625000000000000E+01, 1.8875000000000000E+01, 1.9625000000000000E+01, 2.0375000000000000E+01, 2.0625000000000000E+01, 2.0875000000000000E+01, 2.1625000000000000E+01, 2.1875000000000000E+01, 2.2875000000000000E+01, 2.3000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.0003137955999999E-02, 3.2026914855999997E-02, 5.2916673740000003E-02, 6.4501627202000006E-02, 1.6332031714900000E-01, 3.8234667984699999E-01, 4.7429650916999999E-01, 4.8154647942300000E-01, 1.1956958329000000E-01, 5.3456706457999999E-02, 5.6518765399999997E-04, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -3.9522654810000002E-03, 1.4333551033000000E-02, 2.7432449418999999E-02, 3.6997298482000000E-02, 2.4088081459399999E-01, 3.3449506677599999E-01, 3.3531667482999999E-01, -1.7911598370099999E-01, -5.2924636094199995E-01, -1.4619770135099999E-01, -5.8457398799999999E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 3.1392359813000002E-02, 7.7203857545999996E-02, 1.4799604133399999E-01, 2.1994343937499999E-01, 2.9561598697899999E-01, 4.9531553313499999E-01, 6.0210427433699998E-01, 7.2393612408300001E-01, 9.6511796675399997E-01, 9.8469851483399995E-01, 9.9997236178899995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9731613156152854E-01;
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
    lambda := 5.0000000000000000E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP21EK</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP21EK;
