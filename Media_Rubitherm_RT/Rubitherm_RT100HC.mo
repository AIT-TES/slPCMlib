
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT100HC "Rubitherm GmbH, RT100HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT100HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.6914999999999998E+02, 3.7814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.6614999999999998E+02, 3.7414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.3328493594090350E+05
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
    constant Real[10] data_x =   {9.6000000000000000E+01, 9.8125000000000000E+01, 9.9625000000000000E+01, 1.0037500000000000E+02, 1.0062500000000000E+02, 1.0137500000000000E+02, 1.0162500000000000E+02, 1.0237500000000000E+02, 1.0387500000000000E+02, 1.0500000000000000E+02};
    constant Real[10] data_y =   {0.0000000000000000E+00, 3.4369432826000003E-02, 1.3894002259200000E-01, 3.6829759327900002E-01, 4.8481884301700001E-01, 3.9430149019600003E-01, 2.4919043215600001E-01, 5.8548806538999999E-02, 1.6111832406999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 2.3094933229000001E-02, 2.1639964639299999E-01, 3.9558084416599998E-01, 4.0844190391899998E-01, -5.1475187616600004E-01, -4.6793513123899999E-01, -1.1708696410900001E-01, -9.0373770220000000E-03, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 2.7663869601000000E-02, 1.2085230910200000E-01, 3.0160240980199998E-01, 4.0755078632399999E-01, 7.7831128957499995E-01, 8.5803425595100002E-01, 9.5641087818000003E-01, 9.9193775373600002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9414302076202454E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {9.3000000000000000E+01, 9.5625000000000000E+01, 9.6625000000000000E+01, 9.7375000000000000E+01, 9.7625000000000000E+01, 9.8125000000000000E+01, 9.8375000000000000E+01, 9.8625000000000000E+01, 9.8875000000000000E+01, 9.9375000000000000E+01, 9.9625000000000000E+01, 1.0037500000000000E+02, 1.0100000000000000E+02};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.3970246095000000E-02, 4.7139006942000003E-02, 6.5886005634999995E-02, 5.5865675676000003E-02, 1.1164477668000000E-01, 2.2088882986500000E-01, 4.8084597878500002E-01, 6.8331271987300002E-01, 5.8171944635700001E-01, 3.3004907926600002E-01, 4.4555404670999998E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -8.9645877000000006E-05, 6.1260417779000001E-02, -1.6528376402999999E-02, -2.1833841887000001E-02, 2.2008669149499999E-01, 7.5820603969800004E-01, 9.2894782427599998E-01, 6.5497116192399996E-01, -8.4849076393500000E-01, -7.4160606047300004E-01, -1.3686374243100000E-01, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.8441363031999999E-02, 4.3958117573999998E-02, 9.0123876108999995E-02, 1.0541519410500000E-01, 1.4236085754399999E-01, 1.8123856587600001E-01, 2.6832084163499997E-01, 4.1569870263000003E-01, 7.6429847855099997E-01, 8.7804554909800003E-01, 9.9050386096599996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0029334553878693E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.0000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.5000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.5000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT100HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
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
end Rubitherm_RT100HC;