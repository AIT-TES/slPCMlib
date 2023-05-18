
within slPCMlib.Media_Axiotherm_ATP;
package Axiotherm_ATP_60 "Axiotherm GmbH, ATP 60; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 60";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.2414999999999998E+02, 3.3514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.2314999999999998E+02, 3.3414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.0078484530170853E+05
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
    constant Real[11] data_x =   {5.1000000000000000E+01, 5.4125000000000000E+01, 5.6625000000000000E+01, 5.7625000000000000E+01, 5.8125000000000000E+01, 5.8375000000000000E+01, 5.8625000000000000E+01, 5.9375000000000000E+01, 6.0375000000000000E+01, 6.1375000000000000E+01, 6.2000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.1562659583000000E-02, 4.5996095106000001E-02, 6.5509491239999995E-02, 1.1991673123300001E-01, 1.9462145328100000E-01, 3.3723080101000003E-01, 5.0669650034699998E-01, 1.2016004515600000E-01, 2.8022608700000000E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 3.4327292039999999E-03, 2.9969612760000000E-02, 1.3898317876000000E-02, 1.8635210354000001E-01, 4.6924819825399999E-01, 5.2316266956799995E-01, -2.9776712606099998E-01, -3.5765075135899999E-01, -9.4102192799999998E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 3.0536886256000002E-02, 1.0033838539600000E-01, 1.5676303125300001E-01, 1.9902687951900000E-01, 2.3642833039500000E-01, 3.0185515303200000E-01, 6.5265946792999996E-01, 9.6735564254499995E-01, 9.9943727197800003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.8830972644438819E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {5.0000000000000000E+01, 5.3625000000000000E+01, 5.6125000000000000E+01, 5.7625000000000000E+01, 5.8125000000000000E+01, 5.8375000000000000E+01, 5.8625000000000000E+01, 5.9125000000000000E+01, 5.9375000000000000E+01, 5.9625000000000000E+01, 5.9875000000000000E+01, 6.0875000000000000E+01, 6.1000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 8.4364965489999999E-03, 4.2307196228999998E-02, 4.1771819373999999E-02, 1.0282625108900000E-01, 2.2661225521100001E-01, 5.4236362332200005E-01, 7.6574822808599996E-01, 5.3102107452000002E-01, 2.1594760818400000E-01, 8.6976394550000002E-02, 6.9236726199999998E-04, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -2.8657608130000000E-03, 1.5228675656000000E-02, -1.8628178963000001E-02, 2.3315572907099999E-01, 8.5382427717700005E-01, 1.1047827520500000E+00, -7.0307667166400001E-01, -1.1165279756549999E+00, -8.7821953875199998E-01, -2.5489130223099998E-01, -7.0594851040000003E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.8721090247000000E-02, 7.3581564158999996E-02, 1.4408788248400001E-01, 1.7548119140700000E-01, 2.1402915633500000E-01, 3.1034522742400000E-01, 6.8081092735100002E-01, 8.4766097144999997E-01, 9.4124951806199997E-01, 9.7641662415800001E-01, 9.9996537949599995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0158325486101671E+00;
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
    lambda := 2.0000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 60</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATP_60;