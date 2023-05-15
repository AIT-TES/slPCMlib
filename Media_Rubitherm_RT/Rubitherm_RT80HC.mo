
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT80HC "Rubitherm GmbH, RT80HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT80HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.4414999999999998E+02, 3.5514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.4414999999999998E+02, 3.5314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.0268013694701414E+05
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
    constant Real[12] data_x =   {7.1000000000000000E+01, 7.2875000000000000E+01, 7.4875000000000000E+01, 7.6625000000000000E+01, 7.8375000000000000E+01, 7.8625000000000000E+01, 7.8875000000000000E+01, 7.9125000000000000E+01, 7.9625000000000000E+01, 7.9875000000000000E+01, 8.0875000000000000E+01, 8.2000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 9.9599005480000002E-03, 3.8572737984000002E-02, 8.6888149715000002E-02, 4.1444338726500002E-01, 5.1603427814199998E-01, 5.0553860475000001E-01, 3.8764997142099999E-01, 6.4756738127999994E-02, 2.0607506046999999E-02, 5.2869820470000002E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -1.3452480050000000E-03, 3.0117992784000000E-02, 8.2025463085000000E-02, 3.3682309686900003E-01, 3.3638047404700000E-01, -3.3068451981199998E-01, -6.5275114546400004E-01, -4.2478896541900002E-01, -6.1684115998000000E-02, -2.3420712720000000E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 9.6675635229999998E-03, 4.7462411190999997E-02, 1.4335903408400000E-01, 5.1454196974499999E-01, 6.3008954123500005E-01, 7.6039834873599998E-01, 8.7297953416499996E-01, 9.8061987030599995E-01, 9.8934152880899995E-01, 9.9729101011900001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9342766073907596E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {7.1000000000000000E+01, 7.2375000000000000E+01, 7.4875000000000000E+01, 7.6875000000000000E+01, 7.7375000000000000E+01, 7.7625000000000000E+01, 7.8125000000000000E+01, 7.8375000000000000E+01, 7.9625000000000000E+01, 8.0000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 1.4579498369000000E-02, 4.2481472950000002E-02, 1.3581928523100001E-01, 2.6436165504400000E-01, 3.8064634937699998E-01, 4.8423882330599999E-01, 4.4434922885099998E-01, 4.1402791166000000E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, -1.0298101800000001E-02, 2.6261075902000001E-02, 1.5659661136399999E-01, 3.7488625849700002E-01, 4.0677157187300000E-01, -3.8259117975999998E-02, -3.6739471553500003E-01, -1.7191226489100001E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.1623639894000001E-02, 6.3808711069999996E-02, 1.9840661270700000E-01, 2.9372167292800000E-01, 3.7402786368099999E-01, 5.9908976464300001E-01, 7.1665245370700004E-01, 9.9426255743400005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9808922955132995E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.0000000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT80HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2021-05-19.<br><br>
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
end Rubitherm_RT80HC;