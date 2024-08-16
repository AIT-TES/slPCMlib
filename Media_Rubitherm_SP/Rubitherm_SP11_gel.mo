within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP11_gel "Rubitherm GmbH, SP11_gel; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP11_gel";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7314999999999998E+02, 2.9214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7314999999999998E+02, 2.9214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4300000000000000E+05
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
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {0.0000000000000000E+00, 4.6250000000000000E+00, 6.6250000000000000E+00, 8.6250000000000000E+00, 1.0625000000000000E+01, 1.1375000000000000E+01, 1.1625000000000000E+01, 1.2375000000000000E+01, 1.2625000000000000E+01, 1.3375000000000000E+01, 1.5625000000000000E+01, 1.7375000000000000E+01, 1.9000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.7183263321000000E-02, 3.3328414977000002E-02, 4.8211600067000003E-02, 1.0132230837000000E-01, 2.0414479894199999E-01, 2.6950254830800002E-01, 2.0292062058099999E-01, 1.1847241939500000E-01, 3.8426103443000001E-02, 4.0479458566000000E-02, 4.2488598434999997E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -3.0791180170000000E-03, -4.9752950600000000E-03, -6.3961864510000003E-03, 7.0585731446000000E-02, 2.0516557331900001E-01, 2.1619843335200001E-01, -3.0134839279999998E-01, -2.6528614245900001E-01, -1.6182464112999999E-02, -5.0401009389999998E-03, -3.8023071412000002E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 6.9562201901000001E-02, 1.3179035046199999E-01, 2.1525854427800001E-01, 3.4132876068099999E-01, 4.5149021092600000E-01, 5.1168769056399999E-01, 7.1667863610899996E-01, 7.5737411671999999E-01, 8.0537073508400003E-01, 8.9092979236699998E-01, 9.7338120058499999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0177354426890923E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    17;
    constant Real[17] data_x =   {0.0000000000000000E+00, 2.6250000000000000E+00, 5.6250000000000000E+00, 8.1250000000000000E+00, 9.3750000000000000E+00, 9.6250000000000000E+00, 1.0125000000000000E+01, 1.0375000000000000E+01, 1.0625000000000000E+01, 1.0875000000000000E+01, 1.1375000000000000E+01, 1.2375000000000000E+01, 1.3375000000000000E+01, 1.4625000000000000E+01, 1.6625000000000000E+01, 1.8375000000000000E+01, 1.9000000000000000E+01};
    constant Real[17] data_y =   {0.0000000000000000E+00, 2.3722083807000000E-02, 3.7514800468999999E-02, 7.2117690294000003E-02, 1.9036873602500001E-01, 2.9237515384000001E-01, 2.9821083097500001E-01, 1.9923680351600001E-01, 8.5762613335999993E-02, 4.2957218754000000E-02, 2.9872526859000000E-02, 5.2958611122000003E-02, 2.9379972721999999E-02, 3.3721724927999998E-02, 5.1162490985000003E-02, 5.5850131799999999E-03, 0.0000000000000000E+00};
    constant Real[17] m_k =      {0.0000000000000000E+00, -5.6569723770000004E-03, -4.2050412460000004E-03, 1.6327783857000001E-02, 2.8956698038000001E-01, 3.2554185249800000E-01, -3.5774465935700001E-01, -4.2300328732999998E-01, -3.5913908375300002E-01, -7.6354087693999997E-02, 5.0863244455999997E-02, -4.7808118417000001E-02, 1.9409861295999999E-02, -1.5277521202000001E-02, 1.9708381000999999E-02, -1.6154938657999999E-02, 0.0000000000000000E+00};
    constant Real[17] iy_start = {0.0000000000000000E+00, 3.4646664495999999E-02, 1.2610755198699999E-01, 2.5342074002900000E-01, 3.8287978501900000E-01, 4.4349568977800002E-01, 6.0661597764700004E-01, 6.6961520609699998E-01, 7.0517755025999995E-01, 7.1990653548899997E-01, 7.3558264756400005E-01, 7.8560064241399996E-01, 8.2144058710599999E-01, 8.6573206228599997E-01, 9.3951457929000004E-01, 9.9877122960700004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0076516154875750E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.3300000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.3200000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP11_gel</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_SP11_gel;
