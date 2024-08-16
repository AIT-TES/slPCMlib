within slPCMlib.Media_Rubitherm_PX;
package Rubitherm_PX82 "Rubitherm GmbH, PX82; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "PX82";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.4214999999999998E+02, 3.5714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.4214999999999998E+02, 3.5714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 7.3000000000000000E+04
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
    constant Real[10] data_x =   {6.9000000000000000E+01, 7.0625000000000000E+01, 7.2625000000000000E+01, 7.4375000000000000E+01, 7.5625000000000000E+01, 7.8625000000000000E+01, 8.0375000000000000E+01, 8.1375000000000000E+01, 8.3875000000000000E+01, 8.4000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 3.7549588245000001E-02, 6.8724074546999997E-02, 5.9287683163000003E-02, 8.5076705185000004E-02, 1.0973061078900000E-01, 9.9469650061999995E-02, 1.0228273193100000E-01, 1.0585641150000001E-03, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 8.4851979609999999E-03, 2.3040846180000001E-02, 4.0983925334999997E-02, -1.0411324921000000E-02, 2.2589425454999999E-02, 3.0569366244000001E-02, -5.5590431705999997E-02, -1.1703650306000001E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 2.8718975720000001E-02, 1.3041384299799999E-01, 2.3813417543199999E-01, 3.3531497335600002E-01, 6.0349554564499996E-01, 7.8499662204700005E-01, 8.9334374713200004E-01, 9.9994894176100002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0026925895579653E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {6.9000000000000000E+01, 7.0375000000000000E+01, 7.1625000000000000E+01, 7.3875000000000000E+01, 7.5375000000000000E+01, 7.6875000000000000E+01, 7.8375000000000000E+01, 7.9375000000000000E+01, 8.2125000000000000E+01, 8.3125000000000000E+01, 8.4000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.8674353925000000E-02, 4.5056810603000000E-02, 3.6320143983000000E-02, 8.4946310923999999E-02, 7.8320114765999996E-02, 1.0586765471300000E-01, 8.4005998329999998E-02, 1.0047954547800000E-01, 5.2089314812999998E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -9.2792480920000000E-03, 4.2770808799999997E-02, 1.4280105376000000E-02, -3.0918828318000000E-02, 2.7788776054000000E-02, -4.3294730586000002E-02, 2.5891440585000002E-02, -2.4672051391999999E-02, -8.0569445197000003E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 2.1182358144000001E-02, 6.0499560004000003E-02, 1.6410128532900001E-01, 2.6355773666400001E-01, 3.7503553509799997E-01, 5.2655298094299996E-01, 6.1575282281999999E-01, 9.0137733352000005E-01, 9.8234577671000001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0003199416796253E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 6.9000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 6.9000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 1.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 1.0000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>PX82</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: other<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-03.<br><br>
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
end Rubitherm_PX82;
