
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP31 "Rubitherm GmbH, SP31; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP31";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9714999999999998E+02, 3.1214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9714999999999998E+02, 3.1214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9537681470612742E+05
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
    constant Real[10] data_x =   {2.4000000000000000E+01, 2.9125000000000000E+01, 3.0625000000000000E+01, 3.1375000000000000E+01, 3.1625000000000000E+01, 3.2375000000000000E+01, 3.2625000000000000E+01, 3.4375000000000000E+01, 3.7625000000000000E+01, 3.9000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 3.2245858091000001E-02, 8.5691430314999997E-02, 2.2564735554699999E-01, 3.4029234471100001E-01, 3.8997527292299999E-01, 2.9153964077200001E-01, 3.2676686050999999E-02, 2.3786071573999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.6462927483000000E-02, 6.0604147051999999E-02, 3.5933274542799998E-01, 3.9749991741899998E-01, -3.4210911753099998E-01, -3.2773443130600000E-01, -1.7211910779999999E-02, 5.5268126710000002E-03, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 4.7204794533999998E-02, 1.2842895563400000E-01, 2.3252072967300000E-01, 3.0398620547799998E-01, 6.1653718839700000E-01, 7.0276388332100004E-01, 9.0987668811400002E-01, 9.8255125148199995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0130670910538264E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    16;
    constant Real[16] data_x =   {2.4000000000000000E+01, 2.6375000000000000E+01, 2.7125000000000000E+01, 2.8375000000000000E+01, 2.8625000000000000E+01, 2.9375000000000000E+01, 2.9625000000000000E+01, 3.0125000000000000E+01, 3.0375000000000000E+01, 3.0625000000000000E+01, 3.0875000000000000E+01, 3.1875000000000000E+01, 3.3125000000000000E+01, 3.6125000000000000E+01, 3.8625000000000000E+01, 3.9000000000000000E+01};
    constant Real[16] data_y =   {0.0000000000000000E+00, 6.9055961840000001E-03, 1.7026738598000000E-02, 1.2059243362900000E-01, 1.1855187482500000E-01, 3.3269003074699999E-01, 5.8058637378199995E-01, 5.6414693583200004E-01, 3.1020343468500000E-01, 7.6577052310999996E-02, 2.1062518193999999E-02, 2.8897604465999999E-02, 1.6678284198999999E-02, 1.5271567985000000E-02, 3.8144976929999999E-03, 0.0000000000000000E+00};
    constant Real[16] m_k =      {0.0000000000000000E+00, 1.4976274070000000E-03, 3.1984099630999997E-02, 3.0345020109999999E-02, 2.0209375075000001E-02, 6.5443706550100000E-01, 7.7366528022299996E-01, -9.6688220182700002E-01, -9.8004494281300003E-01, -7.0548919875700000E-01, -6.1502696670000002E-02, 3.9216072764000000E-02, -7.3486496479999996E-03, 2.0324949229999999E-03, -1.7084404130000001E-02, 0.0000000000000000E+00};
    constant Real[16] iy_start = {0.0000000000000000E+00, 7.5102008530000004E-03, 1.5069631530000000E-02, 1.0145340123400000E-01, 1.3145422939199999E-01, 2.7119670745500002E-01, 3.8494380750800000E-01, 7.0798075651400005E-01, 8.1754396832400000E-01, 8.6454772219499998E-01, 8.7341482799100001E-01, 8.9003212261400000E-01, 9.2464362033799996E-01, 9.6560763851599996E-01, 9.9948404365300003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0018366425772991E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.3500000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.3000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP31</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_SP31;