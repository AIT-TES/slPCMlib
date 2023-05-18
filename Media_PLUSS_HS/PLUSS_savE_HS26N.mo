
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS26N "Pluss Advanced Technologies Pvt Ltd, HS26N; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS26N";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.4414999999999998E+02, 2.5414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4314999999999998E+02, 2.4814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.7000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.6000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.4418800275668807E+05
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
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {-2.9000000000000000E+01, -2.7375000000000000E+01, -2.6875000000000000E+01, -2.6625000000000000E+01, -2.6125000000000000E+01, -2.5625000000000000E+01, -2.5375000000000000E+01, -2.4625000000000000E+01, -2.4375000000000000E+01, -2.3625000000000000E+01, -2.3125000000000000E+01, -2.1625000000000000E+01, -2.0375000000000000E+01, -1.9000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 4.5323387879999998E-03, 1.4903250392000000E-02, 4.6884896394000002E-02, 2.8918091355300002E-01, 3.6848768723500003E-01, 2.8103234299899998E-01, 2.2317635373700001E-01, 2.5728984074099998E-01, 1.9676019334000000E-01, 9.9188317807999998E-02, 3.4400344431999998E-02, 3.4875716357000000E-02, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, -8.2435288490000005E-03, 4.2653729028999998E-02, 3.2459948018800000E-01, 4.8717747370600001E-01, -2.4224112799600001E-01, -2.3041947632500001E-01, 8.0691545500999995E-02, 8.3903199516000004E-02, -2.1201577005500000E-01, -1.1232641756700000E-01, -2.5006722634000000E-02, 1.6060649026000001E-02, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 5.5197982049999999E-03, 9.3344154669999999E-03, 1.5615944478000000E-02, 9.6586665995000001E-02, 2.7696035441799999E-01, 3.5843221005499998E-01, 5.3366579279100002E-01, 5.9396150027000005E-01, 7.7888094352799997E-01, 8.5109561155000002E-01, 9.3526947800000004E-01, 9.7338034822200004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0042330777104329E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {-3.0000000000000000E+01, -2.8875000000000000E+01, -2.8625000000000000E+01, -2.8375000000000000E+01, -2.7875000000000000E+01, -2.7625000000000000E+01, -2.7375000000000000E+01, -2.6875000000000000E+01, -2.6625000000000000E+01, -2.6375000000000000E+01, -2.5875000000000000E+01, -2.5625000000000000E+01, -2.5375000000000000E+01, -2.5125000000000000E+01, -2.5000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 4.7799538959000000E-02, 9.1119461274000002E-02, 1.6938738097700001E-01, 2.2805903779700001E-01, 1.8110455461000000E-01, 1.0446980114300000E-01, 1.1107421947100000E-01, 2.0057086401999999E-01, 4.4020065497999999E-01, 6.4266193108699998E-01, 4.9509754236400000E-01, 2.5809992590500003E-01, 6.0437490745000003E-02, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, 1.1401441425000000E-01, 2.3547904178599999E-01, 2.7294079068799998E-01, -1.1040596225200000E-01, -2.4324556199200001E-01, -2.0909181016100001E-01, 1.5441541625899999E-01, 6.7179963101700002E-01, 8.4239816397900003E-01, -3.4432014646600001E-01, -8.8074700508100001E-01, -8.6333500184900003E-01, -7.1523758553100003E-01, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 1.4909210928000001E-02, 3.1694290160999998E-02, 6.4164735671000006E-02, 1.7185168493299999E-01, 2.2385268449199999E-01, 2.5948374665200002E-01, 3.0594291880600000E-01, 3.4231833953899998E-01, 4.2177634625900001E-01, 7.1814814649699998E-01, 8.6361985521100004E-01, 9.5797584477400000E-01, 9.9714496931499996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0031575240964345E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.1220000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.2000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.0000000000000000E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 6.9999999999999996E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS26N</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: multiple options available<br>  The data is taken from: PLUSS datasheet - last access 2022-02-13.<br><br>
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
end PLUSS_savE_HS26N;