
within slPCMlib.Media_Croda_Crodatherm;
package Croda_Crodatherm_24W "Croda International Plc, Crodatherm 24W; data taken from: Croda datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 24W";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8814999999999998E+02, 2.9914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8514999999999998E+02, 2.9714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.7000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.2000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.6020999999999997E+05
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
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {1.5000000000000000E+01, 1.9375000000000000E+01, 2.2375000000000000E+01, 2.2875000000000000E+01, 2.3375000000000000E+01, 2.3875000000000000E+01, 2.5125000000000000E+01, 2.6000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 7.9140503651000005E-02, 2.3877413441600001E-01, 2.6465101232400001E-01, 1.9881121955000000E-01, 1.0865675589899999E-01, 4.0263085008000001E-02, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 7.6728685570000002E-03, 9.2022048549999999E-02, -4.8479363173000001E-02, -1.9539493133999999E-01, -1.0617894456300001E-01, -6.2574507062999996E-02, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 1.6116660770399999E-01, 5.7551036661199995E-01, 7.0452220989900005E-01, 8.2365946760200004E-01, 8.9880084942399996E-01, 9.8635311941799997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0017738610614109E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {1.2000000000000000E+01, 1.6875000000000000E+01, 1.8875000000000000E+01, 2.0375000000000000E+01, 2.1375000000000000E+01, 2.1875000000000000E+01, 2.2375000000000000E+01, 2.2625000000000000E+01, 2.3125000000000000E+01, 2.3625000000000000E+01, 2.3875000000000000E+01, 2.4000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 3.9272370676999999E-02, 6.4146412052999993E-02, 9.7658255131999994E-02, 1.5673428209300000E-01, 1.7065583362200001E-01, 2.6560971157000002E-01, 3.6499708368299999E-01, 3.5170985404499999E-01, 1.1946146463399999E-01, 2.7039799544000001E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.5403260138999999E-02, 1.5728400122000000E-02, 7.6251917764000005E-02, 1.5944487401000001E-02, 8.9556174000999994E-02, 2.8660773698399999E-01, 3.1080596191999998E-01, -3.8514514288899998E-01, -4.5879056159300002E-01, -3.1210400792900000E-01, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 6.4795130000999998E-02, 1.6743138021199999E-01, 2.7671888066799999E-01, 4.0807795391800000E-01, 4.8786781696600001E-01, 5.9214403791400005E-01, 6.7033029886899997E-01, 8.6274218284600002E-01, 9.8129062245200005E-01, 9.9872477410700000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9347449523838105E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.0600000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.4300000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.2000000000000000E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.4300000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 24W</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
  Basic characteristics are the material class: paraffin-based, and encapsulation: none<br>  The data is taken from: Croda datasheet - last access 2023-02-28.<br><br>
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
end Croda_Crodatherm_24W;