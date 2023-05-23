
// within slPCMlib.PLUSS_HS;
package PLUSS_savE_HS78 "Pluss Advanced Technolgies Pvt Ltd, HS78; data taken from: PLUSS datasheet; last access: 2022-02-13."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS78";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.4514999999999998E+02, 3.5614999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.4814999999999998E+02, 3.5114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.4600000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.5900000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.4821681589502440E+05
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
    constant Real[10] data_x =   {7.2000000000000000E+01, 7.3625000000000000E+01, 7.6375000000000000E+01, 7.7125000000000000E+01, 7.7375000000000000E+01, 7.8125000000000000E+01, 7.9625000000000000E+01, 8.0625000000000000E+01, 8.2375000000000000E+01, 8.3000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 1.8105165573000001E-02, 5.0475105605999999E-02, 1.0988793466800000E-01, 1.6061192725699999E-01, 3.6867601190299998E-01, 1.6242850687900001E-01, 4.3999411904000002E-02, 8.0919384669999993E-03, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.1308948010000000E-03, 3.1542714585999998E-02, 1.4686642346400000E-01, 2.9092888854799998E-01, 1.6311028944500000E-01, -2.5016260710999999E-01, -4.9462316621999999E-02, -2.0321025887999999E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.4694100171999999E-02, 9.1034140828000001E-02, 1.4664441923300001E-01, 1.8023813956800000E-01, 3.8800008965600002E-01, 8.7146718334399997E-01, 9.5934666033600002E-01, 9.9810274007699995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0160776844417558E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {7.5000000000000000E+01, 7.6125000000000000E+01, 7.6375000000000000E+01, 7.6875000000000000E+01, 7.7125000000000000E+01, 7.7375000000000000E+01, 7.7875000000000000E+01, 7.8000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 0.0000000000000000E+00, 1.5448486467000000E-02, 7.7634757969900003E-01, 1.1143881623400000E+00, 1.0429937334099999E+00, 1.4843796105900001E-01, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 0.0000000000000000E+00, 1.3772615095499999E-01, 1.7349662355310000E+00, 8.9237680487100002E-01, -1.6722708571130001E+00, -1.7885766572610000E+00, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 0.0000000000000000E+00, 1.2175393500000000E-03, 1.6640658841700001E-01, 4.0789117315700002E-01, 6.9180809201100002E-01, 9.9302973590499999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0031326753780236E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.8900000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.9000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 6.9999999999999996E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.9000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS78</strong>  from manufacturer: <strong>Pluss Advanced Technolgies Pvt Ltd</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    multiple options available<br>  Data taken from: PLUSS datasheet - last access 2022-02-13.<br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-04-20.<br><br>
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
    <li>file creation date: 2023-04-20 </ul>
    </p></html>"));
end PLUSS_savE_HS78;