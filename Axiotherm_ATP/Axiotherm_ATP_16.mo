
// within slPCMlib.Axiotherm_ATP;
package Axiotherm_ATP_16 "Axiotherm GmbH, ATP 16; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 16";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8614999999999998E+02, 2.9214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8514999999999998E+02, 2.9014999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4345475485984306E+05
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
    constant Real[10] data_x =   {1.3000000000000000E+01, 1.4625000000000000E+01, 1.5375000000000000E+01, 1.5625000000000000E+01, 1.5875000000000000E+01, 1.6125000000000000E+01, 1.6625000000000000E+01, 1.6875000000000000E+01, 1.7875000000000000E+01, 1.9000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.1045447672000001E-01, 4.6387471335699998E-01, 5.9827725821199995E-01, 6.1255097453100005E-01, 5.0065527784099995E-01, 1.2070274984899999E-01, 5.0109227855999999E-02, 1.4718291517000000E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 2.1528189917900001E-01, 4.4048691635800002E-01, 4.4814510989000000E-01, -2.4752670797000001E-01, -7.5914240723399995E-01, -7.4001396904200001E-01, -1.2670724481400000E-01, -1.1474298766000000E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.2381450050999999E-01, 3.6651069094400002E-01, 4.9944752382599999E-01, 6.5466688446400001E-01, 7.9670427421099999E-01, 9.5188775677100002E-01, 9.7007336438000002E-01, 9.9292007810899996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0015650135825194E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {1.2000000000000000E+01, 1.3625000000000000E+01, 1.4625000000000000E+01, 1.5625000000000000E+01, 1.5875000000000000E+01, 1.6625000000000000E+01, 1.6875000000000000E+01, 1.7000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 1.3673213147200000E-01, 3.0543972237700001E-01, 4.4217589259399998E-01, 4.0049454450700001E-01, 9.1276936816000007E-02, 2.0033550414000002E-02, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 1.1047702849500000E-01, 2.0497149795399999E-01, 2.0401291514000001E-02, -3.0323347088000002E-01, -4.0081904111500000E-01, -2.2662331762999999E-01, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 8.6857482513999995E-02, 3.0024901762900003E-01, 6.8976651029500002E-01, 7.9687633586100004E-01, 9.8602464589700001E-01, 9.9904217693999997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0008449235024748E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.5600000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.7000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.7000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 16</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    macroencapsulation<br>  Data taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATP_16;