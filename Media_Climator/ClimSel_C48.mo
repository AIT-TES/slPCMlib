
within slPCMlib.Media_Climator;
package ClimSel_C48 "Climator Sweden AB, ClimSel C48; data taken from: Climator Sweden AB datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C48";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1514999999999998E+02, 3.2814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.1314999999999998E+02, 3.2314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {5.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.5000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 9.4159692493327893E+04
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
    constant Real[10] data_x =   {4.2000000000000000E+01, 4.4125000000000000E+01, 4.7625000000000000E+01, 4.8625000000000000E+01, 4.9625000000000000E+01, 5.0375000000000000E+01, 5.1125000000000000E+01, 5.2125000000000000E+01, 5.4625000000000000E+01, 5.5000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.2701610242999999E-02, 9.6617650222000001E-02, 1.4359991332800001E-01, 1.4244811284900000E-01, 1.5338145773299999E-01, 1.8263410789199999E-01, 1.0463301929000000E-01, 3.9787786579999996E-03, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.6088079629000000E-02, 1.4593899236000000E-02, 6.2595836049000000E-02, -3.2227276318000000E-02, 6.1873867874000003E-02, -3.1984590246999997E-02, -7.6351795158000002E-02, -1.6608999433999998E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.8343534006999999E-02, 2.3190303841599999E-01, 3.4979219404799999E-01, 5.0303259487600005E-01, 6.1119126386000000E-01, 7.4359649454700005E-01, 8.9318665563800004E-01, 9.9944016019000004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0153350823725715E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {4.0000000000000000E+01, 4.1375000000000000E+01, 4.3875000000000000E+01, 4.4875000000000000E+01, 4.5375000000000000E+01, 4.5875000000000000E+01, 4.6375000000000000E+01, 4.6625000000000000E+01, 4.7625000000000000E+01, 4.8375000000000000E+01, 4.8875000000000000E+01, 4.9875000000000000E+01, 5.0000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.6860125653000001E-02, 3.8291603591999998E-02, 1.2091555402000000E-01, 1.3422748076999999E-01, 1.2836897802399999E-01, 2.0874747397900001E-01, 3.0682774621400000E-01, 3.7466726228000002E-01, 9.0693417565000004E-02, 6.6450488630000002E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -1.6103463855999998E-02, 3.7477127558999998E-02, 7.8874273401999997E-02, -1.8451396435999998E-02, 5.6541007672000002E-02, 3.1969178459399999E-01, 3.4950497023299998E-01, -3.0439627208199999E-01, -3.8807998910899999E-01, -1.8459861827000001E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.4179594132000000E-02, 5.5361174768000003E-02, 1.3179055182300001E-01, 1.9784207245999999E-01, 2.6216074155699998E-01, 3.4124267068600000E-01, 4.0576693360499999E-01, 8.0243637003000001E-01, 9.8151495364600005E-01, 9.9820934130100003E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0036184713338085E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.3000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.3000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 7.6000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 5.3000000000000003E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C48</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Climator Sweden AB datasheet - last access 2022-10-14.<br><br>
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
end ClimSel_C48;