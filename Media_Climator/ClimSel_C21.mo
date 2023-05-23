
// within slPCMlib.Media_Climator;
package ClimSel_C21 "Climator Sweden AB, ClimSel C21; data taken from: Climator Sweden AB datasheet; last access: 2022-10-14."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C21";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8914999999999998E+02, 3.0114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8914999999999998E+02, 2.9814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {4.4000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 8.3799077618098207E+04
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
    constant Real[12] data_x =   {1.6000000000000000E+01, 1.7625000000000000E+01, 1.9375000000000000E+01, 2.0125000000000000E+01, 2.0375000000000000E+01, 2.1375000000000000E+01, 2.1625000000000000E+01, 2.3375000000000000E+01, 2.4625000000000000E+01, 2.5625000000000000E+01, 2.7125000000000000E+01, 2.8000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.0504616333999999E-02, 1.8717216145999999E-02, 1.9615828581999999E-02, 4.5472308784000000E-02, 2.4225757614400001E-01, 1.8931804486199999E-01, 1.3002313001999999E-01, 1.8144159288600001E-01, 1.1231067996700000E-01, 3.9107402775000001E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -4.7456005059999997E-03, -2.3750566870000000E-02, 4.1961992261000003E-02, 2.7962267777200001E-01, -1.5868616218000001E-01, -1.5311977864800000E-01, 1.3506842587999999E-02, 5.6503798112000000E-02, -1.1277896755000000E-01, -5.9002305434999999E-02, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 9.6151013790000007E-03, 4.0148182440999998E-02, 5.1485033763000003E-02, 5.8409030408000003E-02, 2.3947427318400000E-01, 2.9359385826399997E-01, 5.3137875831699999E-01, 7.2115265620500002E-01, 8.8273767931400005E-01, 9.8660507813599996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0037394595592224E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {1.6000000000000000E+01, 1.7375000000000000E+01, 1.8875000000000000E+01, 1.9875000000000000E+01, 2.0125000000000000E+01, 2.0375000000000000E+01, 2.0625000000000000E+01, 2.0875000000000000E+01, 2.1125000000000000E+01, 2.1375000000000000E+01, 2.1625000000000000E+01, 2.1875000000000000E+01, 2.2625000000000000E+01, 2.4625000000000000E+01, 2.5000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 3.3375722100999997E-02, 4.7023394747000000E-02, 7.4725300059999999E-02, 1.3675477729200000E-01, 2.8853531571500002E-01, 6.5178416134399997E-01, 8.6506392087600004E-01, 7.7605044818000002E-01, 4.4835131905499997E-01, 1.1544255344300000E-01, 3.1507629183999998E-02, 9.5307139020000004E-03, 3.4648371649999999E-03, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, -7.1767852130000002E-03, 2.0931534203000000E-02, 1.1999858767700000E-01, 3.0524441530000002E-01, 9.2381651909399998E-01, 1.2060815414760000E+00, 3.8653896270200000E-01, -1.2946215063669999E+00, -1.3209959733770000E+00, -1.2555902483840000E+00, -9.9511841840000004E-02, 1.2185796052000000E-02, -1.5163908694000000E-02, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 2.4091463994000001E-02, 7.9154632914000003E-02, 1.3180603832099999E-01, 1.5729202821499999E-01, 2.0726254367800001E-01, 3.2340436015099999E-01, 5.1739910518400001E-01, 7.3142715062700003E-01, 8.8470977937200002E-01, 9.5488687097500002E-01, 9.6724206277299996E-01, 9.7740191413199995E-01, 9.9952775227599999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0006204211889229E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.5000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.5000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 9.3000000000000005E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.5000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C21</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    macroencapsulation<br>  Data taken from: Climator Sweden AB datasheet - last access 2022-10-14.<br><br>
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
end ClimSel_C21;