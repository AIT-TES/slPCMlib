within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_50 "Axiotherm GmbH, ATS 50; data taken from: Axiotherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS 50";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1814999999999998E+02, 3.2614999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.1514999999999998E+02, 3.2314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9101017329663254E+05
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
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {4.5000000000000000E+01, 4.7375000000000000E+01, 4.8125000000000000E+01, 4.8625000000000000E+01, 4.9375000000000000E+01, 4.9875000000000000E+01, 5.0125000000000000E+01, 5.0625000000000000E+01, 5.0875000000000000E+01, 5.1875000000000000E+01, 5.3000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.4068340883000000E-02, 2.9226752509999999E-02, 1.7758830613200000E-01, 5.8827596884199995E-01, 6.3204070685699998E-01, 4.7926570612100000E-01, 9.3098495205999998E-02, 3.3177641922999997E-02, 5.8013765209999997E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -7.6677833940000002E-03, 6.9804159015999998E-02, 5.4591952748299999E-01, 4.8407395807300002E-01, -4.8029145584900002E-01, -7.8364872685300002E-01, -5.3559336576500005E-01, -9.5387884583000002E-02, -5.9279831090000004E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 2.0380759166999999E-02, 3.3028577458000002E-02, 7.4957996720000003E-02, 3.6606088993600000E-01, 6.9235725882099997E-01, 8.3333715645799999E-01, 9.7173809538900002E-01, 9.8527660555600005E-01, 9.9735280554399997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0034635834748877E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {4.2000000000000000E+01, 4.3625000000000000E+01, 4.4625000000000000E+01, 4.5625000000000000E+01, 4.6125000000000000E+01, 4.6375000000000000E+01, 4.6625000000000000E+01, 4.6875000000000000E+01, 4.7125000000000000E+01, 4.7375000000000000E+01, 4.7625000000000000E+01, 4.7875000000000000E+01, 4.8375000000000000E+01, 4.9375000000000000E+01, 5.0000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 4.0156126817000003E-02, 8.2185728695000002E-02, 4.1300136299999998E-02, 9.7072938622999996E-02, 2.5228756033299998E-01, 7.2836999352900000E-01, 9.9664078930400002E-01, 7.9638458895700004E-01, 3.3103021494200002E-01, 1.9071866138000000E-02, 4.2681491600000000E-04, 4.2681491600000000E-04, 1.6207910297999999E-02, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, 2.1529199504000000E-02, 5.1668701764999998E-02, -7.4819639849999994E-02, 2.4033587464900000E-01, 1.0574036157930000E+00, 1.5664503109390000E+00, 2.4954042752700001E-01, -1.5998738339900000E+00, -1.5701268745390000E+00, -2.0173219837199999E-01, 0.0000000000000000E+00, 0.0000000000000000E+00, -2.6248414237999999E-02, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 2.7996424730999999E-02, 8.6881025195000003E-02, 1.5944227885600001E-01, 1.8757745547400001E-01, 2.2714333940299999E-01, 3.4753487786199999E-01, 5.7087464986600001E-01, 8.0553301256900001E-01, 9.4684560595400002E-01, 9.8362202042699998E-01, 9.8501399284000002E-01, 9.8522821995099996E-01, 9.9577329684500004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0038407871667536E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.3830000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.3000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATS 50</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATS_50;
