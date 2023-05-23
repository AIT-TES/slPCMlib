
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS3N "Pluss Advanced Technologies Pvt Ltd, HS3N; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS3N";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6714999999999998E+02, 2.7514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6614999999999998E+02, 2.7114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.8000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.9800000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.3027467611267889E+05
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
    constant Real[11] data_x =   {-6.0000000000000000E+00, -4.3750000000000000E+00, -2.8750000000000000E+00, -2.6250000000000000E+00, -2.3750000000000000E+00, -1.8750000000000000E+00, -1.6250000000000000E+00, -1.3750000000000000E+00, -3.7500000000000000E-01, 1.3750000000000000E+00, 2.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 6.5304539480000000E-03, 1.9734495857900000E-01, 3.2752484285900002E-01, 5.3278806621700003E-01, 6.5763551116700003E-01, 5.2650479936200001E-01, 3.1993272837999998E-01, 5.1957045271000001E-02, 2.9615581199999999E-04, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 3.8087991610000001E-03, 3.8477643027699998E-01, 6.3490461316199998E-01, 7.0666461411299997E-01, -3.3518473770399998E-01, -7.0804398050100004E-01, -6.2629518042599996E-01, -8.2370572401000000E-02, -9.8244757099999989E-04, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 4.4781914010000000E-03, 8.6141728143000001E-02, 1.5059640752500000E-01, 2.5800958704100002E-01, 5.7805907114099997E-01, 7.2836536093399995E-01, 8.3398796169500000E-01, 9.7493097132600004E-01, 9.9993929196300002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0023124625140387E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-7.0000000000000000E+00, -5.8750000000000000E+00, -4.6250000000000000E+00, -4.1250000000000000E+00, -3.8750000000000000E+00, -3.6250000000000000E+00, -3.3750000000000000E+00, -3.1250000000000000E+00, -2.8750000000000000E+00, -2.6250000000000000E+00, -2.3750000000000000E+00, -2.1250000000000000E+00, -2.0000000000000000E+00};
    constant Real[13] data_y =   {0.0000000000000000E+00, 9.9496165740000000E-03, 4.8647967032000000E-02, 5.5583371049999997E-02, 1.0288037624400000E-01, 2.2891336409499999E-01, 5.5406471037500005E-01, 8.1713746581299995E-01, 8.6914368288699995E-01, 6.8059368524999997E-01, 3.5834397826100001E-01, 8.4344431447000004E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 2.2168352940000000E-02, -7.0749012219999997E-03, 8.2618169431999997E-02, 2.3463575068900000E-01, 9.4947732231199999E-01, 1.1800647489460001E+00, 8.7910917967199997E-01, -3.8105098071400001E-01, -1.1996513210780000E+00, -1.1810364535450000E+00, -1.0023081301060000E+00, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 3.2659182240000000E-03, 4.3788038305999998E-02, 6.8031659886999996E-02, 8.7090630540999997E-02, 1.2492660323300001E-01, 2.2181526460199999E-01, 3.9517195454199999E-01, 6.1300916926399995E-01, 8.1143506969400003E-01, 9.4149713087999998E-01, 9.9602464268199997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0022486330282869E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.8500000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.0600000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.2000000000000002E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 3.4999999999999998E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS3N</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
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
end PLUSS_savE_HS3N;