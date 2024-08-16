within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS36 "Pluss Advanced Technologies Pvt Ltd, HS36; data taken from: PLUSS datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS36";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.0314999999999998E+02, 3.1414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.0214999999999998E+02, 3.0914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.9800000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.3200000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4204108465131599E+05
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
    constant Real[10] data_x =   {3.0000000000000000E+01, 3.2125000000000000E+01, 3.4625000000000000E+01, 3.5375000000000000E+01, 3.5625000000000000E+01, 3.6375000000000000E+01, 3.6625000000000000E+01, 3.7875000000000000E+01, 3.9625000000000000E+01, 4.1000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 1.9252089994000000E-02, 1.2875526816999999E-01, 2.6563591786800000E-01, 3.6123525327400002E-01, 3.1558068070900003E-01, 2.0573904366699999E-01, 7.8621923243000005E-02, 2.9021356804999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.2960902516000000E-02, 7.2031825430000002E-02, 2.9708956501200001E-01, 3.1878006974299999E-01, -3.7852954168399999E-01, -3.3807982033000000E-01, -1.6688540716000001E-02, -4.4913586193000003E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.5783221550999999E-02, 1.7205696799900000E-01, 3.1121229525400002E-01, 3.9048835025299999E-01, 6.8075246744700002E-01, 7.4656189916799998E-01, 8.8422855048100002E-01, 9.8695453065700001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0131652860510860E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {2.9000000000000000E+01, 3.0375000000000000E+01, 3.3375000000000000E+01, 3.4375000000000000E+01, 3.4625000000000000E+01, 3.5125000000000000E+01, 3.5375000000000000E+01, 3.5625000000000000E+01, 3.5875000000000000E+01, 3.6000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 6.0533754829000003E-02, 1.2204402518700000E-01, 3.1675985674599999E-01, 4.6079859710299997E-01, 4.7042499221000000E-01, 3.3191575287699998E-01, 1.6382157646100001E-01, 3.7271881456000001E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, -1.1282429520000001E-03, 5.3551075283000002E-02, 4.2438602396399999E-01, 4.6545406888500002E-01, -4.7890551621299998E-01, -6.3455789304700005E-01, -6.1212041075800006E-01, -4.3175336165299999E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 4.1730862725999997E-02, 2.7423230279299998E-01, 4.6244335706000000E-01, 5.5927610775199998E-01, 8.1137044350400001E-01, 9.1231926828800003E-01, 9.7407508243300001E-01, 9.9823538623100005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9847227695754459E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.9670000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.8500000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 5.0000000000000000E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 4.6999999999999997E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS36</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
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
end PLUSS_savE_HS36;
