within slPCMlib.Media_Croda_Crodatherm;
package Croda_Crodatherm_47 "Croda International Plc, Crodatherm 47; data taken from: Croda datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 47";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1014999999999998E+02, 3.2214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.1114999999999998E+02, 3.2014999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.9000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8122517309288488E+05
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
    constant Real[11] data_x =   {3.7000000000000000E+01, 3.9125000000000000E+01, 4.1625000000000000E+01, 4.3875000000000000E+01, 4.5375000000000000E+01, 4.5625000000000000E+01, 4.6375000000000000E+01, 4.7125000000000000E+01, 4.7875000000000000E+01, 4.8625000000000000E+01, 4.9000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.3681545060000001E-03, 1.7537251420000000E-02, 8.5679631625000002E-02, 2.8961528135300002E-01, 3.9003756067099998E-01, 3.9131528452100001E-01, 1.2184783657700000E-01, 7.2788888390000004E-03, 7.8371857999999999E-05, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.8035605189999999E-03, 7.0522493900000000E-04, 4.7858082964000001E-02, 3.2115951564500000E-01, 3.4321780451500000E-01, -3.7735987479100003E-01, -2.8680882411899999E-01, -2.2667260780999999E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.8535797280000001E-03, 2.7530401278000001E-02, 1.2459987330300000E-01, 3.5684418223199998E-01, 4.4242924404500000E-01, 7.7207677184599999E-01, 9.6191720850899998E-01, 9.9827386815899999E-01, 9.9998517652800001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0087615326043198E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {3.8000000000000000E+01, 4.1125000000000000E+01, 4.3875000000000000E+01, 4.4625000000000000E+01, 4.5125000000000000E+01, 4.5375000000000000E+01, 4.5625000000000000E+01, 4.6125000000000000E+01, 4.6375000000000000E+01, 4.6625000000000000E+01, 4.6875000000000000E+01, 4.7000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.3360036383000000E-02, 7.5688499849000002E-02, 9.1116000783000001E-02, 1.7387411519700000E-01, 2.9218763478400001E-01, 5.2522814139999996E-01, 6.5248649796000002E-01, 4.8379834185699999E-01, 2.4632947452199999E-01, 5.6972326764999998E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 9.4533365079999994E-03, 3.5282351199999999E-02, 1.6838646249000000E-02, 2.8676793241300003E-01, 6.5813335750799995E-01, 7.8326106929499995E-01, -4.9204217487899998E-01, -8.8377974206599996E-01, -8.5807757135800000E-01, -6.6787395753099998E-01, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.3220864077999999E-02, 1.1969865292500000E-01, 1.8330228250700001E-01, 2.4410542691000001E-01, 3.0059538514599998E-01, 4.0242065460900001E-01, 7.2436667045400005E-01, 8.6886831546300003E-01, 9.6026971983599996E-01, 9.9730090487500000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0029549715656356E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 9.4000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 8.2900000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 2.5000000000000000E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 1.6000000000000000E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 47</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
  Basic characteristics are the material class: paraffin-based, and encapsulation: none<br>  The data is taken from: Croda datasheet - last access 2023-02-28.<br><br>
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
end Croda_Crodatherm_47;
