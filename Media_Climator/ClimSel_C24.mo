within slPCMlib.Media_Climator;
package ClimSel_C24 "Climator Sweden AB, ClimSel C24; data taken from: Climator Sweden AB datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C24";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9214999999999998E+02, 3.0314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8914999999999998E+02, 2.9914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {4.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 8.0699999999999985E+04
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
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {1.9000000000000000E+01, 2.0875000000000000E+01, 2.2875000000000000E+01, 2.4625000000000000E+01, 2.6625000000000000E+01, 2.7375000000000000E+01, 2.7625000000000000E+01, 2.8875000000000000E+01, 3.0000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 3.6669375352000003E-02, 4.4299136189999998E-02, 1.2679552053300000E-01, 3.0507349904100001E-01, 2.0943592043400000E-01, 1.3873707027599999E-01, 2.7587232363000001E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 2.5357179050000001E-02, 2.2827116506000000E-02, 3.9597201520999999E-02, 1.1357136744400000E-01, -2.5849438929700003E-01, -2.3957692734600000E-01, -4.0861897211000003E-02, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 2.6828201507000000E-02, 1.0827431887699999E-01, 2.5305212742599997E-01, 6.5844261322800002E-01, 8.6788369524999998E-01, 9.1111266322999995E-01, 9.8884194227099997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9552939723082057E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {1.6000000000000000E+01, 1.6875000000000000E+01, 1.7875000000000000E+01, 1.8125000000000000E+01, 1.9375000000000000E+01, 2.0625000000000000E+01, 2.2125000000000000E+01, 2.3375000000000000E+01, 2.3625000000000000E+01, 2.3875000000000000E+01, 2.4125000000000000E+01, 2.4625000000000000E+01, 2.4875000000000000E+01, 2.6000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 2.5716522865999999E-02, 1.5318225500000001E-04, 1.5939091700000000E-04, 3.2665204290999998E-02, 4.0771011550000000E-02, 1.2507792018899999E-01, 4.3482336025200002E-01, 5.7132940763600004E-01, 5.5828786178000001E-01, 4.1214306851999999E-01, 3.1385966602999998E-02, 8.0759422639999998E-03, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, -1.0858010031000001E-02, -3.2310894000000000E-04, 3.5444582200000002E-04, -1.8484656366000000E-02, 5.1711006655000002E-02, 1.0092643425099999E-01, 4.2269944478900001E-01, 4.1922514393100002E-01, -3.9475106070299998E-01, -7.9430161994299997E-01, -1.9605117461300001E-01, -1.7814342056999999E-02, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 1.1864694812000001E-02, 2.3841842140000000E-02, 2.3877149622000000E-02, 4.6693517455999999E-02, 8.3207818390999999E-02, 1.9760446579499999E-01, 5.0360652385100002E-01, 6.2856121429000000E-01, 7.7304025206399996E-01, 8.9560850940799996E-01, 9.9337584997600004E-01, 9.9735376917200003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9338167221282658E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.4000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.4000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 7.3999999999999999E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 9.3000000000000005E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C24</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
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
end ClimSel_C24;
