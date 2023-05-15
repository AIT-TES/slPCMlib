
within slPCMlib.Media_Axiotherm_ATP;
package Axiotherm_ATP_70 "Axiotherm GmbH, ATP 70; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 70";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.3514999999999998E+02, 3.4614999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.3614999999999998E+02, 3.4614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.0900000000000000E+05
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
    constant Real[12] data_x =   {6.2000000000000000E+01, 6.5125000000000000E+01, 6.7125000000000000E+01, 6.8625000000000000E+01, 6.9625000000000000E+01, 7.0375000000000000E+01, 7.0625000000000000E+01, 7.0875000000000000E+01, 7.1125000000000000E+01, 7.1625000000000000E+01, 7.1875000000000000E+01, 7.3000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 4.8314490920000003E-03, 9.4261789530000004E-03, 2.5302391877000001E-02, 1.8146335294400001E-01, 5.7718916913200002E-01, 7.4371451936199995E-01, 6.9420511172300003E-01, 4.6225626296700001E-01, 1.9138755981000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 3.3461587380000002E-03, 4.6266124450000000E-03, 9.1683357699999992E-03, 4.0363579962000001E-01, 6.1736118557899999E-01, 6.0727380303300005E-01, -8.0347878989699995E-01, -1.0386796168300001E+00, -1.6443780870200000E-01, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 4.8720905760000002E-03, 1.8834909734000000E-02, 4.4270238571000001E-02, 1.1545381491700001E-01, 3.9254988927500001E-01, 5.5929182102200004E-01, 7.4816511395600005E-01, 8.9533922178799996E-01, 9.9844944295399995E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0095445613032314E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {6.3000000000000000E+01, 6.5875000000000000E+01, 6.7875000000000000E+01, 6.8875000000000000E+01, 6.9125000000000000E+01, 6.9375000000000000E+01, 6.9625000000000000E+01, 6.9875000000000000E+01, 7.0125000000000000E+01, 7.0375000000000000E+01, 7.0625000000000000E+01, 7.0875000000000000E+01, 7.2125000000000000E+01, 7.3000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 1.5062398374000001E-02, 8.1501213693999999E-02, 4.8291995275999998E-02, 9.6780008304999998E-02, 2.4751733731700001E-01, 7.1301610616500000E-01, 9.7034233616700005E-01, 8.1472723636800004E-01, 3.0583588992700000E-01, 2.8365695619000000E-02, 5.4249953279999997E-03, 5.1065008060000002E-03, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, 1.4568406060000001E-02, -3.5082362578000001E-02, 6.9770727901999993E-02, 2.3902008719900000E-01, 1.0176642134730001E+00, 1.5018222642100001E+00, 4.1686761679599998E-01, -1.5280523116110001E+00, -1.4975277164900000E+00, -1.8358608019500000E-01, -1.9149996940000000E-03, 3.5016533910000001E-03, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 1.1729362554999999E-02, 1.2593288057700000E-01, 1.8263271753499999E-01, 2.0005140780899999E-01, 2.3940865707799999E-01, 3.5808601699699999E-01, 5.7623808149199995E-01, 8.1174864571100003E-01, 9.5300785032799995E-01, 9.8827610743200001E-01, 9.9158531438399999E-01, 9.9751881943400000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0096332108330397E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8900000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.0000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.0000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 70</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATP_70;