
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT8HC "Rubitherm GmbH, RT8HC; data taken from: Rubitherm datasheet; last access: 2023-01-02."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT8HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7214999999999998E+02, 2.8514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7214999999999998E+02, 2.8214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.5995607474976432E+05
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
    constant Real[11] data_x =   {-1.0000000000000000E+00, 1.2500000000000000E-01, 1.8750000000000000E+00, 4.1250000000000000E+00, 5.6250000000000000E+00, 7.3750000000000000E+00, 7.8750000000000000E+00, 9.1250000000000000E+00, 9.8750000000000000E+00, 1.1375000000000000E+01, 1.2000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.3817836660000001E-03, 1.5651253067999998E-02, 2.2550582884000001E-02, 7.5668666034000007E-02, 2.5228612659100003E-01, 3.3120070413600000E-01, 1.6921191853199999E-01, 5.3329072959000003E-02, 8.7005907770000002E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 7.6576082980000003E-03, -3.5947048080000000E-03, 1.8314120622000000E-02, 3.6096046590999999E-02, 1.7531981709700001E-01, 7.0812388074000004E-02, -1.9737885621300000E-01, -7.1652890277999995E-02, -2.3511856436000001E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.1069568370000000E-03, 2.0852701171000000E-02, 5.4967338481000001E-02, 1.2609064595100000E-01, 3.8035491603499999E-01, 5.3007312943100005E-01, 8.8167186699800004E-01, 9.6010582539500000E-01, 9.9802439853000002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0112751063594572E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {-1.0000000000000000E+00, 1.8750000000000000E+00, 4.8750000000000000E+00, 6.1250000000000000E+00, 7.1250000000000000E+00, 7.8750000000000000E+00, 8.6250000000000000E+00, 8.8750000000000000E+00, 9.0000000000000000E+00};
    constant Real[9] data_y =   {0.0000000000000000E+00, 2.4547781099999999E-02, 7.3075346144000006E-02, 1.4362701505600001E-01, 3.2924346993300002E-01, 3.6616652797799998E-01, 8.6546115058000000E-02, 1.9074846174000001E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 5.4808214329999999E-03, 2.2985563575000000E-02, 1.1976427734599999E-01, 1.7283909954800000E-01, -2.4715861092200000E-01, -3.8410470658099999E-01, -2.1634610020799999E-01, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 3.1844680611999999E-02, 1.6655716257799999E-01, 2.9069065028600000E-01, 5.2515066541199995E-01, 8.0857565775700002E-01, 9.8662097334099996E-01, 9.9907991744900004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0105497598757063E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8000000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT8HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2023-01-02.<br><br>
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
end Rubitherm_RT8HC;