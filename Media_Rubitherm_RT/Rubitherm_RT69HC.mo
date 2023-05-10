
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT69HC "Rubitherm GmbH, RT69HC; data taken from: Rubitherm datasheet; last access: 2020-10-09."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT69HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.3714999999999998E+02, 3.4414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.3714999999999998E+02, 3.4314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9060702503406294E+05
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
    constant Real[9] data_x =   {6.4000000000000000E+01, 6.6375000000000000E+01, 6.7625000000000000E+01, 6.8375000000000000E+01, 6.8625000000000000E+01, 6.9125000000000000E+01, 6.9375000000000000E+01, 6.9625000000000000E+01, 7.1000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 2.2913050051999999E-02, 8.2959649258999998E-02, 2.7941040478599999E-01, 4.6006589331600001E-01, 6.1417046638899997E-01, 5.3196133529199996E-01, 3.6998934404900002E-01, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 2.7252168481000000E-02, 7.3173403356000005E-02, 5.5514604554799996E-01, 6.2906840679099996E-01, -1.3414800463500001E-01, -5.6915618659000000E-01, -5.4398647632700003E-01, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 1.4527117448000001E-02, 7.5252486775000002E-02, 1.8955441387100000E-01, 2.8242097565899998E-01, 5.6940528559799997E-01, 7.1622918000600000E-01, 8.2984147979699996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0088760569198769E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {6.4000000000000000E+01, 6.5375000000000000E+01, 6.6375000000000000E+01, 6.6875000000000000E+01, 6.7375000000000000E+01, 6.7625000000000000E+01, 6.8125000000000000E+01, 6.8375000000000000E+01, 6.8625000000000000E+01, 6.9125000000000000E+01, 6.9375000000000000E+01, 6.9625000000000000E+01, 6.9875000000000000E+01, 7.0000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 3.0594448768999999E-02, 6.4849704760000002E-02, 1.4095074792699999E-01, 1.4168318239600000E-01, 1.0902730821800000E-01, 1.5269964227899999E-01, 2.5657974904399999E-01, 4.8925566718500002E-01, 6.3896037115700000E-01, 4.7953585543600002E-01, 2.4599410561999999E-01, 5.7118340891000000E-02, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, 1.5314427290000000E-03, 1.3393721150800000E-01, 1.0915185615400000E-01, -8.1496038074999999E-02, -7.9883197025000005E-02, 2.1972561549500000E-01, 6.4446193123100004E-01, 7.8842520186200005E-01, -4.3863809911399998E-01, -8.6801889011000000E-01, -8.4514841266200003E-01, -6.7158179487400005E-01, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 2.0877079555999999E-02, 5.7714756907000002E-02, 1.0989286752700000E-01, 1.8482711665099999E-01, 2.1628512231700001E-01, 2.7571606324199999E-01, 3.2486316085099998E-01, 4.1771940671399999E-01, 7.2659002294300001E-01, 8.6921690819700004E-01, 9.6015789545100005E-01, 9.9729358184899997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0040725463816822E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.4000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.4000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.4000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT69HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
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
end Rubitherm_RT69HC;