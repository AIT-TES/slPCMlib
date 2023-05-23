
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS15N "Pluss Advanced Technologies Pvt Ltd, HS15N; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS15N";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.5314999999999998E+02, 2.6314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.5214999999999998E+02, 2.5914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.8700000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.0212812212875817E+05
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
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-2.0000000000000000E+01, -1.8125000000000000E+01, -1.6375000000000000E+01, -1.5875000000000000E+01, -1.5625000000000000E+01, -1.5375000000000000E+01, -1.4875000000000000E+01, -1.4625000000000000E+01, -1.4375000000000000E+01, -1.3375000000000000E+01, -1.1625000000000000E+01, -1.1125000000000000E+01, -1.0000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.2909177302000001E-02, 5.2802708171000003E-02, 1.1682585252100000E-01, 2.0846105238099999E-01, 3.9384794902300002E-01, 5.4357331273599996E-01, 4.3870484894299999E-01, 2.5977188294900000E-01, 7.2955288406000005E-02, 5.9114613046999998E-02, 9.3559219238999999E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 6.2137726749999997E-03, 2.2460964985000002E-02, 2.1698392745299999E-01, 5.4270448664500004E-01, 6.4543970338099999E-01, -2.2680872346200001E-01, -5.9798851405900000E-01, -5.1525715130700001E-01, -3.5771863852999998E-02, 6.8521940811999998E-02, -5.8074279199999996E-04, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.0413595794000001E-02, 6.4448363827999997E-02, 1.0329415993700000E-01, 1.4275758963499999E-01, 2.1846852481699999E-01, 4.7422986737000000E-01, 6.0054516342399999E-01, 6.8853654183199997E-01, 8.1656194282000005E-01, 9.0664559544500001E-01, 9.4676096457400005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0128072755349669E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {-2.1000000000000000E+01, -1.8375000000000000E+01, -1.7875000000000000E+01, -1.7625000000000000E+01, -1.7375000000000000E+01, -1.7125000000000000E+01, -1.6875000000000000E+01, -1.6625000000000000E+01, -1.6375000000000000E+01, -1.6125000000000000E+01, -1.5625000000000000E+01, -1.4000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 3.9862365688000001E-02, 1.4354827991800001E-01, 3.4497516437400000E-01, 8.6788393732400004E-01, 1.1151848698210001E+00, 8.7089295058500005E-01, 3.6245737990700000E-01, 2.7621805172000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.8666871838999999E-02, 3.7314731301600002E-01, 1.2065263659979999E+00, 1.6557443936930001E+00, -1.4864596137000001E-02, -1.7580006719660000E+00, -1.7127149528290000E+00, -2.7917269794600003E-01, 0.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 4.1716911680000003E-02, 8.0292220760000005E-02, 1.3717588765600000E-01, 2.8686134379400002E-01, 5.4416411974900003E-01, 8.0222287741200005E-01, 9.5658660749199997E-01, 9.9799570518900005E-01, 1.0000000000000000E+00, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0027986337557775E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.0160000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.0700000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.2599999999999998E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 5.3000000000000003E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS15N</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
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
end PLUSS_savE_HS15N;