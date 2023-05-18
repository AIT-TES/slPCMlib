
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT70HC "Rubitherm GmbH, RT70HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT70HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.3814999999999998E+02, 3.4514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.3614999999999998E+02, 3.4414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.1600000000000000E+05
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
    constant Real[11] data_x =   {6.5000000000000000E+01, 6.7125000000000000E+01, 6.8625000000000000E+01, 6.9125000000000000E+01, 6.9375000000000000E+01, 6.9625000000000000E+01, 6.9875000000000000E+01, 7.0375000000000000E+01, 7.0625000000000000E+01, 7.1375000000000000E+01, 7.2000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 9.7464334539999996E-03, 3.0668601699999999E-02, 9.5287718785999995E-02, 2.1509034273800001E-01, 5.1353468998700003E-01, 7.5819692261799998E-01, 6.5367614902899995E-01, 3.6166252252800002E-01, 4.6016593015000003E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 8.9344396970000000E-03, 1.0450864937000001E-02, 2.3172056674200001E-01, 8.8910679726399999E-01, 1.0888230262270000E+00, 8.3064483249099996E-01, -9.7844336466699999E-01, -8.4567907903100004E-01, -1.4237129787600000E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 7.0204713299999999E-03, 3.7163060083000000E-02, 6.4145874386999993E-02, 9.9655478147000001E-02, 1.9004017968699999E-01, 3.5096871208699998E-01, 7.4313084849300004E-01, 8.6984283263899997E-01, 9.9021676362400002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0038512780305910E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    17;
    constant Real[17] data_x =   {6.3000000000000000E+01, 6.5625000000000000E+01, 6.6375000000000000E+01, 6.6625000000000000E+01, 6.7125000000000000E+01, 6.7625000000000000E+01, 6.7875000000000000E+01, 6.8875000000000000E+01, 6.9125000000000000E+01, 6.9375000000000000E+01, 6.9625000000000000E+01, 6.9875000000000000E+01, 7.0125000000000000E+01, 7.0375000000000000E+01, 7.0625000000000000E+01, 7.0875000000000000E+01, 7.1000000000000000E+01};
    constant Real[17] data_y =   {0.0000000000000000E+00, 2.8364162974000000E-02, 1.3234545014800000E-01, 2.5934337297700000E-01, 2.7736841426699999E-01, 3.9476306418999997E-02, 1.1002603111999999E-02, 1.9783361879999999E-02, 5.1238870479999997E-02, 1.3554581883799999E-01, 3.6145413653199998E-01, 5.5631144622499995E-01, 6.0611073765900003E-01, 4.8167564838600002E-01, 2.5588917468500000E-01, 6.0510076025999997E-02, 0.0000000000000000E+00};
    constant Real[17] m_k =      {0.0000000000000000E+00, 1.8561742270999999E-02, 3.3273622793800001E-01, 4.0939132458700001E-01, -4.7249495450500001E-01, -4.0952140485799998E-01, -3.2453918734999997E-02, 5.6570128940999997E-02, 1.5168535446100001E-01, 6.9401529801899997E-01, 8.4228010472000003E-01, 6.9166463528300004E-01, -2.1108091135900001E-01, -8.4290736633600005E-01, -8.3298884623199998E-01, -7.2178985956200004E-01, 0.0000000000000000E+00};
    constant Real[17] iy_start = {0.0000000000000000E+00, 2.6772135656999999E-02, 7.2658684963000000E-02, 1.2159097301400000E-01, 2.7530520834900002E-01, 3.5379865587699999E-01, 3.5817777752899999E-01, 3.6621291753500002E-01, 3.7465924573499998E-01, 3.9533924991800001E-01, 4.5716003158899998E-01, 5.7354625975499995E-01, 7.2469506921399995E-01, 8.6502143825099997E-01, 9.5786825732900005E-01, 9.9713627160800000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0076280124731147E+00;
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
    lambda := 2.0000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT70HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
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
end Rubitherm_RT70HC;