within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP25E2 "Rubitherm GmbH, SP25E2; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP25E2";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9114999999999998E+02, 3.0514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8914999999999998E+02, 2.9914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4280797633639537E+05
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
    constant Real[13] data_x =   {1.8000000000000000E+01, 2.0625000000000000E+01, 2.2625000000000000E+01, 2.3375000000000000E+01, 2.3875000000000000E+01, 2.4375000000000000E+01, 2.4625000000000000E+01, 2.5375000000000000E+01, 2.5875000000000000E+01, 2.7375000000000000E+01, 2.9875000000000000E+01, 3.1375000000000000E+01, 3.2000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.6034705672999999E-02, 1.1404441231300000E-01, 1.4623420729100001E-01, 1.4760789803600000E-01, 2.2136396592999999E-01, 3.0292842476600002E-01, 3.1813268712600001E-01, 1.7984512714699999E-01, 1.9063828357999999E-02, 1.2480096142000000E-02, 1.3090199473000001E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 2.4809968289999998E-03, 8.8414391007000007E-02, -4.5859119629999999E-03, 6.1027965566999999E-02, 2.4888349706999999E-01, 2.7234817432899999E-01, -2.7137956235999999E-01, -2.2318425089700000E-01, -1.6505492954999999E-02, 2.6996526070000000E-03, -2.6046602822000001E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.9652373207000001E-02, 1.2124964989299999E-01, 2.2337699333200001E-01, 2.9558614670299999E-01, 3.8405706828699998E-01, 4.4957627971399999E-01, 7.0837568604400003E-01, 8.3206405453300003E-01, 9.4267054823100005E-01, 9.7214495258699996E-01, 9.9675198485899996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0016032334105387E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {1.6000000000000000E+01, 2.0625000000000000E+01, 2.1875000000000000E+01, 2.2375000000000000E+01, 2.2625000000000000E+01, 2.3375000000000000E+01, 2.4125000000000000E+01, 2.4625000000000000E+01, 2.5125000000000000E+01, 2.6000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 3.4972453932999997E-02, 9.4680827019000005E-02, 2.2259209981500000E-01, 3.5868211768800001E-01, 5.1976358661500000E-01, 2.0628684194600000E-01, 1.6643258804999998E-02, 2.4167383520000002E-03, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, -7.5565729000000000E-05, 1.3142302338000000E-01, 4.5145104740899999E-01, 4.9575994114499999E-01, -2.0792643255300000E-01, -4.9271776600599998E-01, -8.4793616387999995E-02, -1.8995445450000000E-03, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 8.1890601684999995E-02, 1.4649761844500001E-01, 2.1993969669499999E-01, 2.9315687285499997E-01, 6.5950551341100005E-01, 9.4823411976500005E-01, 9.9598255767300004E-01, 9.9905367836500003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0108890069788461E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.6000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.5000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 5.0000000000000000E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 5.0000000000000000E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP25E2</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP25E2;
