
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT64HC "Rubitherm GmbH, RT64HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT64HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.2814999999999998E+02, 3.4114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.2714999999999998E+02, 3.3814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.3666561674476723E+05
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
    constant Real[12] data_x =   {5.5000000000000000E+01, 5.7625000000000000E+01, 5.9375000000000000E+01, 6.0625000000000000E+01, 6.2875000000000000E+01, 6.3375000000000000E+01, 6.3625000000000000E+01, 6.3875000000000000E+01, 6.4625000000000000E+01, 6.4875000000000000E+01, 6.5875000000000000E+01, 6.8000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 4.8581171739999996E-03, 3.2379590287000000E-02, 5.4236998311999997E-02, 2.5955351056100001E-01, 3.6889233846300001E-01, 4.4362026889900003E-01, 4.3447069433000002E-01, 8.1200603102000002E-02, 3.2416890264999998E-02, 4.0971460860000000E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 2.5266217610000002E-03, -1.3628904731000001E-02, 7.1087981962999994E-02, 1.7280032495600001E-01, 2.4879425055800000E-01, 2.4833015539900000E-01, -2.4461758511800000E-01, -4.5762595545000001E-01, -8.9205181260000002E-02, -2.5484913309999999E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 4.9534886579999998E-03, 4.1868496635000002E-02, 8.5218439293000006E-02, 3.9708847399399999E-01, 5.5350224331099995E-01, 6.5564701418200000E-01, 7.6861537470499997E-01, 9.7313473526500005E-01, 9.8548799995100000E-01, 9.9658645961799996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0056935882252931E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    19;
    constant Real[19] data_x =   {5.4000000000000000E+01, 5.6625000000000000E+01, 5.7375000000000000E+01, 5.7625000000000000E+01, 5.8125000000000000E+01, 5.8625000000000000E+01, 5.8875000000000000E+01, 5.9625000000000000E+01, 6.1375000000000000E+01, 6.1875000000000000E+01, 6.2125000000000000E+01, 6.2375000000000000E+01, 6.2625000000000000E+01, 6.2875000000000000E+01, 6.3125000000000000E+01, 6.3375000000000000E+01, 6.3625000000000000E+01, 6.3875000000000000E+01, 6.5000000000000000E+01};
    constant Real[19] data_y =   {0.0000000000000000E+00, 5.2274142811999998E-02, 1.5050775305799999E-01, 2.4434102027300000E-01, 2.1001284517800001E-01, 1.8836680447000001E-02, 3.7190894809999999E-03, 3.5878426280000000E-03, 2.7182756785000001E-02, 2.7583269274999999E-02, 6.1829242739999998E-02, 1.7303905102299999E-01, 5.3927327757200005E-01, 7.4747126108299999E-01, 6.1879242579299998E-01, 2.1179638721700000E-01, 1.1405819773000000E-02, 1.7119430090000000E-03, 0.0000000000000000E+00};
    constant Real[19] m_k =      {0.0000000000000000E+00, 2.4815860825000001E-02, 2.5464295843500001E-01, 2.8652067007299997E-01, -3.8743147409899997E-01, -1.6443650143600000E-01, -1.2251979901000000E-02, 1.3753196239000001E-02, -1.9745921864000000E-02, 4.8239742088999997E-02, 1.6598529155600000E-01, 7.7892690460900005E-01, 1.1925732183250000E+00, 3.4475869165200002E-01, -1.2047736249140000E+00, -1.1375777279779999E+00, -7.2961282219000001E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[19] iy_start = {0.0000000000000000E+00, 5.4573429038999997E-02, 1.2009966282900000E-01, 1.6948278797200000E-01, 2.9761283518300002E-01, 3.5038579843899997E-01, 3.5242059679700000E-01, 3.5394767365599999E-01, 3.8956042661899998E-01, 4.0188374186300002E-01, 4.1248850602400000E-01, 4.3875733377699999E-01, 5.2598296656499999E-01, 6.9189033225800001E-01, 8.7144572719400004E-01, 9.7532545725999997E-01, 9.9776859494000003E-01, 9.9903325265800003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0039247400262648E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.8000000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT64HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_RT64HC;