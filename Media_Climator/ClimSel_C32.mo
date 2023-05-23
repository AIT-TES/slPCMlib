
within slPCMlib.Media_Climator;
package ClimSel_C32 "Climator Sweden AB, ClimSel C32; data taken from: Climator Sweden AB datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C32";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9714999999999998E+02, 3.0814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9614999999999998E+02, 3.0314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {4.9000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.5000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 8.9264022167714240E+04
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
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {2.4000000000000000E+01, 2.5375000000000000E+01, 2.6625000000000000E+01, 2.8625000000000000E+01, 2.9375000000000000E+01, 3.0125000000000000E+01, 3.0875000000000000E+01, 3.1375000000000000E+01, 3.1625000000000000E+01, 3.2125000000000000E+01, 3.2375000000000000E+01, 3.2625000000000000E+01, 3.4125000000000000E+01, 3.5000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 3.3473936257000000E-02, 3.5211657420000000E-02, 5.9412361455000000E-02, 9.2317621284000007E-02, 1.3790735284000000E-01, 9.4157072406000006E-02, 1.7707645482600001E-01, 3.0386443425600002E-01, 4.0715584304000002E-01, 3.4261840247399999E-01, 2.2429873580400000E-01, 1.3115782104000001E-02, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, -2.4314947666000002E-02, 3.5803580054999999E-02, 4.8600714099999998E-04, 8.7466933070999994E-02, -3.1725626120000003E-02, 3.4347674150999999E-02, 3.6986098544599999E-01, 4.3391443362499998E-01, -1.2391318844500000E-01, -4.0095326332400000E-01, -3.6556117038199998E-01, -3.0876184970000001E-02, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 2.7063116504000000E-02, 6.2449924605999997E-02, 1.6971413238800001E-01, 2.2296640484200000E-01, 3.1563754200200000E-01, 4.0024894061400001E-01, 4.6156343726299998E-01, 5.2183497564900005E-01, 7.1275582035700003E-01, 8.0869658692199997E-01, 8.7995329383099996E-01, 9.9620108064099999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0081550074509873E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {2.3000000000000000E+01, 2.4625000000000000E+01, 2.5625000000000000E+01, 2.6875000000000000E+01, 2.7375000000000000E+01, 2.7625000000000000E+01, 2.8375000000000000E+01, 2.9625000000000000E+01, 3.0000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 3.3941744798000002E-02, 3.2784479637000000E-02, 7.7994661748000002E-02, 2.0941916667499999E-01, 3.5715053191899998E-01, 5.3966587081499995E-01, 6.1557579672999997E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 4.6183180286999997E-02, -2.6541360741000001E-02, 1.2745058569300000E-01, 4.9398264255199997E-01, 5.4428413870600001E-01, -2.8040383276699998E-01, -2.5793418169799998E-01, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 1.7465215619999998E-02, 5.7002467570999998E-02, 1.0633032827000000E-01, 1.7073300795700000E-01, 2.4149584061500001E-01, 6.1754124314799996E-01, 9.9145603628800005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0028856225595226E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.5000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.5000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 7.6000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.0800000000000001E+00;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C32</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
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
end ClimSel_C32;