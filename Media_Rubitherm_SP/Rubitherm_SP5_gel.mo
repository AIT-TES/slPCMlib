within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP5_gel "Rubitherm GmbH, SP5_gel; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP5_gel";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6814999999999998E+02, 2.8714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6814999999999998E+02, 2.8714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4500000000000000E+05
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
    constant Real[14] data_x =   {-5.0000000000000000E+00, -2.3750000000000000E+00, -3.7500000000000000E-01, 2.1250000000000000E+00, 3.8750000000000000E+00, 5.1250000000000000E+00, 5.8750000000000000E+00, 6.3750000000000000E+00, 6.8750000000000000E+00, 7.8750000000000000E+00, 9.6250000000000000E+00, 1.1625000000000000E+01, 1.3625000000000000E+01, 1.4000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 8.0216410070000008E-03, 1.8767540215000000E-02, 5.1416679064000000E-02, 1.2792018695700000E-01, 1.9002843942200001E-01, 1.2751634921200000E-01, 1.3036730899400001E-01, 1.5178544132999999E-01, 5.7657969184999998E-02, 2.9587527924999999E-02, 3.5875476516000003E-02, 2.2371957270000001E-03, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, -1.1337433204000000E-02, 5.4963075889999996E-03, 2.3560299034000001E-02, 7.5044492148000005E-02, -5.2117764501000002E-02, -3.4902793175000002E-02, 4.5877996701999997E-02, -2.2494931154000000E-02, -5.9500585389000001E-02, -7.8538135920000005E-03, 1.5144080067000000E-02, -9.2186359710000002E-03, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 1.6971970281000000E-02, 3.8067124389000000E-02, 1.1608292525900001E-01, 2.5930147930699998E-01, 4.7373548879900002E-01, 5.9154553052600001E-01, 6.5408808666899998E-01, 7.2576942317500004E-01, 8.3315354342100001E-01, 8.9606579083700000E-01, 9.5363691400200001E-01, 9.9968977405899995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9609120598849887E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {-5.0000000000000000E+00, -3.3750000000000000E+00, -1.6250000000000000E+00, 6.2500000000000000E-01, 2.6250000000000000E+00, 3.3750000000000000E+00, 3.6250000000000000E+00, 4.1250000000000000E+00, 4.3750000000000000E+00, 4.6250000000000000E+00, 5.6250000000000000E+00, 7.8750000000000000E+00, 9.8750000000000000E+00, 1.1625000000000000E+01, 1.4000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 5.3582439436000003E-02, 4.1859494763000001E-02, 4.5961700406999997E-02, 6.5973174984999994E-02, 1.5197910411800000E-01, 2.3963891492799999E-01, 2.7940627454000000E-01, 2.1496727862900000E-01, 1.2244615828900000E-01, 4.5065770668999999E-02, 2.4344265819000000E-02, 3.1534295265999997E-02, 3.2161958364000001E-02, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, 5.5230983623000003E-02, 2.5799942880000001E-03, -1.5768983412999999E-02, 6.1016634649999997E-03, 2.5362260496799999E-01, 2.8971138677800001E-01, -1.8684129191500001E-01, -3.1921014512599999E-01, -2.7259299604400000E-01, -2.1389908310000001E-03, 5.1930926359999999E-03, 1.7571488620000001E-03, -2.1838866772999999E-02, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 3.1845821671000002E-02, 1.3022725097099999E-01, 2.3834158001100000E-01, 3.4453274147200003E-01, 4.1569871881100001E-01, 4.6518367678200001E-01, 6.0693757004700000E-01, 6.7034714268200002E-01, 7.1290074876400000E-01, 7.7502359781600005E-01, 8.5112495880200001E-01, 9.0899156759699995E-01, 9.7166035844300003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0147786012134070E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.4000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.3000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 5.9999999999999998E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP5_gel</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_SP5_gel;
