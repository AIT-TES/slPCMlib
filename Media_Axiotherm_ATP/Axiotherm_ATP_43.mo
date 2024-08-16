within slPCMlib.Media_Axiotherm_ATP;
package Axiotherm_ATP_43 "Axiotherm GmbH, ATP 43; data taken from: Axiotherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 43";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.0914999999999998E+02, 3.1814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.0814999999999998E+02, 3.1714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.3100000000000000E+05
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
    constant Integer  len_x =    18;
    constant Real[18] data_x =   {3.6000000000000000E+01, 3.9625000000000000E+01, 4.0125000000000000E+01, 4.0375000000000000E+01, 4.0625000000000000E+01, 4.1125000000000000E+01, 4.1375000000000000E+01, 4.1625000000000000E+01, 4.1875000000000000E+01, 4.2125000000000000E+01, 4.2375000000000000E+01, 4.2625000000000000E+01, 4.2875000000000000E+01, 4.3125000000000000E+01, 4.3625000000000000E+01, 4.3875000000000000E+01, 4.4625000000000000E+01, 4.5000000000000000E+01};
    constant Real[18] data_y =   {0.0000000000000000E+00, 5.4944449330000004E-03, 2.7143140913000002E-02, 7.7982662689000004E-02, 2.3323120149099999E-01, 3.6773005681299997E-01, 2.5130417620200002E-01, 9.3501145192999999E-02, 5.6556633763999997E-02, 7.3499632829999995E-02, 1.7471098648100000E-01, 5.4685761223999996E-01, 7.6838343605899995E-01, 6.3206703432500000E-01, 1.5439675466999999E-02, 1.6397592120000000E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[18] m_k =      {0.0000000000000000E+00, 3.0095822500000001E-03, 7.9763591876999995E-02, 4.3097701791199999E-01, 5.6742525407300004E-01, -3.3722620034099998E-01, -5.4961956277799995E-01, -4.2123788617500002E-01, -2.0778309213000001E-02, 1.2369186500600000E-01, 7.9364615649500003E-01, 1.2289812333770000E+00, 3.3966687321000000E-01, -1.2411959672249999E+00, -1.2572145802599999E-01, -5.9746411740000003E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[18] iy_start = {0.0000000000000000E+00, 6.7087651079999998E-03, 1.3314147936000000E-02, 2.4703275334000002E-02, 6.3156472052000007E-02, 2.3340425401600001E-01, 3.1242844859600000E-01, 3.5515169932800000E-01, 3.7193762269399999E-01, 3.8754862583900002E-01, 4.1527461326600001E-01, 5.0380683319499997E-01, 6.7400402709100005E-01, 8.5855204185900003E-01, 9.9814123257200005E-01, 9.9966285330500004E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0068636889762204E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    18;
    constant Real[18] data_x =   {3.5000000000000000E+01, 3.7125000000000000E+01, 3.8625000000000000E+01, 3.9375000000000000E+01, 3.9625000000000000E+01, 4.0125000000000000E+01, 4.0375000000000000E+01, 4.0625000000000000E+01, 4.0875000000000000E+01, 4.1875000000000000E+01, 4.2125000000000000E+01, 4.2375000000000000E+01, 4.2625000000000000E+01, 4.3125000000000000E+01, 4.3375000000000000E+01, 4.3625000000000000E+01, 4.3875000000000000E+01, 4.4000000000000000E+01};
    constant Real[18] data_y =   {0.0000000000000000E+00, 9.0267386160000004E-03, 2.5909088159000001E-02, 1.1824860338400001E-01, 2.4665708062899999E-01, 2.8892293747999997E-01, 1.7479725417799999E-01, 5.2733994240999997E-02, 1.7748725869999999E-02, 2.7722845140000001E-02, 6.3984968101999995E-02, 1.5281310931600001E-01, 3.7023600243499999E-01, 5.8117338904200000E-01, 4.5516298764099999E-01, 2.3967363925400001E-01, 5.6415778912000003E-02, 0.0000000000000000E+00};
    constant Real[18] m_k =      {0.0000000000000000E+00, 8.6186398809999995E-03, 4.8817366789999997E-03, 3.3482163347700000E-01, 4.2371988562500001E-01, -4.3279906911600002E-01, -4.7203376687199999E-01, -4.0535705344700002E-01, -4.7390669250000003E-02, 7.1424725088999996E-02, 1.7520221807299999E-01, 6.4308343074100005E-01, 7.9074913242699996E-01, -2.5427890067400000E-01, -8.0223098287399996E-01, -7.8981368046699996E-01, -6.7043949690600002E-01, 0.0000000000000000E+00};
    constant Real[18] iy_start = {0.0000000000000000E+00, 6.3935605169999996E-03, 3.3490474040000003E-02, 7.2362516242000005E-02, 1.1783893064400000E-01, 2.7067441371299999E-01, 3.2926406577599998E-01, 3.5756118076700000E-01, 3.6455729772500001E-01, 3.7748453119699998E-01, 3.8848641984400001E-01, 4.1332749047500000E-01, 4.7840636707599998E-01, 7.3990594531200005E-01, 8.7325848629799996E-01, 9.6067545873799998E-01, 9.9732781353699995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0072251205818217E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.4400000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.6000000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATP 43</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATP_43;
