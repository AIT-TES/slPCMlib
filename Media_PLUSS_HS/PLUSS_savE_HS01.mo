within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS01 "Pluss Advanced Technologies Pvt Ltd, HS01; data taken from: PLUSS datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS01";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6914999999999998E+02, 2.7814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6814999999999998E+02, 2.7614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.9000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.3172344373282418E+05
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
    constant Integer  len_x =    16;
    constant Real[16] data_x =   {-4.0000000000000000E+00, -2.8750000000000000E+00, -1.3750000000000000E+00, -1.1250000000000000E+00, -8.7500000000000000E-01, -6.2500000000000000E-01, -3.7500000000000000E-01, -1.2500000000000000E-01, 1.2500000000000000E-01, 3.7500000000000000E-01, 6.2500000000000000E-01, 8.7500000000000000E-01, 1.3750000000000000E+00, 2.8750000000000000E+00, 4.3750000000000000E+00, 5.0000000000000000E+00};
    constant Real[16] data_y =   {0.0000000000000000E+00, 1.5273706720000000E-03, 1.3723994580000000E-03, 2.4589833160000001E-03, 1.1552180948000000E-02, 6.8491744051000003E-02, 4.2288030994900000E-01, 8.3492496100900004E-01, 9.7657545298500004E-01, 7.3440557922899996E-01, 3.0223321127600000E-01, 1.3450356964300000E-01, 5.1810786779000001E-02, 3.6182635177000003E-02, 2.2839336799999999E-04, 0.0000000000000000E+00};
    constant Real[16] m_k =      {0.0000000000000000E+00, 3.5138415139999999E-03, -8.5606041739999996E-03, 6.4357176630000000E-03, 4.1389393243999997E-02, 7.6734206724599996E-01, 1.5383604733880001E+00, 1.5069265068110000E+00, -3.3507952166900001E-01, -1.4186623967609999E+00, -1.0460079304080001E+00, -3.1288521268300001E-01, 1.7234921896000000E-02, -4.6405262645999999E-02, -7.7786636099999996E-04, 0.0000000000000000E+00};
    constant Real[16] iy_start = {0.0000000000000000E+00, 4.8937485900000003E-04, 4.9356960650000000E-03, 5.3371934740000003E-03, 6.9092026479999997E-03, 1.3144256113000000E-02, 7.0647490931999998E-02, 2.2830404401000001E-01, 4.6473605108999999E-01, 6.8462497981799997E-01, 8.1248059286300001E-01, 8.6334053374200004E-01, 9.0310901444500002E-01, 9.8116890116400002E-01, 9.9995387006799996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0016975455192136E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-5.0000000000000000E+00, -4.1250000000000000E+00, -2.3750000000000000E+00, -1.3750000000000000E+00, -8.7500000000000000E-01, -6.2500000000000000E-01, -3.7500000000000000E-01, -1.2500000000000000E-01, 1.2500000000000000E-01, 6.2500000000000000E-01, 8.7500000000000000E-01, 1.3750000000000000E+00, 3.0000000000000000E+00};
    constant Real[13] data_y =   {0.0000000000000000E+00, 3.2489959700000002E-04, 5.2800521474000002E-02, 3.4368771001999997E-02, 9.5555394141999997E-02, 2.6762327680100001E-01, 8.2022197370600003E-01, 1.1295876922129999E+00, 9.0103270673699998E-01, 1.8314076158999999E-02, 0.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.0378416150000000E-03, 6.0959878580000001E-02, -5.5308610975999997E-02, 2.5433292072000002E-01, 1.1924310432809999E+00, 1.8072917757160001E+00, 2.9429877186199999E-01, -1.8258105144540000E+00, -1.8318480656800001E-01, 0.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 7.6133571000000006E-05, 3.1353088858999997E-02, 8.4771633208000005E-02, 1.1087258830799999E-01, 1.5149414965200000E-01, 2.8463344574900001E-01, 5.3692397272100001E-01, 8.0251397002699998E-01, 9.9866119752700000E-01, 1.0000000000000000E+00, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0027191156185309E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 9.2400000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.0100000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 2.2000000000000002E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 5.5000000000000004E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS01</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
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
end PLUSS_savE_HS01;
