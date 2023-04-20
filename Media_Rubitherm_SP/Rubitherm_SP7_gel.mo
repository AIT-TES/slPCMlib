
// within slPCMlib.Rubitherm_SP;
package Rubitherm_SP7_gel "Rubitherm GmbH, SP7_gel; data taken from: Rubitherm datasheet; last access: 2020-07-12."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP7_gel";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7214999999999998E+02, 2.9014999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7114999999999998E+02, 2.9014999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.3700000000000000E+05
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
    constant Real[18] data_x =   {-1.0000000000000000E+00, 6.2500000000000000E-01, 2.6250000000000000E+00, 3.8750000000000000E+00, 6.1250000000000000E+00, 7.1250000000000000E+00, 7.3750000000000000E+00, 7.6250000000000000E+00, 8.1250000000000000E+00, 8.3750000000000000E+00, 8.6250000000000000E+00, 9.1250000000000000E+00, 9.8750000000000000E+00, 1.1375000000000000E+01, 1.2625000000000000E+01, 1.4625000000000000E+01, 1.6625000000000000E+01, 1.7000000000000000E+01};
    constant Real[18] data_y =   {0.0000000000000000E+00, 2.5016760030000001E-02, 6.6201257720000004E-02, 1.0442606425100000E-01, 3.3391535264000002E-02, 9.2370920218000002E-02, 1.5776220300300001E-01, 2.8204961855200000E-01, 3.3172156940500003E-01, 2.3089485882399999E-01, 1.0311480985500000E-01, 2.8886698658000001E-02, 2.0120260530000000E-02, 0.0000000000000000E+00, 2.7718672848999999E-02, 5.1533654608000003E-02, 3.6663870260000000E-03, 0.0000000000000000E+00};
    constant Real[18] m_k =      {0.0000000000000000E+00, 3.6380730100000000E-04, 4.5589219309999998E-03, 2.1103506261999999E-02, -4.5802488999999997E-05, 1.6419212974000000E-01, 3.4838929222999998E-01, 4.1243419753999999E-01, -3.3009718117999998E-01, -4.5908102835499998E-01, -3.7391215186200000E-01, -3.9612049517000002E-02, -3.4975825420000002E-02, 0.0000000000000000E+00, 4.1047198970000000E-03, 2.6926237382000000E-02, -1.5143016497999999E-02, 0.0000000000000000E+00};
    constant Real[18] iy_start = {0.0000000000000000E+00, 2.0149733667999999E-02, 1.0954203333199999E-01, 2.1353273153999999E-01, 3.7671976785099998E-01, 4.2568044109699998E-01, 4.5584352441300002E-01, 5.1022645325299998E-01, 6.7833499855199997E-01, 7.4899604257100005E-01, 7.9010712865099997E-01, 8.1601904646099999E-01, 8.3409292960899994E-01, 8.4258456285500005E-01, 8.5929438211400000E-01, 9.3059866166000005E-01, 9.9949243610399996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9524217072484900E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    17;
    constant Real[17] data_x =   {-2.0000000000000000E+00, -1.2500000000000000E-01, 1.6250000000000000E+00, 3.6250000000000000E+00, 4.6250000000000000E+00, 5.3750000000000000E+00, 5.6250000000000000E+00, 6.1250000000000000E+00, 6.3750000000000000E+00, 6.6250000000000000E+00, 6.8750000000000000E+00, 8.3750000000000000E+00, 1.0625000000000000E+01, 1.2625000000000000E+01, 1.4375000000000000E+01, 1.5625000000000000E+01, 1.7000000000000000E+01};
    constant Real[17] data_y =   {0.0000000000000000E+00, 2.9518880765000002E-02, 2.9361896108000000E-02, 5.4023842138000003E-02, 6.0628789952000002E-02, 1.6164618803899999E-01, 2.8238524244199997E-01, 3.2551026958899998E-01, 2.2544598516000000E-01, 1.0034969930700000E-01, 5.0404610433999998E-02, 4.5967104629000000E-02, 3.3890873687000000E-02, 3.1311064558999999E-02, 6.2571120521000001E-02, 4.8203232414000000E-02, 0.0000000000000000E+00};
    constant Real[17] m_k =      {0.0000000000000000E+00, 2.7122660840999999E-02, 2.0034077730000000E-03, 3.0885026199999999E-02, -4.8363813149999999E-03, 3.3255979352699999E-01, 3.9535590294200001E-01, -3.3259419206399998E-01, -4.5096126411499998E-01, -3.6750483233400000E-01, -9.5127644082000001E-02, -2.6839327242000000E-02, 2.4603218164999999E-02, 1.1105296638000001E-02, -1.6141838158000001E-02, -8.0613910649999995E-03, 0.0000000000000000E+00};
    constant Real[17] iy_start = {0.0000000000000000E+00, 1.9431675720000000E-02, 7.6493249962000007E-02, 1.4914441302600001E-01, 2.0854215626600001E-01, 2.7506585478000001E-01, 3.2941432614900001E-01, 4.9404449278700002E-01, 5.6248729574199996E-01, 6.0217219915599995E-01, 6.1933623768099999E-01, 6.7791804453399995E-01, 7.4503296223299997E-01, 8.1368775243100000E-01, 9.0145065391599999E-01, 9.6860883938900000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.8498656347846314E-01;
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
    lambda := 1.3000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP7_gel</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    macroencapsulation<br>  Data taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP7_gel;