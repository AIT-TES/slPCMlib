
// within slPCMlib.Media_Climator;
package ClimSel_C28 "Climator Sweden AB, ClimSel C28; data taken from: Climator Sweden AB datasheet; last access: 2022-10-14."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C28";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9514999999999998E+02, 3.0714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9414999999999998E+02, 3.0214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.3000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4357229203206539E+05
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
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {2.2000000000000000E+01, 2.4875000000000000E+01, 2.6375000000000000E+01, 2.7375000000000000E+01, 2.8375000000000000E+01, 2.9625000000000000E+01, 3.0375000000000000E+01, 3.0875000000000000E+01, 3.2375000000000000E+01, 3.4000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.5514712086000001E-02, 4.6758529022999998E-02, 1.0570105353299999E-01, 1.4174636774499999E-01, 3.1629129027000003E-01, 2.5883147988700000E-01, 1.4895694675900001E-01, 2.3141667700999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.4318617776000001E-02, 7.3886760171000002E-02, 5.7674592690000002E-03, 1.3256993244199999E-01, 1.2748024143200001E-01, -2.2153267925600001E-01, -1.7383708485300001E-01, -1.2420309697000000E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 2.6822419630999999E-02, 6.9870705465999994E-02, 1.5180067014400001E-01, 2.6499006395899999E-01, 5.5200887378399999E-01, 7.8410665004699998E-01, 8.8508914344000000E-01, 9.8392588672299996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0002877063473021E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {2.1000000000000000E+01, 2.2375000000000000E+01, 2.4875000000000000E+01, 2.6375000000000000E+01, 2.6625000000000000E+01, 2.6875000000000000E+01, 2.7125000000000000E+01, 2.7375000000000000E+01, 2.7625000000000000E+01, 2.7875000000000000E+01, 2.8875000000000000E+01, 2.9000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.6442215973999998E-02, 6.1571084259999999E-02, 3.9935850118900001E-01, 5.8018042651699997E-01, 6.4648888192200005E-01, 5.6839776734600000E-01, 3.8109995105700001E-01, 1.6893621521200000E-01, 7.3909782063999993E-02, 7.3761231900000004E-04, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -1.8747227203999999E-02, 5.9970456525999999E-02, 5.4718775616600002E-01, 5.8977667397400002E-01, -3.0415796257000000E-02, -6.8501929324900002E-01, -7.9718060893099996E-01, -6.9805112489300003E-01, -2.0528359198599999E-01, -7.6077043169999999E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.4324154636000000E-02, 7.1105460991000000E-02, 3.2663510630300002E-01, 4.4942542240400002E-01, 6.0671912392100003E-01, 7.6271318659900000E-01, 8.8254059394999995E-01, 9.5109692809099999E-01, 9.7901572804299997E-01, 9.9996363636200003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0046617972586853E+00;
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
    lambda := 9.7999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.5000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C28</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    macroencapsulation<br>  Data taken from: Climator Sweden AB datasheet - last access 2022-10-14.<br><br>
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
end ClimSel_C28;