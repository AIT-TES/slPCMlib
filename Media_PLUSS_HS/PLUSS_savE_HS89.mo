
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS89 "Pluss Advanced Technologies Pvt Ltd, HS89; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS89";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.5614999999999998E+02, 3.6414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.5514999999999998E+02, 3.6614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.8000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.6500000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.6040298766699349E+05
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
    constant Real[12] data_x =   {8.3000000000000000E+01, 8.6625000000000000E+01, 8.7625000000000000E+01, 8.8375000000000000E+01, 8.8875000000000000E+01, 8.9375000000000000E+01, 8.9625000000000000E+01, 8.9875000000000000E+01, 9.0375000000000000E+01, 9.0625000000000000E+01, 9.0875000000000000E+01, 9.1000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.3729952958999999E-02, 8.5438388552999997E-02, 1.8574206185799999E-01, 2.2974809082600001E-01, 3.7807118561100000E-01, 5.2170875343099998E-01, 5.7057610319700003E-01, 3.5123198143599998E-01, 1.7181006609899999E-01, 3.8903867674000003E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 8.1081405070000008E-03, 1.5502605347100001E-01, 8.5734069435000004E-02, 1.6165574864599999E-01, 4.2327778448300002E-01, 4.5640604076199998E-01, -4.2195602931999997E-02, -6.8202915265099995E-01, -6.5839366956500001E-01, -4.4916426146100003E-01, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.6092779500000001E-02, 5.3634591820999999E-02, 1.5913964038500000E-01, 2.6198054502899998E-01, 4.0927273103900003E-01, 5.2217657246299998E-01, 6.6205724691599999E-01, 9.0715006876500004E-01, 9.7275814303899999E-01, 9.9814342726500005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0053774916977649E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    16;
    constant Real[16] data_x =   {8.2000000000000000E+01, 8.3375000000000000E+01, 8.4625000000000000E+01, 8.5625000000000000E+01, 8.6125000000000000E+01, 8.6375000000000000E+01, 8.6625000000000000E+01, 8.6875000000000000E+01, 8.7125000000000000E+01, 8.7375000000000000E+01, 8.7625000000000000E+01, 8.7875000000000000E+01, 8.8375000000000000E+01, 8.9375000000000000E+01, 9.1625000000000000E+01, 9.3000000000000000E+01};
    constant Real[16] data_y =   {0.0000000000000000E+00, 3.1101161039000000E-02, 5.4536366938000000E-02, 1.9919893941000001E-02, 6.3639297846999995E-02, 2.1491124328100000E-01, 7.9851044307499996E-01, 1.1391584371179999E+00, 8.8127475985899995E-01, 2.7231782097700002E-01, 2.4021370720000001E-03, 0.0000000000000000E+00, 0.0000000000000000E+00, 1.1567366783000001E-02, 1.4100399507999999E-02, 0.0000000000000000E+00};
    constant Real[16] m_k =      {0.0000000000000000E+00, 2.8161449630000001E-03, 3.6065114746000000E-02, -6.0679615531000000E-02, 1.8584937797900000E-01, 1.1535529753780001E+00, 1.9268359773770001E+00, 3.8341351326299999E-01, -1.9280130954700001E+00, -1.6661487544069999E+00, -1.7136639015000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00, -1.4682035431000001E-02, 1.0045919687999999E-02, 0.0000000000000000E+00};
    constant Real[16] iy_start = {0.0000000000000000E+00, 2.0855537718999999E-02, 6.9855113873000002E-02, 1.1496616130600000E-01, 1.3065762500700001E-01, 1.6031853109700001E-01, 2.8248358744099999E-01, 5.3174101128800000E-01, 7.9528725292299995E-01, 9.3755746882900004E-01, 9.6320700004100002E-01, 9.6341717918900005E-01, 9.6341717918900005E-01, 9.7039664883800003E-01, 9.8876782469600000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9604453287263095E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.6300000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.5400000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 6.9999999999999996E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.9000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS89</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
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
end PLUSS_savE_HS89;