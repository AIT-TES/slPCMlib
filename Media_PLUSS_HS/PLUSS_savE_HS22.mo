
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS22 "Pluss Advanced Technolgies Pvt Ltd, HS22; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS22";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9114999999999998E+02, 3.0114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9014999999999998E+02, 2.9614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.2700000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.5300000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.6516110065422437E+05
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
    constant Real[10] data_x =   {1.8000000000000000E+01, 1.9875000000000000E+01, 2.1125000000000000E+01, 2.2125000000000000E+01, 2.3125000000000000E+01, 2.3875000000000000E+01, 2.5125000000000000E+01, 2.6375000000000000E+01, 2.7125000000000000E+01, 2.8000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.9396517798999999E-02, 8.0602779481000003E-02, 2.4190688209600000E-01, 2.5053906821900002E-01, 1.2893380807499999E-01, 1.1143379842100000E-01, 8.7183856998999998E-02, 6.2669858918999996E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.4949951981000000E-02, 9.4238908590000003E-02, 1.3899155605200000E-01, -1.4254882734200000E-01, -8.5372186050000007E-02, -1.7761941228000001E-02, 2.1361780440000001E-03, -8.9235250905999999E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 2.3328401270000000E-02, 8.2129532945000000E-02, 2.4066780558100001E-01, 5.1208645000099995E-01, 6.5260634141899998E-01, 7.9494199614100003E-01, 9.1726861915799995E-01, 9.7813565022100002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0064296225003264E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {1.7000000000000000E+01, 1.7625000000000000E+01, 1.8375000000000000E+01, 1.8875000000000000E+01, 2.0125000000000000E+01, 2.0875000000000000E+01, 2.1375000000000000E+01, 2.1625000000000000E+01, 2.2125000000000000E+01, 2.2375000000000000E+01, 2.2625000000000000E+01, 2.2875000000000000E+01, 2.3000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 9.5950069885000003E-02, 1.7372019023200000E-01, 1.2875541447700001E-01, 1.3576545273400001E-01, 1.2350745919000000E-01, 2.2858446750899999E-01, 3.6451982371899999E-01, 4.0663042927600002E-01, 2.9353509175300002E-01, 1.4697732740500000E-01, 3.3695073334000000E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.9703158917700000E-01, -7.9853268717000006E-02, -4.4307337322999997E-02, -2.0096249994999999E-02, 7.1289538107000006E-02, 3.7949887381100000E-01, 4.3901935358600003E-01, -3.6562059775500000E-01, -5.4863048649099999E-01, -5.3010424177799997E-01, -3.9244960136000001E-01, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.3581093827999999E-02, 1.3773718035999999E-01, 2.1264885193700001E-01, 3.7489405346100002E-01, 4.6787903382899998E-01, 5.4951728939400002E-01, 6.2337816668299995E-01, 8.3302228375999998E-01, 9.2153550917200000E-01, 9.7652752376100005E-01, 9.9840435047099996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0004448633283596E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.6510000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.5400000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 1.1299999999999999E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.5400000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS22</strong>  from manufacturer: <strong>Pluss Advanced Technolgies Pvt Ltd</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: multiple options available<br>  The data is taken from: PLUSS datasheet - last access 2022-02-13.<br><br>
  <br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-05-15.<br><br>
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
    <li>file creation date: 2023-05-15 </ul>
    </p></html>"));
end PLUSS_savE_HS22;