within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP_minus_11 "Rubitherm GmbH, SP-11; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP-11";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.5514999999999998E+02, 2.6814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.5414999999999998E+02, 2.6314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.6649015346177510E+05
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
    constant Real[13] data_x =   {-1.8000000000000000E+01, -1.5875000000000000E+01, -1.3375000000000000E+01, -1.2125000000000000E+01, -1.1625000000000000E+01, -1.1375000000000000E+01, -1.1125000000000000E+01, -1.0875000000000000E+01, -1.0375000000000000E+01, -1.0125000000000000E+01, -9.1250000000000000E+00, -7.6250000000000000E+00, -5.0000000000000000E+00};
    constant Real[13] data_y =   {0.0000000000000000E+00, 3.5955819179999998E-03, 3.5378113560000000E-03, 1.0883609786600000E-01, 3.7630802454099999E-01, 6.7648998335499999E-01, 8.1500131583799995E-01, 7.0706543049500004E-01, 1.2613603938500001E-01, 4.1648653797000001E-02, 3.0522841435000000E-02, 2.4135891658000001E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 2.5507709959999999E-03, 2.7951005040000002E-03, 2.6046521862499999E-01, 8.4495221667499998E-01, 9.7791796925600005E-01, 8.8108650997999999E-02, -1.1408800585840000E+00, -1.1155010157190000E+00, -1.1558013272800000E-01, 2.7170746932000001E-02, -2.9910230200000001E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.8592550060000002E-03, 1.1645084004000000E-02, 4.8312715474999997E-02, 1.5737652908300001E-01, 2.8822927458900000E-01, 4.7922059861100003E-01, 6.7579808726900004E-01, 8.8348326292500001E-01, 8.9924186781700000E-01, 9.2342164247400005E-01, 9.7005161156099995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9958386161069734E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-1.9000000000000000E+01, -1.6875000000000000E+01, -1.5375000000000000E+01, -1.3125000000000000E+01, -1.2375000000000000E+01, -1.1625000000000000E+01, -1.1375000000000000E+01, -1.1125000000000000E+01, -1.0375000000000000E+01, -1.0125000000000000E+01, -1.0000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.5188824157000000E-02, 3.1841902720000002E-02, 4.5088870824000003E-02, 1.6052020873699999E-01, 4.3242460127400001E-01, 5.5794561028599998E-01, 5.8487314407699997E-01, 1.6591307682000001E-01, 3.7325060088000003E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.8988607288000000E-02, -1.5079258680000001E-02, 5.7175167706000000E-02, 2.7396593837599997E-01, 4.3967364351700000E-01, 4.4531065979300000E-01, -1.3529433786100001E-01, -6.6141377748999997E-01, -4.2902745615400001E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 9.0829132400000005E-03, 5.1161742514999997E-02, 1.0778912787700000E-01, 1.7540221108800000E-01, 3.9214231036699998E-01, 5.1715120882099996E-01, 6.6449138340299996E-01, 9.7377080066999999E-01, 9.9820800956099998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0100348420952103E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.1000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.2000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP-11</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_SP_minus_11;
