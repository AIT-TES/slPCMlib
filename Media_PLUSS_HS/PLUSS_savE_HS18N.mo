
// within slPCMlib.PLUSS_HS;
package PLUSS_savE_HS18N "Pluss Advanced Technolgies Pvt Ltd, HS18N; data taken from: PLUSS datasheet; last access: 2022-02-13."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS18N";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.5014999999999998E+02, 2.6114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4914999999999998E+02, 2.5814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.9000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.4800000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.2338980983476364E+05
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
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-2.3000000000000000E+01, -2.1375000000000000E+01, -1.8625000000000000E+01, -1.7625000000000000E+01, -1.7375000000000000E+01, -1.6625000000000000E+01, -1.6375000000000000E+01, -1.5125000000000000E+01, -1.3625000000000000E+01, -1.2625000000000000E+01, -1.2000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 4.2603277799000000E-02, 1.4768645066800001E-01, 2.8247137588600002E-01, 3.3708850014600000E-01, 2.5552529512700001E-01, 1.6957139344300001E-01, 5.0363504167999999E-02, 1.2023255266000000E-02, 2.4489433099999998E-04, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.2982366218000000E-02, 8.5918640924999998E-02, 1.8520682987000001E-01, 1.8853111630000000E-01, -3.0602989186700003E-01, -2.7911345652300001E-01, -2.5412538011000001E-02, -3.6055452597999998E-02, -8.2909311500000001E-04, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 3.2191571065999999E-02, 2.5081695919699998E-01, 4.6044280881999999E-01, 5.3892664266099999E-01, 7.8768695082600004E-01, 8.4140675895299999E-01, 9.4725652510000002E-01, 9.9670760840700001E-01, 9.9994978346200003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0136406445168094E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {-2.4000000000000000E+01, -2.2375000000000000E+01, -2.1125000000000000E+01, -2.0875000000000000E+01, -2.0625000000000000E+01, -2.0375000000000000E+01, -2.0125000000000000E+01, -1.9875000000000000E+01, -1.9625000000000000E+01, -1.9375000000000000E+01, -1.9125000000000000E+01, -1.8875000000000000E+01, -1.7625000000000000E+01, -1.5000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 3.6509141910000001E-02, 8.4589589240000006E-03, 2.8645455140000001E-02, 1.3499681459000001E-01, 6.8454047980199995E-01, 1.1619941019290001E+00, 1.0873062867620000E+00, 5.3583440898300005E-01, 5.9802884579000001E-02, 6.9735154310000000E-03, 1.1588479050000000E-03, 1.4719888617000000E-02, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, 3.0776642486000001E-02, 1.6412008113000000E-02, 9.4556170847999998E-02, 1.3906230174809999E+00, 2.0697765124270000E+00, 1.6197377689790000E+00, -1.9766344722180000E+00, -2.0460239156390001E+00, -7.1763461494900005E-01, -2.5758275476999999E-02, -2.6239553620000002E-03, -1.5422011224000000E-02, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 2.2862717423999999E-02, 5.2800860657000001E-02, 5.7026644552000003E-02, 7.0714519547000004E-02, 1.6949630198800000E-01, 4.0236682676899999E-01, 7.0188715876200003E-01, 9.0488812875799995E-01, 9.7234002416400001E-01, 9.7707764659399998E-01, 9.7797258530900000E-01, 9.8954878059399998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9875515199836429E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.0830000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.0950000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.0000000000000000E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.0950000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS18N</strong>  from manufacturer: <strong>Pluss Advanced Technolgies Pvt Ltd</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    multiple options available<br>  Data taken from: PLUSS datasheet - last access 2022-02-13.<br><br>
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
end PLUSS_savE_HS18N;