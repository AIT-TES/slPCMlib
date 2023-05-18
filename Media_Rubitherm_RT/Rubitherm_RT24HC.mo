
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT24HC "Rubitherm GmbH, RT24HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT24HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8814999999999998E+02, 3.0114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8814999999999998E+02, 2.9914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7781810081823650E+05
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
    constant Integer  len_x =    7;
    constant Real[7] data_x =   {1.5000000000000000E+01, 1.9125000000000000E+01, 2.2875000000000000E+01, 2.4625000000000000E+01, 2.6125000000000000E+01, 2.7375000000000000E+01, 2.8000000000000000E+01};
    constant Real[7] data_y =   {0.0000000000000000E+00, 2.8898207927999998E-02, 1.1098119017300000E-01, 2.6276772100800000E-01, 1.3120892307300000E-01, 1.0804014139000000E-02, 0.0000000000000000E+00};
    constant Real[7] m_k =      {0.0000000000000000E+00, 8.5713260440000005E-03, 5.5650566429000001E-02, 1.0240110522199999E-01, -1.5664209546600000E-01, -3.1995104374000001E-02, 0.0000000000000000E+00};
    constant Real[7] iy_start = {0.0000000000000000E+00, 4.7997456585999999E-02, 2.5749562431099998E-01, 5.7623912367399999E-01, 9.2427140707900002E-01, 9.9763824999899997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0115656577720546E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {1.5000000000000000E+01, 1.6875000000000000E+01, 2.1125000000000000E+01, 2.2375000000000000E+01, 2.2625000000000000E+01, 2.3125000000000000E+01, 2.3375000000000000E+01, 2.3625000000000000E+01, 2.3875000000000000E+01, 2.4125000000000000E+01, 2.4625000000000000E+01, 2.4875000000000000E+01, 2.6000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.5003945454000000E-02, 7.9851635034999993E-02, 9.5161276148000004E-02, 8.2700252545000005E-02, 1.5545017804500000E-01, 2.9045276559900002E-01, 5.9201983240000000E-01, 7.0408955963400000E-01, 5.3903705171700000E-01, 2.4632992789999999E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.1131815240000001E-03, 3.8882274219000000E-02, -3.3641352302000001E-02, -3.7383588609999997E-02, 2.8704892048899999E-01, 7.3519043714800003E-01, 9.2171139714799999E-01, -1.8900933776699999E-01, -1.0765445054980001E+00, -2.4941320800700001E-01, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.3139158072999998E-02, 1.8927985700800001E-01, 3.0822992761099999E-01, 3.3050529696399999E-01, 3.8333889091300000E-01, 4.3679833685600000E-01, 5.4624988139899999E-01, 7.1422341057100003E-01, 8.7440355483800003E-01, 9.9821804811100001E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0010420305027024E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.0000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.0000000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT24HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2022-11-18.<br><br>
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
end Rubitherm_RT24HC;