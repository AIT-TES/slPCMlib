
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT90HC "Rubitherm GmbH, RT90HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT90HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.5814999999999998E+02, 3.6714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.5714999999999998E+02, 3.6514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.2707059932403779E+05
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
    constant Real[11] data_x =   {8.5000000000000000E+01, 8.6125000000000000E+01, 8.7125000000000000E+01, 8.8375000000000000E+01, 8.9875000000000000E+01, 9.0375000000000000E+01, 9.0625000000000000E+01, 9.1375000000000000E+01, 9.1625000000000000E+01, 9.3125000000000000E+01, 9.4000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.3983769993000000E-02, 7.2189291620000000E-03, 6.4584754113000006E-02, 1.1328410633500000E-01, 2.3518397635499999E-01, 3.6213773127299997E-01, 4.3697839334799998E-01, 3.3472055718900001E-01, 2.9783901889000000E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -1.7203592035999998E-02, 8.1674305929999996E-03, 8.4316636099999998E-04, 1.2638065849600000E-01, 3.9955736773000000E-01, 4.4393861012100000E-01, -3.6206645325699999E-01, -3.5207729444500002E-01, -6.8379815627000004E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.5258623884999999E-02, 2.8704579865000000E-02, 7.4395757048000000E-02, 1.8392398994799999E-01, 2.6510144375900002E-01, 3.3930844755599998E-01, 6.7572910731000002E-01, 7.7184535436799995E-01, 9.9135875794700001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9694955039123057E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {8.4000000000000000E+01, 8.5375000000000000E+01, 8.6875000000000000E+01, 8.7625000000000000E+01, 8.8375000000000000E+01, 8.8625000000000000E+01, 8.9125000000000000E+01, 8.9375000000000000E+01, 8.9625000000000000E+01, 8.9875000000000000E+01, 9.0125000000000000E+01, 9.0375000000000000E+01, 9.0625000000000000E+01, 9.1875000000000000E+01, 9.2000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 1.9815415680999999E-02, 1.0137062059000000E-02, 5.0681194826000003E-02, 6.3709236237000005E-02, 4.4268801583999998E-02, 8.0374429418999999E-02, 1.7275022745999999E-01, 4.1986922434300000E-01, 6.3657839273899997E-01, 7.1298803843000003E-01, 6.1030818734000003E-01, 3.8065952320699997E-01, 4.9803738990000001E-03, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, -2.7858782158000001E-02, 1.5563327183999999E-02, 7.9574742623000003E-02, -5.1429540017000001E-02, -5.1945262446000000E-02, 1.6465242404200001E-01, 7.8123723036000003E-01, 9.2467893317000005E-01, 7.8537608812399995E-01, -7.0405528100999995E-02, -7.6559094127399996E-01, -6.9981486065300003E-01, -5.3573865678000002E-02, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 1.8035707884999998E-02, 3.2377027840000003E-02, 5.2209070165000003E-02, 1.0131001360200000E-01, 1.1482749264300000E-01, 1.4151046843000001E-01, 1.6997660474099999E-01, 2.4340221009700000E-01, 3.7635619909200002E-01, 5.4973414819900002E-01, 7.1898653035500004E-01, 8.4267539058399998E-01, 9.9975817050700000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0012991359859884E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.5000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.5000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.5000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT90HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2022-12-06.<br><br>
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
end Rubitherm_RT90HC;