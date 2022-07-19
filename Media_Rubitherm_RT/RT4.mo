within slPCMlib.Media_Rubitherm_RT;
package RT4 "Rubitherm RT4; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT4";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.721500000000000e+02, 2.791500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.701500000000000e+02, 2.781500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.416943141660181e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 2.721500000000000e+02
             "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy  href = 0.0
             "reference enthalpy at Tref";

  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  protected
    constant Integer pieces =   data_H.pieces;
    constant Integer order[:] = data_H.order;
    constant Real breaks[:] =   data_H.breaks;
    constant Real coefs[:,:] =  data_H.coefs;
  algorithm
    (xi, dxi) := BasicUtilities.quartQuintSplineEval(T-273.15,
                 pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
  protected
    constant Integer pieces =   data_C.pieces;
    constant Integer order[:] = data_C.order;
    constant Real breaks[:] =   data_C.breaks;
    constant Real coefs[:,:] =  data_C.coefs;
  algorithm
    (xi, dxi) := BasicUtilities.quartQuintSplineEval(T-273.15,
                     pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complSolidification;

  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  9;
    constant Integer[9] order =  {1, 6, 6, 5, 5, 5, 5, 6, 1};
    constant Real[10] breaks = {-1.010000000000000e+02, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 1.060000000000000e+02};
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.769414262105793e-13, 1.159557840248382e-15, 1.937590912252106e-18, 0.000000000000000e+00, 1.798386285948729e-03, -4.905369520501085e-04},  {1.307849334076724e-03, 6.048709717632087e-03, 1.062580857873689e-02, 8.173123818485124e-03, 1.633877148992017e-03, -2.789648103431268e-03, 1.002465140053756e-03},  {2.600218563454533e-02, 5.042175724969700e-02, 3.208893899463805e-02, 6.861454181215647e-03, 2.722613732642028e-03, 6.071371101404251e-04, 0.000000000000000e+00},  {1.187040869028785e-01, 1.491101382638977e-01, 7.508035503554147e-02, 2.382328021318800e-02, 5.758299283344153e-03, -1.625144396401241e-03, 0.000000000000000e+00},  {3.708510153024486e-01, 3.856481641260332e-01, 1.648485474111588e-01, 3.060503338255221e-02, -2.367422698662053e-03, -5.391463004347319e-03, 0.000000000000000e+00},  {9.441938745191834e-01, 7.707333532795596e-01, 1.885444813233695e-01, -3.277928745556918e-02, -2.932473772039865e-02, 1.133946341972913e-02, 0.000000000000000e+00},  {1.852707147365873e+00, 9.888828197769297e-01, 2.765282683156318e-02, -3.668360413987248e-02, 2.737257937824700e-02, -1.089298226063585e-02, 1.806155461662151e-03},  {2.850844942413767e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  10;
    constant Integer[10] order =  {1, 6, 5, 5, 5, 6, 5, 5, 6, 1};
    constant Real[11] breaks = {-1.030000000000000e+02, -3.000000000000000e+00, -2.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 1.050000000000000e+02};
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.130690243999148e-29, 4.013172009700235e-32, 1.280627570230507e-18, 0.000000000000000e+00, 1.785161105591671e-03, -7.179585833215242e-04},  {1.067202522270148e-03, 4.618054028019708e-03, 7.082232306093878e-03, 3.492439389486223e-03, -1.843573221864510e-03, 4.064931216040647e-04, 0.000000000000000e+00},  {1.482284814560951e-02, 2.391800952925887e-02, 1.056304235940603e-02, 1.830777180688304e-04, 1.888923861558133e-04, 5.326968003916111e-04, 0.000000000000000e+00},  {5.020856693889067e-02, 4.901238094885936e-02, 1.757259783446351e-02, 6.265615266608196e-03, 2.852376388113869e-03, -1.482088550594492e-03, 0.000000000000000e+00},  {1.244294488263411e-01, 1.069534852171044e-01, 3.866281645702637e-02, 2.854235313118766e-03, -4.558066364858592e-03, 1.756660261531779e-02, -6.700036284851205e-03},  {2.792084857791987e-01, 2.222423539785575e-01, 9.504260608764098e-02, 2.628727030983823e-02, -1.722559756103770e-02, 8.042800741961335e-03, 0.000000000000000e+00},  {6.135979193361590e-01, 4.625009905493909e-01, 1.509788390705415e-01, 3.781288748530078e-02, 2.298840614876898e-02, -2.066620375627976e-02, 0.000000000000000e+00},  {1.267212838833882e+00, 8.665199369597134e-01, 1.956859008562613e-01, -7.689552548242087e-02, -8.034261263262980e-02, 8.734274775083009e-02, -2.375807507476804e-02},  {2.235765211210867e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.800000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.700000000000000e+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 2.000000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 2.000000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------

annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT4</strong>.<br><br>
  Information taken from: data_sheet - last access 01.12.2019.<br><br>
  It also contains the phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  These functions are modelled by piece-wise splines using <strong>variable order quartic and quintic</strong> method,
  see also 
  <blockquote>
  <p>
  Barz, T., Krämer, J., & Emhofer, J. (2020). Identification of Phase
  Fraction–Temperature Curves from Heat Capacity Data for Numerical
  Modeling of Heat Transfer in Commercial Paraffin Waxes.
  Energies, 13(19), 5149.
  <a href>doi.org/10.3390/en13195149</a>.
  </p>
  </blockquote>
  </p></html>",
  revisions="<html>
  <ul>
  <li>file creation date: 19-Jul-2022  </ul>
  </html>"));
end RT4;
