within slPCMlib.Media_Rubitherm_RT;
package RT18HC "Rubitherm RT18HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT18HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.881500000000000e+02, 2.931500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.861500000000000e+02, 2.921500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.144326842067587e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 2.881500000000000e+02
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
    constant Integer pieces =  7;
    constant Integer[7] order =  {1, 6, 5, 6, 5, 6, 1};
    constant Real[8] breaks = {-8.500000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 1.200000000000000e+02};
    constant Real[7,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 5.809983120114349e-14, 1.122493530223062e-16, -4.090168532635013e-19, -1.156482317317871e-18, 3.662988226052195e-03, -1.725901837592933e-03},  {1.937086388517473e-03, 7.959530104756375e-03, 1.074135469662806e-02, 2.111845508663292e-03, -7.573586433633019e-03, 8.581878134817141e-03, 0.000000000000000e+00},  {2.375810839974932e-02, 4.839282096354158e-02, 5.745415396899123e-02, 5.763628112230264e-02, 3.533580424045268e-02, 1.582159764517500e-02, -1.610337599940753e-02},  {2.222953903408049e-01, 4.600409214596641e-01, 3.590431592392523e-01, 3.512795454771286e-02, -1.271068475247852e-01, 3.784233543727492e-02, 0.000000000000000e+00},  {9.872429134999240e-01, 9.642953906684589e-01, 8.020929210642880e-02, -9.487608117867875e-02, 6.210482966158937e-02, -2.122103937566787e-02, 2.933357814449998e-03},  {1.980688663196505e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  8;
    constant Integer[8] order =  {1, 6, 5, 5, 6, 5, 6, 1};
    constant Real[9] breaks = {-8.700000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 1.190000000000000e+02};
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 7.777526421142647e-16, -3.729217026907731e-19, 2.726392448195468e-18, 1.445602896647339e-19, 6.065615288230704e-04, -4.644598573607718e-05},  {5.601155430877736e-04, 2.754131729698285e-03, 5.368925502189555e-03, 5.136695573509163e-03, 2.336117858074194e-03, -1.332987803320291e-03, 0.000000000000000e+00},  {1.482299840323868e-02, 3.158160187029301e-02, 2.146584133795931e-02, 1.151288972603027e-03, -4.328821158527260e-03, 4.320856855538742e-03, 0.000000000000000e+00},  {6.901376628110550e-02, 8.225615110759660e-02, 4.215534985999225e-02, 2.704457289388143e-02, 1.727546311916645e-02, -1.812662239324916e-02, 7.737882031780979e-03},  {2.273565629002740e-01, 2.725966022102869e-01, 1.617438538008583e-01, 6.963784207367518e-02, 4.271058162963532e-02, -3.434505681058853e-02, 0.000000000000000e+00},  {7.397003858041413e-01, 8.041148784986262e-01, 2.834703016938105e-01, -1.029703995136688e-01, -1.290147024233073e-01, 1.341028817927465e-01, -3.609998043602836e-02},  {1.693303365416320e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT18HC</strong>.<br><br>
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
end RT18HC;
