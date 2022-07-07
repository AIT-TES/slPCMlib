within slPCMlib.Media_Rubitherm_RT;
package RT55 "Rubitherm RT55; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT55";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.211500000000000e+02, 3.311500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.211500000000000e+02, 3.311500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.092116632658725e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+25
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
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 5, 5, 5, 5, 5, 6, 5, 5, 6, 1};
    constant Real[13] breaks = {-5.200000000000000e+01, 4.800000000000000e+01, 4.900000000000000e+01, 5.000000000000000e+01, 5.100000000000000e+01, 5.200000000000000e+01, 5.300000000000000e+01, 5.400000000000000e+01, 5.500000000000000e+01, 5.600000000000000e+01, 5.700000000000000e+01, 5.800000000000000e+01, 1.580000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.796189748103960e-14, 2.055858429391362e-16, 3.173117582142664e-19, -2.891205793294678e-19, 9.319786679628692e-04, -3.139831995385271e-04},  {6.179954684425097e-04, 2.775994142609030e-03, 4.610038686551075e-03, 3.040122688858150e-03, -4.985465326356029e-05, 5.460520954336038e-05, 0.000000000000000e+00},  {1.104890154274057e-02, 2.119004701694835e-02, 1.397733092897777e-02, 3.386756171237512e-03, 2.231713944532416e-04, 2.697404749602157e-04, 0.000000000000000e+00},  {5.009594752931765e-02, 6.154636534123765e-02, 2.817403255901200e-02, 6.976846498652636e-03, 1.571873769254320e-03, -1.946860628545134e-04, 0.000000000000000e+00},  {1.481703796346197e-01, 1.441390347179526e-01, 5.658895404195054e-02, 1.131748094712478e-02, 5.984434549817534e-04, -8.410237517043880e-04, 0.000000000000000e+00},  {3.599732690449250e-01, 2.894580407046541e-01, 8.572182009617178e-02, 5.301017250007919e-03, -3.606675303540186e-03, 1.236521381135697e-03, 0.000000000000000e+00},  {7.380839931733544e-01, 4.685606383385430e-01, 9.235003383631143e-02, 3.239529847204141e-03, 2.575931602138299e-03, 2.095426175901370e-03, -1.697925368026489e-03},  {1.305207627605426e+00, 6.735726006327032e-01, 1.130095942293703e-01, 5.390106542412509e-04, -1.241581803875219e-02, 3.390453447332612e-03, 0.000000000000000e+00},  {2.083303468530322e+00, 8.684978161357751e-01, 7.403625243290649e-02, -1.521972702744139e-02, 4.536449197910868e-03, -2.150710679408691e-03, 0.000000000000000e+00},  {3.013003548590064e+00, 9.783033833138159e-01, 3.408865974396005e-02, -1.858103702988483e-02, -6.217104199132586e-03, 1.054799446827152e-02, -3.101524542815000e-03},  {4.008043920344279e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 5, 5, 5, 5, 6, 5, 5, 6, 5, 1};
    constant Real[13] breaks = {-5.200000000000000e+01, 4.800000000000000e+01, 4.900000000000000e+01, 5.000000000000000e+01, 5.100000000000000e+01, 5.200000000000000e+01, 5.300000000000000e+01, 5.400000000000000e+01, 5.500000000000000e+01, 5.600000000000000e+01, 5.700000000000000e+01, 5.800000000000000e+01, 1.580000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -9.925256005767123e-14, -1.107443632508757e-15, 2.140035311269049e-19, 0.000000000000000e+00, 9.099172238547851e-04, -3.257351861266657e-04},  {5.841820376277596e-04, 2.595175002417798e-03, 4.213144446646817e-03, 2.584468516014537e-03, -3.364416726260602e-04, 2.814838576269982e-04, 0.000000000000000e+00},  {9.922012187707848e-03, 1.883652204138007e-02, 1.276273853520398e-02, 4.053540401780277e-03, 1.070977615508931e-03, -1.389191016742324e-04, 0.000000000000000e+00},  {4.650687167990688e-02, 6.011193527078518e-02, 2.996003441685598e-02, 6.948259847073677e-03, 3.763821071377686e-04, 2.936396453897657e-04, 0.000000000000000e+00},  {1.441971229671493e-01, 1.438505103011960e-01, 5.599950305480105e-02, 1.139018472952241e-02, 1.844580334086597e-03, -1.431532295803374e-03, 0.000000000000000e+00},  {3.558503690909519e-01, 2.902407304566596e-01, 8.692221628985371e-02, 4.453183107835058e-03, -5.313081144930274e-03, 5.017127682753878e-03, -1.416482011183976e-03},  {7.357540634719399e-01, 4.727791341268401e-01, 9.732732540355662e-02, 5.042495131973220e-03, -1.474672898920525e-03, 8.595969660905420e-04, 0.000000000000000e+00},  {1.310287942201480e+00, 6.809605635646436e-01, 1.122027430668586e-01, 7.739773197196539e-03, 2.823311931532185e-03, -5.430996888285126e-03, 0.000000000000000e+00},  {2.108583337073425e+00, 9.127236325746605e-01, 9.805196536479009e-02, -3.527694795952598e-02, -2.433167250989344e-02, 2.658094390983769e-02, -7.007037903557344e-03},  {3.079324220549737e+00, 9.965325215140648e-01, 6.934956971869696e-03, -6.934956971869698e-03, 3.467478485934850e-03, -6.934956971869700e-04, 0.000000000000000e+00},  {4.078630724852549e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT55</strong>.<br><br>
  Information taken from: data_sheet - last access 01.12.2019.<br><br>
  It also contains the phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  false</li>
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
  <li>file creation date: 07-Jul-2022  </ul>
  </html>"));
end RT55;
