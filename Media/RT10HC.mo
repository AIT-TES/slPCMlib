within slPCMlib.Media;
package RT10HC "Rubitherm RT10HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT10HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={
        2.781500000000000e+02,2.841500000000000e+02}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
        2.781500000000000e+02,2.841500000000000e+02}
      "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={
        2.000000000000000e+03,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={
        2.000000000000000e+03,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=1.344092867525681e+05
      "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature Tref=273.15 + 25
      "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy href=0.0
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
    (xi, dxi) :=BasicUtilities.quartQuintSplineEval(
        T - 273.15,
        pieces,
        order,
        breaks,
        coefs[:, :]);
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
    (xi, dxi) :=BasicUtilities.quartQuintSplineEval(
        T - 273.15,
        pieces,
        order,
        breaks,
        coefs[:, :]);
  end phaseFrac_complSolidification;
  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  8;
    constant Integer[8] order =  {1, 6, 5, 5, 6, 6, 5, 1};
    constant Real[9] breaks = {-9.500000000000000e+01, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.110000000000000e+02};
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.091944350084455e-13, -1.945498336301235e-15, -9.810239719716309e-19, 0.000000000000000e+00, 1.983481923188432e-03, -6.987751075583268e-04},  {1.284706815418964e-03, 5.724758970384065e-03, 9.353192618507514e-03, 5.859317080717780e-03, -5.642169974327426e-04, 1.008417916998913e-04, 0.000000000000000e+00},  {2.175860027929547e-02, 4.025643641830360e-02, 2.455425979306315e-02, 4.610867007985726e-03, -6.000803893328627e-05, 2.068658596953311e-03, 0.000000000000000e+00},  {9.318881405666797e-02, 1.133008178574180e-01, 5.871339855295370e-02, 2.505742082178538e-02, 1.028328494583327e-02, 9.043161651328768e-02, -4.733950589970579e-02},  {3.436358468482401e-01, 5.151540643804321e-01, 3.898089473306015e-01, 2.371660774387954e-02, -2.476512209833151e-01, 1.613838881363044e-01, -3.530974089140154e-02},  {1.150738392564741e+00, 9.703778936726664e-01, 5.924421265436620e-02, -5.924421265436761e-02, 2.962210632718380e-02, -5.924421265436761e-03, 0.000000000000000e+00},  {2.144813971299153e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  8;
    constant Integer[8] order =  {1, 6, 5, 6, 5, 5, 6, 1};
    constant Real[9] breaks = {-9.500000000000000e+01, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.110000000000000e+02};
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 6.701996963605337e-14, 4.668392704075928e-16, 2.537725222093146e-18, 5.782411586589357e-19, 2.390817899663141e-03, -9.891114709613133e-04},  {1.401706428769317e-03, 6.019420672619853e-03, 9.071506932212210e-03, 4.125949577405141e-03, -2.882582566103999e-03, 8.692032501435619e-04, 0.000000000000000e+00},  {1.860520429504608e-02, 2.935596925555720e-02, 1.284589276923923e-02, 1.287651814424766e-03, 1.463433684613811e-03, 4.997302701221627e-04, 1.494385756607880e-03},  {6.555226784561113e-02, 7.622941086599712e-02, 5.290253937053604e-02, 4.202640438625924e-02, 2.637787138434282e-02, -9.813465620927731e-03, 0.000000000000000e+00},  {2.532750282318186e-01, 3.645578601987230e-01, 2.391143246260944e-01, 4.940323371435321e-02, -2.268945672029583e-02, -5.548864632561348e-03, 0.000000000000000e+00},  {8.781121254181321e-01, 8.724940605499795e-01, 1.956986391217655e-01, -9.684323949244358e-02, -5.043377988310257e-02, 6.939999575421514e-02, -1.977107992586487e-02},  {1.848656721542681e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT10HC</strong>.
  It also contains the phase transition functions for complete melting
  (and solidification), which are modelled by piece-wise splines,
  see
  <blockquote>
  <p>
  Barz, T., Krämer, J., & Emhofer, J. (2020). Identification of Phase
  Fraction–Temperature Curves from Heat Capacity Data for Numerical
  Modeling of Heat Transfer in Commercial Paraffin Waxes.
  Energies, 13(19), 5149.
  <a href>doi.org/10.3390/en13195149</a>.
  </p>
  </blockquote>
  <p>
  </p></html>",
  revisions="<html>
  <ul>
  <li>01-Jun-2022  </ul>
  </html>"));
end RT10HC;
