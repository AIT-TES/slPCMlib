within slPCMlib.Media;
package RT69HC "Rubitherm RT69HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT69HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={
        3.381500000000000e+02,3.441500000000000e+02}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
        3.371500000000000e+02,3.431500000000000e+02}
      "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={
        2.000000000000000e+03,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={
        2.000000000000000e+03,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=1.943941166457007e+05
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
    constant Integer[8] order =  {1, 6, 6, 5, 5, 6, 5, 1};
    constant Real[9] breaks = {-3.500000000000000e+01, 6.500000000000000e+01, 6.600000000000000e+01, 6.700000000000000e+01, 6.800000000000000e+01, 6.900000000000000e+01, 7.000000000000000e+01, 7.100000000000000e+01, 1.710000000000000e+02};
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.346970213900130e-14, -2.887112011460491e-17, 2.684812800961073e-18, -5.782411586589357e-19, 1.679393190189528e-03, -5.758916459755023e-04},  {1.103501544200529e-03, 4.941616075080863e-03, 8.155557212262715e-03, 5.276098982385229e-03, -2.414087386848978e-04, 2.490894382270961e-04, -9.132398648981648e-05},  {1.939313052698172e-02, 3.681289576421974e-02, 2.365643631223273e-02, 4.974878680120268e-03, -3.658213448966645e-04, 4.403382220367924e-03, 0.000000000000000e+00},  {8.887490215902571e-02, 1.196040301512996e-01, 8.041996648689279e-02, 4.754541550421285e-02, 2.165108975694295e-02, -8.403385723055384e-03, 0.000000000000000e+00},  {3.496920183353185e-01, 4.676676400501903e-01, 2.689288943106352e-01, 5.011591730143081e-02, -2.036583885833397e-02, -5.269825933783667e-02, 2.252088605210612e-02},  {1.085861257853511e+00, 9.460438447658840e-01, 1.079123104681491e-01, -1.079123104681491e-01, 5.395615523407458e-02, -1.079123104681492e-02, 0.000000000000000e+00},  {2.075070026806654e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  8;
    constant Integer[8] order =  {1, 6, 5, 5, 5, 6, 6, 1};
    constant Real[9] breaks = {-3.600000000000000e+01, 6.400000000000000e+01, 6.500000000000000e+01, 6.600000000000000e+01, 6.700000000000000e+01, 6.800000000000000e+01, 6.900000000000000e+01, 7.000000000000000e+01, 1.700000000000000e+02};
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 3.601632883088569e-12, 3.462042883051432e-16, 2.559766931375142e-18, -1.156482317317871e-18, 4.781546078636267e-03, -2.072582476163658e-03},  {2.708963606074590e-03, 1.147223553973406e-02, 1.672672364390814e-02, 6.363811263089513e-03, -7.181006749273537e-03, 3.264964889567938e-03, 0.000000000000000e+00},  {3.335569219310070e-02, 5.161791406775563e-02, 2.538176583321484e-02, 1.028943316167475e-02, 9.143817698566150e-03, -3.381436900253700e-03, 0.000000000000000e+00},  {1.264071860540584e-01, 1.529178315121241e-01, 7.729860250709900e-02, 1.305033495340235e-02, -7.763366802702347e-03, 1.898065807258566e-03, 0.000000000000000e+00},  {3.638086540312401e-01, 3.251029032114023e-01, 8.885006462367762e-02, 9.775258151791049e-04, 1.726962233590479e-03, 1.123887930731337e-01, -5.812272053371119e-02},  {8.347321824545121e-01, 7.258511010004539e-01, 3.541915381964270e-01, -3.068110519334565e-02, -3.081698804064088e-01, 2.557402358831308e-01, -6.470208660061634e-02},  {1.766961985334153e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 9.400000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 8.400000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT69HC</strong>.
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
end RT69HC;
