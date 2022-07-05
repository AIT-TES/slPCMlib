within slPCMlib.Media;
package RT54HC "Rubitherm RT54HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT54HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={
        3.241500000000000e+02,3.311500000000000e+02}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
        3.231500000000000e+02,3.281500000000000e+02}
      "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={
        2.000000000000000e+03,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={
        2.000000000000000e+03,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=1.450353656438634e+05
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
    constant Integer pieces =  9;
    constant Integer[9] order =  {1, 5, 5, 6, 6, 5, 5, 6, 1};
    constant Real[10] breaks = {-4.900000000000000e+01, 5.100000000000000e+01, 5.200000000000000e+01, 5.300000000000000e+01, 5.400000000000000e+01, 5.500000000000000e+01, 5.600000000000000e+01, 5.700000000000000e+01, 5.800000000000000e+01, 1.580000000000000e+02};
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -7.859227461957215e-14, -1.128924002298123e-16, 2.144004763415389e-17, 5.782411586589357e-19, 1.839953498960722e-03, 0.000000000000000e+00},  {1.839953498882038e-03, 9.199767494723312e-03, 1.839953498960717e-02, 1.839953498960721e-02, 9.199767494803607e-03, -5.164289750582203e-03, 0.000000000000000e+00},  {5.187426871704114e-02, 1.121750636690788e-01, 7.715384742142845e-02, 3.555707462999695e-03, -1.662168125810741e-02, 1.598421836536363e-01, -7.736968452679406e-02},  {3.106097051392829e-01, 5.456559669758316e-01, 4.259674508962355e-01, -1.190287156894775e-02, -3.779560308918367e-01, 3.052074182852428e-01, -7.653323904685376e-02},  {1.121048399788954e+00, 9.168957876389009e-01, 2.659824798799336e-02, -2.317593220941911e-03, 8.247483157068007e-05, -3.703055468226948e-05, 0.000000000000000e+00},  {2.062270286471795e+00, 9.632842505049292e-01, 1.977001176776901e-02, -2.357999441481885e-03, -1.026779418406673e-04, -2.608067012390614e-04, 0.000000000000000e+00},  {3.042603064659932e+00, 9.940355304424620e-01, 9.471878779888730e-03, -5.376778221235169e-03, -1.406711448035974e-03, 2.738402624799330e-03, -8.190201117307117e-04},  {4.041246366726080e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  7;
    constant Integer[7] order =  {1, 6, 5, 6, 5, 6, 1};
    constant Real[8] breaks = {-5.000000000000000e+01, 5.000000000000000e+01, 5.100000000000000e+01, 5.200000000000000e+01, 5.300000000000000e+01, 5.400000000000000e+01, 5.500000000000000e+01, 1.550000000000000e+02};
    constant Real[7,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -4.844500635737350e-13, -1.529992776148875e-16, -1.041945074341719e-18, 5.782411586589357e-19, 2.811064923098048e-03, -1.263320579324869e-03},  {1.547744343288575e-03, 6.475401139055487e-03, 9.160840541107289e-03, 2.844237644483097e-03, -4.894484074382796e-03, 3.347459689698212e-03, 0.000000000000000e+00},  {1.848119928324986e-02, 3.048915730564884e-02, 2.180124592524193e-02, 1.674089824393396e-02, 1.184281437410827e-02, -1.944548158318221e-02, 1.140767355221289e-02},  {9.131750710121353e-02, 1.429042347817780e-01, 1.197411143530646e-01, 9.781081095280265e-02, 8.573050974139051e-02, -5.581287967595845e-02, 0.000000000000000e+00},  {4.816912972542909e-01, 7.396765369327490e-01, 3.694278089002312e-01, -1.173959468412198e-01, -1.933338886384017e-01, 1.898858949630874e-01, -5.040637241180233e-02},  {1.419545330158935e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.500000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 8.000000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT54HC</strong>.
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
end RT54HC;
