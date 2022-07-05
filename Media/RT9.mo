within slPCMlib.Media;
package RT9 "Rubitherm RT9, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT9";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={
        2.751500000000000e+02,2.861500000000000e+02}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
        2.751500000000000e+02,2.841500000000000e+02}
      "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={
        2.000000000000000e+03,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={
        2.000000000000000e+03,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=1.063272686605944e+05
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
    constant Integer pieces =  13;
    constant Integer[13] order =  {1, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 6, 1};
    constant Real[14] breaks = {-9.800000000000000e+01, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.130000000000000e+02};
    constant Real[13,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.709507411798086e-14, 1.607778920057637e-16, -2.707408706601551e-18, 5.782411586589357e-19, 1.919172226556605e-03, -4.915369813854850e-04},  {1.427635245188374e-03, 6.646639244508039e-03, 1.181866754478412e-02, 9.360982637856343e-03, 2.222806412000747e-03, -2.541788242400887e-03, 5.327327177977872e-04},  {2.946767555973452e-02, 5.774560299042144e-02, 3.581156227331549e-02, 3.488980217806210e-03, -2.495144033036877e-03, 9.931271884315923e-04, 0.000000000000000e+00},  {1.250118041966723e-01, 1.348207280004778e-01, 4.123891061282874e-02, 3.439675969974623e-03, 2.470491909121084e-03, -1.011647745748616e-03, 0.000000000000000e+00},  {3.059699629433260e-01, 2.324413060438109e-01, 5.626441251999306e-02, 3.205166148972800e-03, -2.587746819621995e-03, 1.171996568680703e-03, 0.000000000000000e+00},  {5.964650974051614e-01, 3.500946250956322e-01, 6.207339573598653e-02, 4.574144557291850e-03, 3.272236023781520e-03, -1.318826662599459e-03, 0.000000000000000e+00},  {1.015160672155254e+00, 4.944586610216066e-01, 8.224097892455658e-02, 4.474822026423340e-03, -3.321897289215775e-03, 8.093203456690304e-04, 0.000000000000000e+00},  {1.593822557184294e+00, 6.631240975214786e-01, 8.382726472522231e-02, -7.195636737494550e-04, 7.247044391293771e-04, -1.005801905980919e-03, 0.000000000000000e+00},  {2.339773258290394e+00, 8.264897441772886e-01, 7.595878127894103e-02, -7.878764977041160e-03, -4.304305090775219e-03, 1.294786087359932e-03, 0.000000000000000e+00},  {3.231333499766167e+00, 9.440277218777803e-01, 3.944451667676587e-02, -1.214812446654270e-02, 2.169625346024442e-03, -4.327695358547155e-04, 0.000000000000000e+00},  {4.204394469664340e+00, 9.929870355365324e-01, 1.169019999473749e-02, -7.797318440992079e-03, 5.777666750864223e-06, 2.334573398896932e-03, -7.785763107490350e-04},  {5.202836161509517e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 5, 5, 6, 5, 5, 5, 6, 1};
    constant Real[12] breaks = {-9.800000000000000e+01, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.110000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.801853950947819e-15, 6.060060597417378e-17, 6.661440225800383e-18, 0.000000000000000e+00, 4.004441224539636e-03, -1.420854201917403e-03},  {2.583587022624102e-03, 1.149708091119780e-02, 1.873159921663545e-02, 1.162732820704829e-02, -1.290606906062875e-03, -4.260933442350641e-05, 0.000000000000000e+00},  {4.310637911701927e-02, 7.846678966924481e-02, 4.544384905716804e-02, 6.038807238561730e-03, -1.503653578180407e-03, 4.925378570389481e-05, 0.000000000000000e+00},  {1.716014252895173e-01, 1.817025641150631e-01, 5.503088716080971e-02, 5.167307828790517e-04, -1.257384649660933e-03, 9.742567227473717e-04, 0.000000000000000e+00},  {4.085684794213557e-01, 2.931562758004130e-01, 5.877933883895499e-02, 5.229759411709033e-03, 3.613898964075925e-03, 8.720051945090720e-04, -1.116254662207030e-03},  {7.691035029688107e-01, 4.385223255690582e-01, 8.812824287052295e-02, 6.080313968962848e-03, -8.769894996484168e-03, 2.798060374113150e-03, 0.000000000000000e+00},  {1.295862550754984e+00, 6.119304751016202e-01, 8.173041853963793e-02, -1.018662275842328e-03, 5.220406874081581e-03, -1.632897778852377e-03, 0.000000000000000e+00},  {1.992092291215629e+00, 7.850524639554329e-01, 9.366789516807664e-02, 3.533987431960222e-03, -2.944082020180306e-03, -2.376592411651801e-03, 0.000000000000000e+00},  {2.869025963339267e+00, 9.593309264484854e-01, 6.283944122635741e-02, -3.200826476527902e-02, -1.482704407843931e-02, 2.146411469233515e-02, -6.166235292215762e-03},  {3.859658901570511e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.800000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.600000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT9</strong>.
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
end RT9;
