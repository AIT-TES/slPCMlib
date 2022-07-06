within slPCMlib.Media_Rubitherm_RT;
package RT60 "Rubitherm RT60, data taken from data_sheet"
   extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT60";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={
        3.261500000000000e+02,3.351500000000000e+02}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
        3.251500000000000e+02,3.341500000000000e+02}
      "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={
        2.000000000000000e+03,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={
        2.000000000000000e+03,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=9.477650239292810e+04
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
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 5, 5, 5, 6, 5, 5, 6, 1};
    constant Real[12] breaks = {-4.700000000000000e+01, 5.300000000000000e+01, 5.400000000000000e+01, 5.500000000000000e+01, 5.600000000000000e+01, 5.700000000000000e+01, 5.800000000000000e+01, 5.900000000000000e+01, 6.000000000000000e+01, 6.100000000000000e+01, 6.200000000000000e+01, 1.620000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 4.441246673195701e-15, 1.210555389019960e-16, 2.294038792954656e-18, -1.156482317317871e-18, 3.138644127560932e-03, -1.197709531675586e-03},  {1.940934595889910e-03, 8.506963447757008e-03, 1.342079830047568e-02, 7.432250642097599e-03, -2.272422337329130e-03, 8.370173776138175e-04, 0.000000000000000e+00},  {2.986554202650489e-02, 5.274070951375084e-02, 3.045318997893180e-02, 6.712735068919251e-03, 1.912664550739957e-03, -8.206178119997574e-04, 0.000000000000000e+00},  {1.208642233268470e-01, 1.373328638213310e-01, 5.386120437013165e-02, 6.157215151881502e-03, -2.190424509258830e-03, 8.539284543925340e-04, 0.000000000000000e+00},  {3.168790106153248e-01, 2.590348622521647e-01, 6.772958731414847e-02, 5.934801658771520e-03, 2.079217762703840e-03, -4.839169761689836e-04, 0.000000000000000e+00},  {6.511735626269443e-01, 4.181957280267461e-01, 9.317012910499621e-02, 9.412502947897038e-03, -3.403671181410779e-04, 2.003867552963325e-03, -1.283793462077959e-03},  {1.172331629679328e+00, 6.337286036002103e-01, 1.201472088383046e-01, 2.413840763406791e-03, -9.577931284493841e-03, 1.926244704262009e-03, 0.000000000000000e+00},  {1.920969596301018e+00, 8.525840419503793e-01, 8.918359046418216e-02, -1.663543733194848e-02, 5.329223681620488e-05, -5.005640516341569e-04, 0.000000000000000e+00},  {2.845654519568813e+00, 9.787552595719895e-01, 3.459139137289230e-02, -2.142790890102524e-02, -2.449528021354579e-03, 8.387995087391239e-03, -2.632696494373441e-03},  {3.840879032184333e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 5, 5, 6, 5, 5, 5, 6, 1};
    constant Real[12] breaks = {-4.800000000000000e+01, 5.200000000000000e+01, 5.300000000000000e+01, 5.400000000000000e+01, 5.500000000000000e+01, 5.600000000000000e+01, 5.700000000000000e+01, 5.800000000000000e+01, 5.900000000000000e+01, 6.000000000000000e+01, 6.100000000000000e+01, 1.610000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -8.725769018050897e-14, -1.592732853574890e-15, 3.976223245371303e-18, 0.000000000000000e+00, 5.922174179221695e-03, -2.669624661771406e-03},  {3.252549517361443e-03, 1.359312292538295e-02, 1.917737186564418e-02, 5.829248556788857e-03, -1.043349903046261e-02, 3.623686646050228e-03, 0.000000000000000e+00},  {3.504248048076505e-02, 4.582004943545186e-02, 1.030098981373760e-02, 3.321188954407083e-04, 7.684934199788530e-03, -2.269619611413847e-03, 0.000000000000000e+00},  {9.691095321376990e-02, 8.681002449132157e-02, 3.471075558465222e-02, 8.375659580456364e-03, -3.663163857280703e-03, 1.262250065064864e-03, 0.000000000000000e+00},  {2.244064790779842e-01, 1.730171092982015e-01, 5.048125183298582e-02, 6.345504801982195e-03, 2.648086468043614e-03, -4.031694661653104e-03, 1.450900660318761e-03},  {4.543176374778630e-01, 2.921554038959486e-01, 6.685284833544460e-02, 5.638917264000844e-03, 4.253123064559518e-03, -1.585878423585194e-03, 0.000000000000000e+00},  {8.216320516142315e-01, 4.518609524991477e-01, 9.342955427895221e-02, 6.792625286386971e-03, -3.676269053366454e-03, 2.434550920170371e-03, 0.000000000000000e+00},  {1.372473465545522e+00, 6.565656153036032e-01, 1.160953250196182e-01, 1.643305827462486e-02, 8.496485547485397e-03, -1.001661331242470e-02, 0.000000000000000e+00},  {2.160047336378429e+00, 9.219583157945533e-01, 1.162072800041585e-01, -4.974713265968055e-02, -4.158658101463810e-02, 4.819340460961465e-02, -1.329202946889567e-02},  {3.141780593643541e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT60</strong>.
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
end RT60;
