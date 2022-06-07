within slPCMlib.Media;
package RT2HC "Rubitherm RT2HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT2HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {2.681500000000000e+02, 2.781500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {2.681500000000000e+02, 2.771500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.318094403826373e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.SIunits.Temp_K            Tref = 273.15+25
             "reference temperature";
    constant Modelica.SIunits.SpecificEnthalpy  href = 0.0
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
    (Xi, dXi) := BasicUtilities.splineEval(T-273.15,
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
    (Xi, dXi) := BasicUtilities.splineEval(T-273.15,
                     pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complSolidification;

  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 6, 5, 5, 5, 5, 5, 5, 5, 6, 1};
    constant Real[13] breaks = {-1.050000000000000e+02, -5.000000000000000e+00, -4.000000000000000e+00, -3.000000000000000e+00, -2.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 1.050000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 7.539702675046558e-15, 1.668378183311257e-16, -1.323096840003974e-18, 2.891205793294678e-19, 1.011544890151086e-03, -4.177090138396426e-04},  {5.938358763191488e-04, 2.551470367724819e-03, 3.849813693916371e-03, 1.761268624718007e-03, -1.207910756839210e-03, 1.054343789497829e-03, -3.019654181395820e-04},  {8.300856177197382e-03, 1.416318704100582e-02, 7.900111649919680e-03, 1.433755129547818e-03, -4.656730814437948e-04, 2.719806246379032e-04, 0.000000000000000e+00},  {3.160421754086481e-02, 3.376188652690201e-02, 1.212714479627938e-02, 2.290869050151671e-03, 8.942300417457212e-04, -9.100655151540748e-06, 0.000000000000000e+00},  {8.066924730079206e-02, 6.842020016114009e-02, 2.427412564569329e-02, 5.776782665619148e-03, 8.487267659880174e-04, 1.891223334906624e-04, 0.000000000000000e+00},  {1.801782048727233e-01, 1.386393181807935e-01, 4.858805757338555e-02, 1.106291306447784e-02, 1.794338433441329e-03, -3.011012672167036e-04, 0.000000000000000e+00},  {3.799617308576048e-01, 2.746760199186721e-01, 8.953181469529987e-02, 1.522925412607612e-02, 2.888320973578110e-04, 3.133010801314226e-04, 0.000000000000000e+00},  {7.600009527751421e-01, 5.021492454775854e-01, 1.400855804589893e-01, 1.951759331682155e-02, 1.855337498014924e-03, -5.504766988847696e-03, 0.000000000000000e+00},  {1.418103942537706e+00, 8.207707013938474e-01, 1.547227155090665e-01, -2.810872657959568e-02, -2.566849744622355e-02, 1.018190069666733e-02, 0.000000000000000e+00},  {2.350002036111468e+00, 9.941254663716417e-01, 1.820455805961150e-02, -2.896370939781660e-02, 2.524100603711308e-02, -1.150369201034548e-02, 2.151830267640956e-03},  {3.349257495439313e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 5, 5, 5, 5, 6, 5, 6, 1};
    constant Real[12] breaks = {-1.050000000000000e+02, -5.000000000000000e+00, -4.000000000000000e+00, -3.000000000000000e+00, -2.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 1.040000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.742829345710174e-13, -2.582221489978185e-15, -5.389042826672207e-19, 0.000000000000000e+00, 6.338255581415761e-04, -2.045619243485536e-04},  {4.292636336161569e-04, 1.941756244438704e-03, 3.269826716184897e-03, 2.247017094444690e-03, 1.006989254795768e-04, -6.653210073751164e-05, 0.000000000000000e+00},  {7.922030513426513e-03, 1.529259615837311e-02, 9.949750545021306e-03, 1.984491788987879e-03, -2.319615782079814e-04, 2.306794675231269e-04, 0.000000000000000e+00},  {3.514758689512396e-02, 4.137112364016179e-02, 1.681825111796830e-02, 3.363440151387224e-03, 9.214357594076530e-04, -1.966799766310377e-04, 0.000000000000000e+00},  {9.742515758741790e-02, 8.780028948473757e-02, 3.047038636226555e-02, 5.082383422707460e-03, -6.196412374753523e-05, -4.057101755696779e-04, 0.000000000000000e+00},  {2.203105425578113e-01, 1.617118051045515e-01, 4.128865013220592e-02, 7.774251720205421e-04, -2.090515001595925e-03, 4.081685340239959e-03, 0.000000000000000e+00},  {4.260795933052333e-01, 2.586677475798436e-01, 7.189468904109163e-02, 3.323221856803643e-02, 1.831791169960387e-02, 1.781213638941822e-02, -1.299243065275667e-02},  {8.130118659304704e-01, 5.865315261951610e-01, 2.647337190456571e-01, 2.477661620550067e-02, -8.750786614465511e-02, 2.537912559700348e-02, 0.000000000000000e+00},  {1.626924986829138e+00, 9.671929763094377e-01, 6.780762676426408e-02, -7.146359240308500e-02, 3.938776184036227e-02, -1.007113175136431e-02, 7.311931277639533e-04},  {2.620509820716517e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT2HC</strong>.
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
end RT2HC;
