within slPCMlib.Media;
package RT100HC "Rubitherm RT100HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT100HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {3.691500000000000e+02, 3.781500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {3.641500000000000e+02, 3.731500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.059071117908937e+05
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
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 6, 5, 5, 5, 5, 5, 5, 6, 1};
    constant Real[12] breaks = {-4.000000000000000e+00, 9.600000000000000e+01, 9.700000000000000e+01, 9.800000000000000e+01, 9.900000000000000e+01, 1.000000000000000e+02, 1.010000000000000e+02, 1.020000000000000e+02, 1.030000000000000e+02, 1.040000000000000e+02, 1.050000000000000e+02, 2.050000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 4.138911834356316e-13, 3.628911420365913e-15, -1.124461037564610e-18, 5.782411586589357e-19, 1.601317816620632e-03, -7.056357021043549e-04},  {8.956821149337965e-04, 3.772774870887650e-03, 5.428642634644532e-03, 1.900464124119273e-03, -2.577946448462166e-03, 8.759983472954695e-03, -3.598564391841648e-03},  {1.458103637723613e-02, 3.222819773239535e-02, 2.928372516815147e-02, 7.217225222984601e-03, -1.275649496131342e-02, 9.332131051115626e-03, 0.000000000000000e+00},  {7.988582059056976e-02, 1.080819991480081e-01, 6.771774158038128e-02, 4.951255588888720e-02, 3.390416029426471e-02, -1.838505962922329e-02, 0.000000000000000e+00},  {3.207172178728878e-01, 4.357464930062054e-01, 2.358297747203968e-01, 1.278600773713094e-03, -5.802113785185176e-02, 1.725000461618212e-02, 0.000000000000000e+00},  {9.528009531375334e-01, 7.654073164412276e-01, 6.403879609224300e-02, -5.830590447187278e-02, 2.822888522905882e-02, -5.549408019565318e-03, 0.000000000000000e+00},  {1.746620638408625e+00, 8.037356960284943e-01, 3.000313855324337e-03, -8.844437512906744e-04, 4.818451312322356e-04, 2.351719599382499e-03, 0.000000000000000e+00},  {2.555305769271768e+00, 8.207689710070936e-01, 2.675524938267055e-02, 2.456013276746325e-02, 1.224044312814473e-02, -9.880053583404528e-03, 0.000000000000000e+00},  {3.429750511973735e+00, 9.475213726702973e-01, 7.507777061988270e-02, -2.527863055400316e-02, -3.715982478887791e-02, 3.731144899730327e-02, -9.959828013175895e-03},  {4.417262820905163e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 6, 5, 5, 5, 5, 5, 5, 6, 1};
    constant Real[12] breaks = {-9.000000000000000e+00, 9.100000000000000e+01, 9.200000000000000e+01, 9.300000000000000e+01, 9.400000000000000e+01, 9.500000000000000e+01, 9.600000000000000e+01, 9.700000000000000e+01, 9.800000000000000e+01, 9.900000000000000e+01, 1.000000000000000e+02, 2.000000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.454528086174847e-12, -9.470112008004124e-15, 4.496520647995262e-18, -1.156482317317871e-18, 3.623531399774491e-03, -1.641917206905754e-03},  {1.981614190404743e-03, 8.266153754968697e-03, 1.160655589414915e-02, 3.396969859629833e-03, -6.511101104713855e-03, 4.603351768825868e-03, -1.137797407270311e-03},  {2.220574695599413e-02, 3.181574510383254e-02, 1.169741542395963e-02, 6.301349836268741e-04, -5.613033696391774e-04, 1.150965208263616e-03, 0.000000000000000e+00},  {6.693870430603761e-02, 6.061059346534251e-02, 2.172965223964115e-02, 9.894573587706320e-03, 5.193522671678900e-03, -2.746502044139297e-03, 0.000000000000000e+00},  {1.616205442262672e-01, 1.407951991740725e-01, 5.510948859144173e-02, 3.203643833028950e-03, -8.538987549017584e-03, 2.546743916534719e-03, 0.000000000000000e+00},  {3.547366321923275e-01, 2.392028772423480e-01, 3.895393396176912e-02, -5.484867197694196e-03, 4.194732033656011e-03, -2.168196698419503e-03, 0.000000000000000e+00},  {6.294351115339870e-01, 3.065940882153785e-01, 2.598575758642776e-02, -1.038790604726518e-02, -6.646251458441501e-03, 1.241712234246222e-02, 0.000000000000000e+00},  {9.573979221725487e-01, 3.629024911246755e-01, 7.911575411860425e-02, 8.719831154359103e-02, 5.543936025386961e-02, -3.974066131216913e-02, 0.000000000000000e+00},  {1.502313177901120e+00, 8.057830684509971e-01, 2.759402371509179e-01, -8.845086056262033e-02, -1.432639463069761e-01, 1.411464152143670e-01, -3.749787531765725e-02},  {2.455970216530149e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.000000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 8.500000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT100HC</strong>.
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
end RT100HC;
