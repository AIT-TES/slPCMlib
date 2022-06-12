within slPCMlib.Media;
package RT62HC "Rubitherm RT62HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT62HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {3.321500000000000e+02, 3.381500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {3.311500000000000e+02, 3.361500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.942042828941463e+05
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
    (xi, dxi) := BasicUtilities.splineEval(T-273.15,
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
    (xi, dxi) := BasicUtilities.splineEval(T-273.15,
                     pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complSolidification;
  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  8;
    constant Integer[8] order =  {1, 6, 6, 5, 5, 6, 5, 1};
    constant Real[9] breaks = {-4.100000000000000e+01, 5.900000000000000e+01, 6.000000000000000e+01, 6.100000000000000e+01, 6.200000000000000e+01, 6.300000000000000e+01, 6.400000000000000e+01, 6.500000000000000e+01, 1.650000000000000e+02};
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -6.525163553303419e-13, -1.679974163304125e-15, -2.059018448443851e-18, 5.782411586589357e-19, 2.930698552155855e-03, -1.266025438772819e-03},  {1.664673112728838e-03, 7.057340127554294e-03, 1.031660393996474e-02, 3.986476746102154e-03, -4.336888820813021e-03, 2.283851320225506e-03, 6.215934312699372e-04},  {2.159364985703245e-02, 3.745124015110248e-02, 2.841711592469672e-02, 2.190930329050387e-02, 1.640626924936357e-02, -1.952807880340307e-03, 0.000000000000000e+00},  {1.238247705923588e-01, 2.158744194683812e-01, 1.730545624889883e-01, 6.800630148455505e-02, 6.642229847662033e-03, -1.550504222924261e-02, 0.000000000000000e+00},  {5.718972416527026e-01, 7.150461571451719e-01, 2.618764237362015e-01, -6.047520141722286e-02, -7.088298129855099e-02, 5.813550541203828e-02, -1.353874038064474e-02},  {1.462058404849696e+00, 9.832865599480299e-01, 3.342688010393867e-02, -3.342688010393868e-02, 1.671344005196934e-02, -3.342688010393868e-03, 0.000000000000000e+00},  {2.458715716839301e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  7;
    constant Integer[7] order =  {1, 6, 5, 5, 6, 6, 1};
    constant Real[8] breaks = {-4.200000000000000e+01, 5.800000000000000e+01, 5.900000000000000e+01, 6.000000000000000e+01, 6.100000000000000e+01, 6.200000000000000e+01, 6.300000000000000e+01, 1.630000000000000e+02};
    constant Real[7,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.519720267873975e-15, 1.081421753142474e-18, -1.516956235032029e-18, 0.000000000000000e+00, 1.515462100880499e-03, -5.835458500841632e-04},  {9.319162507948160e-04, 4.076035403898648e-03, 6.401433257542538e-03, 3.483704007121723e-03, -1.175877246859954e-03, 7.973596300393411e-04, 0.000000000000000e+00},  {1.451457130253711e-02, 2.661330310309551e-02, 1.777087809814140e-02, 6.753791320075317e-03, 2.810920903336752e-03, 2.673627014324387e-04, 0.000000000000000e+00},  {6.873082742861852e-02, 9.499693038011690e-02, 5.757140449271225e-02, 2.067110194774625e-02, 4.147734410498945e-03, 1.665763190178374e-01, -8.669669135215297e-02},  {3.259976263253773e-01, 6.014454298265496e-01, 5.097839366950243e-01, -3.090859727494356e-02, -4.634210407826087e-01, 3.800094118085700e-01, -9.577506788401607e-02},  {1.227131698713953e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.500000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT62HC</strong>.
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
end RT62HC;
