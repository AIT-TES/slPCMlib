within slPCMlib.Media;
package RT90HC "Rubitherm RT90HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT90HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {3.601500000000000e+02, 3.681500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {3.601500000000000e+02, 3.681500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.141249992234551e+05
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
    constant Integer pieces =  10;
    constant Integer[10] order =  {1, 6, 6, 5, 5, 5, 5, 6, 5, 1};
    constant Real[11] breaks = {-1.300000000000000e+01, 8.700000000000000e+01, 8.800000000000000e+01, 8.900000000000000e+01, 9.000000000000000e+01, 9.100000000000000e+01, 9.200000000000000e+01, 9.300000000000000e+01, 9.400000000000000e+01, 9.500000000000000e+01, 1.950000000000000e+02};
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 2.242703954606984e-13, 2.505540032053055e-15, 5.981793775126519e-19, -2.312964634635743e-18, 6.742154374676128e-03, -2.964105342974736e-03},  {3.778049031928167e-03, 1.592613981571192e-02, 2.295996360214218e-02, 8.139436887266555e-03, -1.075080827124040e-02, 8.503065901228318e-03, -2.370689619516000e-03},  {4.618515734752073e-02, 7.155233638582918e-02, 3.234373935604205e-02, 2.753070424268130e-03, -3.795823057838816e-03, 5.930393848892735e-03, 0.000000000000000e+00},  {1.549688743047140e-01, 1.589677033837813e-01, 7.713195077074042e-02, 4.687371668184021e-02, 2.585614618662485e-02, -1.387425041699006e-02, 0.000000000000000e+00},  {4.499241409107108e-01, 4.879060876322340e-01, 2.341474737661085e-01, 1.155579725843907e-02, -4.351510589832542e-02, 1.080286597917584e-02, 0.000000000000000e+00},  {1.150821259648343e+00, 8.708223332422801e-01, 1.157528899432309e-01, -5.447596654310419e-02, 1.049922399755379e-02, 9.006337728359709e-05, 0.000000000000000e+00},  {2.093509803665587e+00, 9.813474263760777e-01, 1.622096807207719e-02, -1.157843678005308e-02, 1.094954088397177e-02, -6.919910637489789e-03, 1.685588084520365e-03},  {3.085214979664691e+00, 9.983661910356711e-01, 3.267617928656587e-03, -3.267617928656594e-03, 1.633808964328297e-03, -3.267617928656594e-04, 0.000000000000000e+00},  {4.084888217871826e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  10;
    constant Integer[10] order =  {1, 6, 5, 6, 5, 5, 5, 5, 6, 1};
    constant Real[11] breaks = {-1.300000000000000e+01, 8.700000000000000e+01, 8.800000000000000e+01, 8.900000000000000e+01, 9.000000000000000e+01, 9.100000000000000e+01, 9.200000000000000e+01, 9.300000000000000e+01, 9.400000000000000e+01, 9.500000000000000e+01, 1.950000000000000e+02};
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.602981865336943e-15, 3.907294273922313e-18, -1.302431424640771e-18, -2.312964634635743e-18, 8.736739166473944e-03, -3.939421989676175e-03},  {4.797317176795168e-03, 2.004716389431769e-02, 2.827606181959680e-02, 8.578951871215929e-03, -1.540763401277291e-02, 7.948471556129537e-03, 0.000000000000000e+00},  {5.424033230528220e-02, 8.044796487671749e-02, 4.105182891790250e-02, 2.643313138141958e-02, 2.433472376787478e-02, 3.667633733978414e-02, -2.148974121466599e-02},  {2.416945773743147e-01, 3.936331513391734e-01, 3.107768208472615e-01, 6.074057555744031e-02, -1.146297077531944e-01, 3.075221347185403e-02, 0.000000000000000e+00},  {9.229676308368497e-01, 8.926507560525782e-01, 1.127424357189564e-01, -9.025612073679691e-02, 3.913135960607576e-02, -7.135350654757565e-03, 0.000000000000000e+00},  {1.870100710822906e+00, 9.682159504305817e-01, 5.408724597444568e-03, -5.084188860069547e-03, 3.454606332287929e-03, -5.485039793804496e-04, 0.000000000000000e+00},  {2.841547299343770e+00, 9.748567384775120e-01, 5.398756217159005e-03, 3.249196675277673e-03, 7.120864353856813e-04, -9.759088305307463e-04, 0.000000000000000e+00},  {3.824788168318573e+00, 9.933706425265471e-01, 9.659776549998654e-03, -3.661545888487064e-03, -4.167457717268050e-03, 4.432429940360558e-03, -1.199646132302316e-03},  {4.823222367597421e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 9.500000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT90HC</strong>.
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
end RT90HC;
