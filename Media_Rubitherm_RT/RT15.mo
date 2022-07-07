within slPCMlib.Media_Rubitherm_RT;
package RT15 "Rubitherm RT15; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT15";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.801500000000000e+02, 2.931500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.801500000000000e+02, 2.911500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 8.275752595434210e+04
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
    constant Integer pieces =  15;
    constant Integer[15] order =  {1, 6, 5, 5, 5, 6, 5, 5, 5, 5, 5, 5, 5, 6, 1};
    constant Real[16] breaks = {-9.300000000000000e+01, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 1.200000000000000e+02};
    constant Real[15,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 6.962062647661589e-15, 1.682520729656082e-16, 3.672465419918411e-18, 1.156482317317871e-18, 4.586832190645534e-03, -1.830281797429759e-03},  {2.756550393222910e-03, 1.195247016865446e-02, 1.841409494500908e-02, 9.262685957860152e-03, -4.520066008218724e-03, 1.269385190655635e-03, 0.000000000000000e+00},  {3.913512064718352e-02, 6.483537985265359e-02, 3.177560867583348e-02, 3.876273831541603e-03, 1.826859945059448e-03, -8.617139365364811e-04, 0.000000000000000e+00},  {1.405875290157352e-01, 1.430142887965002e-01, 4.574845047545015e-02, 2.566574246414585e-03, -2.481709737622957e-03, 9.462254442630523e-04, 0.000000000000000e+00},  {3.303813582407402e-01, 2.370152007574682e-01, 4.802016923158670e-02, 2.101989738553280e-03, 2.249417483692304e-03, -1.738365365060186e-03, 4.108865725268735e-04},  {6.184406566595073e-01, 3.421326709809315e-01, 5.660228828670162e-02, 1.933737473258103e-03, -2.791107537055241e-04, -1.293143371717608e-04, 0.000000000000000e+00},  {1.018700928309521e+00, 4.593754452734287e-01, 5.943569281252518e-02, -4.758489132816012e-04, -9.256824395643279e-04, 3.733791870175450e-04, 0.000000000000000e+00},  {1.536483914229647e+00, 5.749834503354662e-01, 5.618784330546990e-02, -4.447868013634630e-04, 9.412134955233967e-04, -1.557217810448225e-04, 0.000000000000000e+00},  {2.167995912783698e+00, 6.890110216191855e-01, 5.894354606407168e-02, 1.762849370281898e-03, 1.626045902992841e-04, -3.311689765500743e-04, 0.000000000000000e+00},  {2.917544765450986e+00, 8.111812353366213e-01, 6.189603195121234e-02, -8.984220340217087e-04, -1.493240292451087e-03, -3.733448819895754e-04, 0.000000000000000e+00},  {3.787857025530357e+00, 9.244383475572270e-01, 4.650787527454489e-02, -1.060483202372181e-02, -3.359964702398964e-03, 1.708449412763812e-03, 0.000000000000000e+00},  {4.746546901048772e+00, 9.807419902893745e-01, 1.161808511662379e-02, -6.960196705679549e-03, 5.182282361420096e-03, -1.588895687984235e-03, 0.000000000000000e+00},  {5.735540166422526e+00, 9.958822214113421e-01, 5.942232288263351e-03, -2.120024139841509e-03, -2.762196078501080e-03, 2.845764104753316e-03, -7.644416296843667e-04},  {6.734563722378859e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  13;
    constant Integer[13] order =  {1, 6, 5, 5, 5, 5, 5, 6, 5, 5, 5, 6, 1};
    constant Real[14] breaks = {-9.300000000000000e+01, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.180000000000000e+02};
    constant Real[13,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -4.204920851237877e-15, -9.864948490132651e-17, -6.548255684427638e-19, 2.891205793294678e-19, 8.223994974494736e-04, -1.522557392575475e-04},  {6.701437581876221e-04, 3.198463051696641e-03, 5.940158885631400e-03, 5.178880189343787e-03, 1.828161398384156e-03, -8.191560411005541e-04, 0.000000000000000e+00},  {1.599665124214305e-02, 3.383228677902971e-02, 2.425420743296227e-02, 4.299965371874869e-03, -2.267618807118615e-03, 7.891376923460007e-04, 0.000000000000000e+00},  {7.690462971123727e-02, 9.011581099383389e-02, 3.143976762933518e-02, 3.120867066860417e-03, 1.678069654611388e-03, -5.917885670460939e-04, 0.000000000000000e+00},  {2.026673564888321e-01, 1.661112832363029e-01, 4.495290108712387e-02, 3.915260014845031e-03, -1.280873180619081e-03, 4.117475262576703e-04, 0.000000000000000e+00},  {4.167776751727424e-01, 2.646981103638965e-01, 5.313091731052116e-02, 2.909242554945408e-03, 7.778644506692699e-04, -4.815573550188538e-04, 0.000000000000000e+00},  {7.378122524977558e-01, 3.803913436773577e-01, 6.171025812918445e-02, 1.205126807433949e-03, -1.629922324424999e-03, 1.774821079995047e-03, -5.752985608555976e-04},  {1.180688581306446e+00, 5.063298650951694e-01, 6.466483699205280e-02, 9.276770925724729e-04, -1.385295337283727e-03, 1.453377626826230e-03, 0.000000000000000e+00},  {1.752679042775784e+00, 6.401682771419878e-01, 7.366987251433015e-02, 9.920272011699868e-03, 5.881592796847424e-03, -3.732086481705673e-03, 0.000000000000000e+00},  {2.478586970758944e+00, 8.221347769846065e-01, 1.013993805134575e-01, -3.874221617967160e-03, -1.277883961168094e-02, 2.738651057308228e-03, 0.000000000000000e+00},  {3.388206718084668e+00, 9.758887699974459e-01, 4.049018856255288e-02, -2.760306949160863e-02, 9.144156748601995e-04, 7.549388307594428e-03, -2.577423814188822e-03},  {4.382868987321324e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT15</strong>.<br><br>
  Information taken from: data_sheet - last access 01.12.2019.<br><br>
  It also contains the phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
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
end RT15;
