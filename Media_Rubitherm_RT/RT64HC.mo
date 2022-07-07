within slPCMlib.Media_Rubitherm_RT;
package RT64HC "Rubitherm RT64HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT64HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.331500000000000e+02, 3.401500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.311500000000000e+02, 3.381500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.341401856521931e+05
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
    constant Integer pieces =  9;
    constant Integer[9] order =  {1, 6, 5, 5, 5, 6, 5, 6, 1};
    constant Real[10] breaks = {-4.000000000000000e+01, 6.000000000000000e+01, 6.100000000000000e+01, 6.200000000000000e+01, 6.300000000000000e+01, 6.400000000000000e+01, 6.500000000000000e+01, 6.600000000000000e+01, 6.700000000000000e+01, 1.670000000000000e+02};
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.035204161868890e-13, 7.921497050362154e-16, 6.087220848470331e-19, 0.000000000000000e+00, 4.006539094938877e-03, -1.705062736990743e-03},  {2.301476358052447e-03, 9.802319052843326e-03, 1.448944989452832e-02, 5.964136209573907e-03, -5.543245580166765e-03, 1.561580612272777e-03, 0.000000000000000e+00},  {2.857571654710401e-02, 4.230854821130248e-02, 1.473819116497710e-02, -5.930399883653883e-04, 2.264657481197118e-03, 5.916942344036100e-03, 0.000000000000000e+00},  {9.321101576025141e-02, 1.086491522212059e-01, 8.571643952742522e-02, 6.763501337678407e-02, 3.184936920137762e-02, -1.880113454054765e-02, 0.000000000000000e+00},  {3.682598555464965e-01, 5.163788755094536e-01, 2.917063494605688e-01, 7.021144776818101e-03, -6.215630350136061e-02, 7.192041258508141e-03, 4.641197674224536e-03},  {1.133043160724709e+00, 9.360371870933066e-01, 8.137034048130751e-02, -7.685970315905219e-02, 4.342186790454814e-02, -9.999531967489947e-03, 0.000000000000000e+00},  {2.107013321077329e+00, 9.918885703595003e-01, 1.132711875654024e-02, -3.167551215759082e-03, -6.575791932901587e-03, 6.210898911048996e-03, -1.631913508156226e-03},  {3.105064652447601e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  9;
    constant Integer[9] order =  {1, 6, 5, 5, 5, 6, 5, 6, 1};
    constant Real[10] breaks = {-4.200000000000000e+01, 5.800000000000000e+01, 5.900000000000000e+01, 6.000000000000000e+01, 6.100000000000000e+01, 6.200000000000000e+01, 6.300000000000000e+01, 6.400000000000000e+01, 6.500000000000000e+01, 1.650000000000000e+02};
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.296145749127707e-14, -6.047643963137829e-18, 6.641810590317428e-18, -1.156482317317871e-18, 3.520711840444073e-03, -1.398802828661500e-03},  {2.121909011759611e-03, 9.210742230204147e-03, 1.422507597451823e-02, 7.231061831210723e-03, -3.378483227702141e-03, 3.546556902399950e-03, 0.000000000000000e+00},  {3.295686272239052e-02, 6.357293127407336e-02, 5.111293112593705e-02, 2.918269794440162e-02, 1.435430128429761e-02, -9.987487416381897e-03, 0.000000000000000e+00},  {1.811922369347183e-01, 2.608266554144715e-01, 1.249119585011086e-01, -1.327497108222692e-02, -3.558313579761187e-02, 1.480080783863671e-02, 0.000000000000000e+00},  {5.328735518090962e-01, 4.024971551724332e-01, 1.959630885512366e-02, -7.599435886307450e-03, 3.842090339557166e-02, -4.645604340800498e-02, 1.870339413049146e-02},  {9.580358340684039e-01, 4.525152265491033e-01, 4.331389944695341e-02, 5.559162622575860e-02, 8.669159831291867e-02, -4.734889985102701e-02, 0.000000000000000e+00},  {1.548799284752111e+00, 8.159397981171126e-01, 2.567493694914713e-01, -7.113097903283677e-02, -1.500529009422163e-01, 1.413816144636242e-01, -3.712367809172696e-02},  {2.504562508757539e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.800000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.800000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT64HC</strong>.<br><br>
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
end RT64HC;
