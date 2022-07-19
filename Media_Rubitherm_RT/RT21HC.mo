within slPCMlib.Media_Rubitherm_RT;
package RT21HC "Rubitherm RT21HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "RT21HC";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.871500000000000e+02, 2.981500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.861500000000000e+02, 2.961500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.656564717778307e+05
             "scalar phase transition enthalpy";
      
    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 2.871500000000000e+02
             "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy  href = 0.0
             "reference enthalpy at Tref";
      
  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  protected
    constant Integer pieces   = data_H.pieces;
    constant Integer order[:] = data_H.order;
    constant Real breaks[:]   = data_H.breaks;
    constant Real coefs[:,:]  = data_H.coefs;
  algorithm 
    (xi, dxi) := BasicUtilities.quartQuintSplineEval(T-273.15, 
                 pieces, order, breaks, coefs[:,:]);     
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
  protected
    constant Integer pieces   = data_C.pieces;
    constant Integer order[:] = data_C.order;
    constant Real breaks[:]   = data_C.breaks;
    constant Real coefs[:,:]  = data_C.coefs;
  algorithm
    (xi, dxi) := BasicUtilities.quartQuintSplineEval(T-273.15, 
                     pieces, order, breaks, coefs[:,:]);     
  end phaseFrac_complSolidification;
      
  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;    
    constant Integer pieces  = 13; 
    constant Integer[13] order  = {1, 6, 5, 5, 5, 5, 6, 5, 5, 5, 6, 5, 1}; 
    constant Real[14] breaks = {-8.600000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 1.250000000000000e+02}; 
    constant Real[13,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 5.650431199991921e-16, -1.279099860848623e-18, 4.263666202828742e-19, 0.000000000000000e+00, 1.288385026784240e-03, -4.570424334754814e-04},  {8.313425933093232e-04, 3.699670533068009e-03, 6.028213765710183e-03, 3.743001598332777e-03, -4.137113682110190e-04, 2.588879602145481e-05, 0.000000000000000e+00},  {1.391440591823073e-02, 2.545970136675114e-02, 1.503383831165694e-02, 2.347044085703248e-03, -2.842673881037450e-04, 6.189303412500542e-05, 0.000000000000000e+00},  {5.653261532836332e-02, 6.174090586538477e-02, 2.098829658139427e-02, 1.828904874538322e-03, 2.519778252128205e-05, 2.906866116087152e-05, 0.000000000000000e+00},  {1.411449890933628e-01, 1.094503480876814e-01, 2.691688451174564e-02, 2.220382616232163e-03, 1.705410883256397e-04, 1.027525399781596e-03, 0.000000000000000e+00},  {2.809306707971292e-01, 1.757650563120786e-01, 4.487653288821193e-02, 1.317780096735068e-02, 5.308168087233621e-03, -2.642833787325452e-03, 6.037470080498611e-04},  {5.180191422727285e-01, 3.166925104511658e-01, 9.888681156115911e-02, 2.005707560402787e-02, 1.150204271354280e-03, -2.538193221459350e-03, 0.000000000000000e+00},  {9.522675509389762e-01, 5.665472113636729e-01, 1.405773317867749e-01, -7.240395251485064e-04, -1.154076183594247e-02, 2.860115722775943e-03, 0.000000000000000e+00},  {1.649987408451109e+00, 8.136672876318856e-01, 9.776179942343402e-02, -1.828592964115896e-02, 2.759816777937242e-03, -1.334043485990518e-03, 0.000000000000000e+00},  {2.544556339157217e+00, 9.587021472370746e-01, 4.612247630767542e-02, -2.058709738931517e-02, -3.910400652015346e-03, 7.415878143125611e-03, -2.085361231222100e-03},  {3.530213981572540e+00, 9.981114284047187e-01, 3.777143190562430e-03, -3.777143190562431e-03, 1.888571595281216e-03, -3.777143190562431e-04, 0.000000000000000e+00},  {4.529836267253484e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer pieces  = 12; 
    constant Integer[12] order  = {1, 6, 5, 5, 5, 5, 5, 6, 5, 5, 6, 1}; 
    constant Real[13] breaks = {-8.700000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 1.230000000000000e+02}; 
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -3.250520126426846e-14, -3.968623305646785e-17, -7.080229699957057e-20, 2.891205793294678e-19, 8.692852370986259e-04, -3.351494981280588e-04},  {5.341357389380224e-04, 2.335529196691830e-03, 3.665609899065337e-03, 1.989862408425082e-03, -6.808162864277535e-04, 1.792077217081689e-04, 0.000000000000000e+00},  {8.023528678400687e-03, 1.380910968292845e-02, 7.342376622855752e-03, 1.058674479795757e-03, 2.152223221130913e-04, -7.822188326941808e-05, 0.000000000000000e+00},  {3.037068990282431e-02, 3.213966624013214e-02, 1.102751516222739e-02, 1.137344935553944e-03, -1.758870942339991e-04, 5.002092406460465e-04, 0.000000000000000e+00},  {7.499953838714983e-02, 5.940422919754536e-02, 1.838631980994569e-02, 5.435888965078410e-03, 2.325159108996233e-03, -8.210921442415472e-04, 0.000000000000000e+00},  {1.597300433244740e-01, 1.176797114274528e-01, 4.043401991674286e-02, 6.525603958647869e-03, -1.780301612211503e-03, 2.208561908825473e-04, 0.000000000000000e+00},  {3.228099332059885e-01, 2.121076376424513e-01, 5.153758402824292e-02, 1.612959418627356e-03, -6.760206577987662e-04, 1.036876055893386e-02, -4.267968833917478e-03},  {5.934928853625278e-01, 3.435535911148166e-01, 9.198841141790878e-02, 1.723710569842128e-02, -1.285175037189166e-02, 6.228978703255797e-03, 0.000000000000000e+00},  {1.039649221925039e+00, 5.589796230746236e-01, 1.288790133143806e-01, 2.811989124341261e-02, 1.829314314438732e-02, -1.644097949289242e-02, 0.000000000000000e+00},  {1.757479913208950e+00, 8.920649985467317e-01, 1.585877509820182e-01, -6.311733110796229e-02, -6.391175432007477e-02, 7.006460278844852e-02, -1.909408397481119e-02},  {2.732074096123300e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT21HC</strong>.<br><br>
  Information taken from: data_sheet - last access 01.12.2019.<br><br>
  It also contains the phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  false</li>
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
  <li>file creation date: 19-Jul-2022  </ul>
  </html>"));
  end RT21HC;
