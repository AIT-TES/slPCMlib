within slPCMlib.Media_Rubitherm_RT;
package RT5 "Rubitherm RT5; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "RT5";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.721500000000000e+02, 2.801500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.711500000000000e+02, 2.801500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.242052446212872e+05
             "scalar phase transition enthalpy";
      
    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 2.721500000000000e+02
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
    constant Integer pieces  = 10; 
    constant Integer[10] order  = {1, 6, 6, 5, 5, 5, 5, 5, 6, 1}; 
    constant Real[11] breaks = {-1.010000000000000e+02, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 1.070000000000000e+02}; 
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.222930562390897e-13, -1.491436324588392e-15, 2.016287414424062e-18, 0.000000000000000e+00, 4.501084331050176e-03, -1.922323497343506e-03},  {2.578760833482887e-03, 1.097148067096661e-02, 1.617599085034769e-02, 6.564373363631655e-03, -6.329430804901707e-03, 3.629438410150234e-03, -6.456428296734195e-04},  {3.294497049400395e-02, 5.197219431565119e-02, 2.450226776823339e-02, 4.628177652058771e-03, 2.133118800748167e-03, -6.096532704255811e-04, 0.000000000000000e+00},  {1.155710757602699e-01, 1.203454716591456e-01, 4.508898082464281e-02, 7.064120150795628e-03, -9.151475513797377e-04, 5.166490796581104e-04, 0.000000000000000e+00},  {2.876711499231322e-01, 2.306384489536056e-01, 6.595694676533247e-02, 8.570020741857781e-03, 1.668097846910814e-03, 9.697744878669871e-04, 0.000000000000000e+00},  {5.954744387187058e-01, 3.997836685368087e-01, 1.113733409510405e-01, 2.494015700817091e-02, 6.516970286245750e-03, -5.557904040736496e-03, 0.000000000000000e+00},  {1.132530671460235e+00, 6.956291824046886e-01, 1.697165932856627e-01, -4.571002254211046e-03, -2.127254991743673e-02, 4.246535526165444e-03, 0.000000000000000e+00},  {1.976279430505104e+00, 9.574918401745548e-01, 7.083364228006422e-02, -4.719584666230350e-02, -3.987228660950563e-05, 1.419065182797865e-02, -4.727559123552251e-03},  {2.966832286715236e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer pieces  = 11; 
    constant Integer[11] order  = {1, 6, 5, 5, 5, 5, 5, 6, 5, 6, 1}; 
    constant Real[12] breaks = {-1.020000000000000e+02, -2.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 1.070000000000000e+02}; 
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 2.251861794113628e-14, 3.872847188610760e-17, 3.281261813747667e-18, 0.000000000000000e+00, 2.944672900006421e-03, -1.159631818113830e-03},  {1.785041081915151e-03, 7.765573591372053e-03, 1.205225172835680e-02, 6.254092637787603e-03, -2.671112771675352e-03, 7.645375526940027e-04, 0.000000000000000e+00},  {2.595038382045026e-02, 4.377059163821739e-02, 2.243322853860752e-02, 3.215017078026224e-03, 1.151574991794662e-03, -2.696408099478623e-04, 0.000000000000000e+00},  {9.625115525714820e-02, 1.015401958669497e-01, 3.629132162397554e-02, 5.124908945726248e-03, -1.966290579446498e-04, 3.156104729089621e-04, 0.000000000000000e+00},  {2.393265631087640e-01, 1.902891020848523e-01, 5.364237884257601e-02, 7.494497443037271e-03, 1.381423306600161e-03, -6.467116798782973e-04, 0.000000000000000e+00},  {4.914872531059514e-01, 3.223494869261196e-01, 7.794729421250582e-02, 6.553073870654942e-03, -1.852135092791326e-03, 1.239153575114545e-03, 0.000000000000000e+00},  {8.977241265975549e-01, 4.966905244675128e-01, 9.888524101886814e-02, 1.153606925063507e-02, 4.343632782781396e-03, -1.786804900508956e-03, -5.330959504572770e-04},  {1.506859693266386e+00, 7.343111451829816e-01, 1.336907572055130e-01, 3.806323675255575e-04, -1.258683097662254e-02, 9.408399945571647e-04, 0.000000000000000e+00},  {2.363596237040341e+00, 9.571914327628736e-01, 6.872006839392611e-02, -4.055829159339294e-02, -7.882631003836715e-03, 1.847359228108726e-02, -5.632355360106639e-03},  {3.353908052520892e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT5</strong>.<br><br>
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
  <li>file creation date: 19-Jul-2022  </ul>
  </html>"));
  end RT5;
