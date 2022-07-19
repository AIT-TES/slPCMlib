within slPCMlib.Media_Rubitherm_RT;
package RT35HC "Rubitherm RT35HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "RT35HC";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.041500000000000e+02, 3.111500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.041500000000000e+02, 3.101500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.055606478888954e+05
             "scalar phase transition enthalpy";
      
    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 3.041500000000000e+02
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
    constant Integer pieces  = 9; 
    constant Integer[9] order  = {1, 6, 5, 6, 5, 5, 5, 6, 1}; 
    constant Real[10] breaks = {-6.900000000000000e+01, 3.100000000000000e+01, 3.200000000000000e+01, 3.300000000000000e+01, 3.400000000000000e+01, 3.500000000000000e+01, 3.600000000000000e+01, 3.700000000000000e+01, 3.800000000000000e+01, 1.380000000000000e+02}; 
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.240980301983814e-13, -7.631257380724631e-16, 7.316829884788960e-19, 2.891205793294678e-19, 1.556843503193726e-03, -5.504037706263503e-04},  {1.006439732342516e-03, 4.481794892092473e-03, 7.312378472541612e-03, 4.560359619410255e-03, -4.718390434266247e-04, -1.032350238934175e-04, 0.000000000000000e+00},  {1.678589864906681e-02, 3.038409940230950e-02, 1.713007283127872e-02, 1.640653206769566e-03, -9.880141628937124e-04, 1.881827924564478e-03, 1.743879121376852e-03},  {6.857841697247222e-02, 8.548656238454491e-02, 6.110041354052224e-02, 5.138445822837653e-02, 3.457931228058145e-02, -1.683463939615334e-02, 0.000000000000000e+00},  {2.842945240103441e-01, 4.159848162921154e-01, 2.543832679476066e-01, 2.135531338916895e-02, -4.959388470018524e-02, 1.182199381903310e-02, 0.000000000000000e+00},  {9.382460307580829e-01, 8.495517226489878e-01, 1.391058381043321e-01, -5.880028722124099e-02, 9.516084394980258e-03, 5.954524188363268e-05, 0.000000000000000e+00},  {1.877678933927026e+00, 9.897246009831898e-01, 2.039693522932671e-02, -2.014049722248364e-02, 9.813810604398421e-03, -1.808899316773646e-03, -5.128760136867962e-05},  {2.875613596603315e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer pieces  = 8; 
    constant Integer[8] order  = {1, 6, 6, 5, 5, 6, 5, 1}; 
    constant Real[9] breaks = {-6.900000000000000e+01, 3.100000000000000e+01, 3.200000000000000e+01, 3.300000000000000e+01, 3.400000000000000e+01, 3.500000000000000e+01, 3.600000000000000e+01, 3.700000000000000e+01, 1.370000000000000e+02}; 
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -3.143361842836314e-13, -1.359236556995747e-15, -2.028891940762838e-18, 0.000000000000000e+00, 2.660987622903916e-03, -1.253426344186147e-03},  {1.407561278402071e-03, 5.784380048997627e-03, 7.808481066245201e-03, 1.541349345316212e-03, -5.496457048272632e-03, 1.536179138618721e-02, -4.878360622305588e-03},  {2.152874545457010e-02, 5.157835522130870e-02, 5.989629133984571e-02, 3.560622256798601e-02, -1.862909451920407e-03, 2.210629147924184e-03, 0.000000000000000e+00},  {1.689573342797143e-01, 2.817911135378500e-01, 1.776437938115272e-01, 5.026087623954621e-02, 9.190236287700512e-03, -1.468154200847887e-02, 0.000000000000000e+00},  {6.731618121478594e-01, 7.512145649887730e-01, 2.367524201715838e-01, -5.979359869444038e-02, -6.421747375469382e-02, 5.550463388288316e-02, -1.329988472870120e-02},  {1.579322474013264e+00, 9.861925752706520e-01, 2.761484945840738e-02, -2.761484945840802e-02, 1.380742472920401e-02, -2.761484945840803e-03, 0.000000000000000e+00},  {2.576560989067279e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT35HC</strong>.<br><br>
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
  end RT35HC;
