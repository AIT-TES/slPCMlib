within slPCMlib.Media_Rubitherm_SP;
package SP50 "Rubitherm SP50; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP50";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.191500000000000e+02, 3.271500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.161500000000000e+02, 3.221500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.200000000000000e+05
             "scalar phase transition enthalpy";
      
    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+2.500000000000000e+01
             "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy  href = 0.0
             "reference enthalpy at Tref";
      
  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  protected
    constant Integer len_x    = data_H.len_x;
    constant Real data_x[:]   = data_H.data_x;
    constant Real data_y[:]   = data_H.data_y;
    constant Real m_k[:]      = data_H.m_k;
    constant Real iy_start[:] = data_H.iy_start;
    constant Real iy_scaler   = data_H.iy_scaler;
  algorithm 
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15, 
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
  protected
    constant Integer len_x    = data_C.len_x;
    constant Real data_x[:]   = data_C.data_x;
    constant Real data_y[:]   = data_C.data_y;
    constant Real m_k[:]      = data_C.m_k;
    constant Real iy_start[:] = data_C.iy_start;
    constant Real iy_scaler   = data_C.iy_scaler;
  algorithm
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15, 
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complSolidification;
  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 9; 
    constant Real[9] data_x   = {4.600000000000000e+01, 4.700000000000000e+01, 4.800000000000000e+01, 4.900000000000000e+01, 5.000000000000000e+01, 5.100000000000000e+01, 5.200000000000000e+01, 5.300000000000000e+01, 5.400000000000000e+01}; 
    constant Real[9] data_y   = {0.000000000000000e+00, 9.741795621122724e-03, 1.957849382416917e-02, 4.899461037294182e-02, 5.414325847204291e-01, 3.129040153964464e-01, 4.805288278259229e-02, 1.929561728229790e-02, 0.000000000000000e+00}; 
    constant Real[9] m_k      = {9.741795621122724e-03, 9.789246912084585e-03, 6.619164939859322e-03, 8.799976062920696e-02, 0.000000000000000e+00, -2.466898509689184e-01, -8.513907936141946e-02, -1.393413208562771e-02, -1.929561728229790e-02}; 
    constant Real[9] iy_start = {0.000000000000000e+00, 4.855195011125819e-03, 1.974348677517560e-02, 4.718192745533083e-02, 3.489985073335319e-01, 7.956435103582067e-01, 9.622562282899724e-01, 9.899297687144510e-01, 9.999999999999999e-01}; 
    constant Real    iy_scaler = 9.975860568134682e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 7; 
    constant Real[7] data_x   = {4.300000000000000e+01, 4.400000000000000e+01, 4.500000000000000e+01, 4.600000000000000e+01, 4.700000000000000e+01, 4.800000000000000e+01, 4.900000000000000e+01}; 
    constant Real[7] data_y   = {0.000000000000000e+00, 4.822703670653017e-02, 7.187538564589666e-02, 9.529563662679361e-02, 2.245894705013105e-01, 5.600124705194683e-01, 0.000000000000000e+00}; 
    constant Real[7] m_k      = {4.822703670653017e-02, 3.593769282294833e-02, 2.069468215024688e-02, 6.714390169459990e-02, 2.323584169463374e-01, 0.000000000000000e+00, -5.600124705194683e-01}; 
    constant Real[7] iy_start = {0.000000000000000e+00, 2.392495510177146e-02, 8.228818213005962e-02, 1.581573717614557e-01, 2.972803735502134e-01, 6.890852770214579e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.517585697132896e-01; 
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.400000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.300000000000000e+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 6.000000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 6.000000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------
      
annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP50</strong>.<br><br>
  Information taken from: data sheet - last access 02.12.2019.<br><br>
  It also contains the phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  These functions are modelled by piece-wise splines using <strong>pchip</strong> method,
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
  end SP50;
