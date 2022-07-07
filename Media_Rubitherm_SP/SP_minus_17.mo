within slPCMlib.Media_Rubitherm_SP;
package SP_minus_17 "Rubitherm SP-17; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP_minus_17";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.511500000000000e+02, 2.631500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.491500000000000e+02, 2.571500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.954951439970473e+05
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
    constant Integer  len_x    = 13; 
    constant Real[13] data_x   = {-2.200000000000000e+01, -2.100000000000000e+01, -2.000000000000000e+01, -1.900000000000000e+01, -1.800000000000000e+01, -1.700000000000000e+01, -1.600000000000000e+01, -1.500000000000000e+01, -1.400000000000000e+01, -1.300000000000000e+01, -1.200000000000000e+01, -1.100000000000000e+01, -1.000000000000000e+01}; 
    constant Real[13] data_y   = {0.000000000000000e+00, 6.722798809851667e-03, 9.993806167102015e-03, 3.338922930517069e-02, 2.312628560537255e-01, 2.113463917222146e-01, 2.522029054447774e-01, 1.278015080823168e-01, 6.893935054571151e-02, 3.103702632392603e-02, 2.389691614681072e-02, 3.407211398393906e-03, 0.000000000000000e+00}; 
    constant Real[13] m_k      = {6.722798809851667e-03, 3.443737587705410e-03, 5.809418710133686e-03, 6.994542921831175e-02, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, -9.163177744953292e-02, -4.838224087919537e-02, -1.825881621828361e-02, -1.120027635279837e-02, -9.829786366624745e-03, -2.803052908281700e-03}; 
    constant Real[13] iy_start = {0.000000000000000e+00, 3.631771530348665e-03, 1.178646056207950e-02, 2.812034457610943e-02, 1.661655898123739e-01, 3.872946767925822e-01, 6.188854837355595e-01, 8.163669070569304e-01, 9.110580405153874e-01, 9.584982845009903e-01, 9.853557240819868e-01, 9.988828422455913e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.992068086753999e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 9; 
    constant Real[9] data_x   = {-2.400000000000000e+01, -2.300000000000000e+01, -2.200000000000000e+01, -2.100000000000000e+01, -2.000000000000000e+01, -1.900000000000000e+01, -1.800000000000000e+01, -1.700000000000000e+01, -1.600000000000000e+01}; 
    constant Real[9] data_y   = {0.000000000000000e+00, 1.027380959052955e-02, 9.592407323277664e-02, 1.202247974985261e-01, 1.137873133306928e-01, 1.860936568613920e-01, 2.689214305178280e-01, 2.047749189684584e-01, 0.000000000000000e+00}; 
    constant Real[9] m_k      = {6.455721114701123e-03, 3.013775267685557e-02, 5.497549395399826e-02, 0.000000000000000e+00, 0.000000000000000e+00, 7.756705859356756e-02, 0.000000000000000e+00, -1.344607152589140e-01, -2.047749189684584e-01}; 
    constant Real[9] iy_start = {0.000000000000000e+00, 3.108681434319119e-03, 5.325510596032431e-02, 1.639621066048014e-01, 2.789441838711326e-01, 4.199388807452613e-01, 6.498630927777682e-01, 8.936254865864045e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.827019367407412e-01; 
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := -3.000000000000000e+04;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.100000000000000e+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP-17</strong>.<br><br>
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
  end SP_minus_17;
