within slPCMlib.Media_Rubitherm_SP;
package SP29Eu "Rubitherm SP29Eu; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP29Eu";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.991500000000000e+02, 3.081500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.961500000000000e+02, 3.031500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.000000000000000e+05
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
    constant Integer  len_x    = 10; 
    constant Real[10] data_x   = {2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01, 2.900000000000000e+01, 3.000000000000000e+01, 3.100000000000000e+01, 3.200000000000000e+01, 3.300000000000000e+01, 3.400000000000000e+01, 3.500000000000000e+01}; 
    constant Real[10] data_y   = {0.000000000000000e+00, 6.628130454433402e-03, 5.956530986325864e-02, 6.818199592913335e-01, 6.962170466465549e-02, 1.416082161814685e-01, 2.042467254112671e-02, 6.897656844076670e-03, 1.343435015972519e-02, 0.000000000000000e+00}; 
    constant Real[10] m_k      = {4.319592909768242e-03, 9.115562167936681e-03, 1.585497120786469e-01, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, -4.058104709115011e-02, 0.000000000000000e+00, 0.000000000000000e+00, -1.343435015972519e-02}; 
    constant Real[10] iy_start = {0.000000000000000e+00, 2.910095649661532e-03, 2.352347258922551e-02, 4.068614364741133e-01, 7.820272124663629e-01, 8.874861468962468e-01, 9.717596629146623e-01, 9.820238878109623e-01, 9.921748729784655e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.985226904166510e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 8; 
    constant Real[8] data_x   = {2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01, 2.900000000000000e+01, 3.000000000000000e+01}; 
    constant Real[8] data_y   = {0.000000000000000e+00, 6.896551724137931e-03, 3.448275862068965e-02, 3.793103448275862e-01, 4.000000000000000e-01, 1.448275862068966e-01, 3.448275862068965e-02, 0.000000000000000e+00}; 
    constant Real[8] m_k      = {6.896551724137931e-03, 7.630196736060922e-03, 8.240612474945797e-02, 6.206896551724150e-02, 0.000000000000000e+00, -1.827586206896552e-01, -7.241379310344828e-02, -3.448275862068965e-02}; 
    constant Real[8] iy_start = {0.000000000000000e+00, 3.375499125584861e-03, 1.778414198816071e-02, 2.256586483390607e-01, 6.191294387170676e-01, 9.057846506300116e-01, 9.859679266895763e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.965635738831616e-01; 
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.550000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.500000000000000e+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.000000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 5.000000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------
      
annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP29Eu</strong>.<br><br>
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
  end SP29Eu;
