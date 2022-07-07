within slPCMlib.Media_Rubitherm_SP;
package SP58 "Rubitherm SP58; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP58";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.251500000000000e+02, 3.341500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.241500000000000e+02, 3.301500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.500000000000000e+05
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
    constant Real[10] data_x   = {5.200000000000000e+01, 5.300000000000000e+01, 5.400000000000000e+01, 5.500000000000000e+01, 5.600000000000000e+01, 5.700000000000000e+01, 5.800000000000000e+01, 5.900000000000000e+01, 6.000000000000000e+01, 6.100000000000000e+01}; 
    constant Real[10] data_y   = {0.000000000000000e+00, 1.053954941069164e-02, 2.124572197405386e-02, 5.847002541624152e-02, 1.331712428739610e-01, 2.564814084384445e-01, 2.868850796809143e-01, 2.176436930148033e-01, 1.556327919099859e-02, 0.000000000000000e+00}; 
    constant Real[10] m_k      = {1.053954941069164e-02, 1.062286098702693e-02, 2.396523800277494e-02, 5.596276044995359e-02, 9.900569151110147e-02, 7.685691840347661e-02, 0.000000000000000e+00, -1.356609002449578e-01, -4.621955174876617e-02, -6.610141355209676e-03}; 
    constant Real[10] iy_start = {0.000000000000000e+00, 5.255321479130435e-03, 2.001499876175579e-02, 5.715333633588823e-02, 1.492554326866083e-01, 3.456468184738137e-01, 6.233379459866460e-01, 8.865312670159475e-01, 9.955255392388309e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.985728986316630e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 7; 
    constant Real[7] data_x   = {5.100000000000000e+01, 5.200000000000000e+01, 5.300000000000000e+01, 5.400000000000000e+01, 5.500000000000000e+01, 5.600000000000000e+01, 5.700000000000000e+01}; 
    constant Real[7] data_y   = {0.000000000000000e+00, 1.005246103969156e-02, 3.025843509202920e-02, 1.511921393082126e-01, 5.447579729460870e-01, 2.637389916139793e-01, 0.000000000000000e+00}; 
    constant Real[7] m_k      = {1.005246103969156e-02, 1.270692672791929e-02, 5.927112703302484e-02, 2.572497689270289e-01, 0.000000000000000e+00, -2.723789864730435e-01, -2.637389916139793e-01}; 
    constant Real[7] iy_start = {0.000000000000000e+00, 4.697839357328752e-03, 2.060988880338592e-02, 9.318117274147716e-02, 4.543532245522806e-01, 8.717760554352136e-01, 9.999999999999999e-01}; 
    constant Real    iy_scaler = 9.776930010820354e-01; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP58</strong>.<br><br>
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
  end SP58;
