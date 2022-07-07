within slPCMlib.Media_Rubitherm_SP;
package SP_minus_11UK "Rubitherm SP-11UK; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP_minus_11UK";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.571500000000000e+02, 2.661500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.551500000000000e+02, 2.631500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.932149386100603e+05
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
    constant Real[10] data_x   = {-1.600000000000000e+01, -1.500000000000000e+01, -1.400000000000000e+01, -1.300000000000000e+01, -1.200000000000000e+01, -1.100000000000000e+01, -1.000000000000000e+01, -9.000000000000000e+00, -8.000000000000000e+00, -7.000000000000000e+00}; 
    constant Real[10] data_y   = {0.000000000000000e+00, 3.398875681582424e-03, 2.037014937769684e-02, 3.722605657490882e-02, 3.564673171668293e-01, 4.899991682389316e-01, 7.540431825361220e-02, 1.371197617793149e-02, 3.422138528428993e-03, 0.000000000000000e+00}; 
    constant Real[10] m_k      = {3.227747618502605e-03, 9.672272142621003e-03, 5.063907247457337e-03, 5.031353009239516e-02, 2.263865558320114e-01, 0.000000000000000e+00, -1.829989009286035e-01, -2.765696830761441e-02, -6.855988088965747e-03, -3.422138528428993e-03}; 
    constant Real[10] iy_start = {0.000000000000000e+00, 1.161750338158111e-03, 1.342349834004777e-02, 3.843693803543129e-02, 2.205099756804495e-01, 6.623639081739776e-01, 9.601505392316054e-01, 9.917460165350618e-01, 9.985758740554460e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.994461497412643e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 9; 
    constant Real[9] data_x   = {-1.800000000000000e+01, -1.700000000000000e+01, -1.600000000000000e+01, -1.500000000000000e+01, -1.400000000000000e+01, -1.300000000000000e+01, -1.200000000000000e+01, -1.100000000000000e+01, -1.000000000000000e+01}; 
    constant Real[9] data_y   = {0.000000000000000e+00, 3.394143752840892e-03, 1.695919855238951e-02, 3.387235535877418e-02, 1.253992690421433e-01, 1.249736458965131e-01, 5.080877114480444e-01, 1.873136759492940e-01, 0.000000000000000e+00}; 
    constant Real[9] m_k      = {3.394143752840892e-03, 8.479599276194755e-03, 1.372890535272361e-02, 4.884681173047185e-02, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, -2.540438557240222e-01, -1.873136759492940e-01}; 
    constant Real[9] iy_start = {0.000000000000000e+00, 1.253365039969011e-03, 1.084023640817388e-02, 3.297770546087964e-02, 1.153746075105360e-01, 2.386026851455428e-01, 5.501816517030529e-01, 9.132821578897569e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.843562964084897e-01; 
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.100000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.200000000000000e+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP-11UK</strong>.<br><br>
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
  end SP_minus_11UK;
