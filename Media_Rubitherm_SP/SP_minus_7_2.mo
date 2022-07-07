within slPCMlib.Media_Rubitherm_SP;
package SP_minus_7_2 "Rubitherm SP-7_2; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP_minus_7_2";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.641500000000000e+02, 2.711500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.631500000000000e+02, 2.691500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.867197414377566e+05
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
    constant Integer  len_x    = 8; 
    constant Real[8] data_x   = {-9.000000000000000e+00, -8.000000000000000e+00, -7.000000000000000e+00, -6.000000000000000e+00, -5.000000000000000e+00, -4.000000000000000e+00, -3.000000000000000e+00, -2.000000000000000e+00}; 
    constant Real[8] data_y   = {0.000000000000000e+00, 2.699496263836573e-02, 1.354120499965367e-01, 4.993378152889021e-01, 2.482739373611794e-01, 6.872947881817944e-02, 2.125175589683180e-02, 0.000000000000000e+00}; 
    constant Real[8] m_k      = {2.699496263836573e-02, 6.770602499826835e-02, 2.361714263252682e-01, 0.000000000000000e+00, -2.153041682353613e-01, -1.135110907321738e-01, -3.436473940908972e-02, -2.125175589683180e-02}; 
    constant Real[8] iy_start = {0.000000000000000e+00, 1.006442815316408e-02, 7.696019260015971e-02, 4.126663507265151e-01, 8.028455024286316e-01, 9.522637092124849e-01, 9.905050457098672e-01, 9.999999999999999e-01}; 
    constant Real    iy_scaler = 9.959955402921060e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 7; 
    constant Real[7] data_x   = {-1.000000000000000e+01, -9.000000000000000e+00, -8.000000000000000e+00, -7.000000000000000e+00, -6.000000000000000e+00, -5.000000000000000e+00, -4.000000000000000e+00}; 
    constant Real[7] data_y   = {0.000000000000000e+00, 3.798200666469613e-02, 9.150994787817839e-02, 3.073624474343970e-01, 5.396119691440867e-01, 2.353362887863820e-02, 0.000000000000000e+00}; 
    constant Real[7] m_k      = {3.798200666469613e-02, 4.575497393908919e-02, 1.346902203848504e-01, 2.240510106329542e-01, 0.000000000000000e+00, -7.033384111691308e-02, -6.134817643426157e-03}; 
    constant Real[7] iy_start = {0.000000000000000e+00, 1.827606589380864e-02, 7.540075928958295e-02, 2.666869793781783e-01, 7.072255084793915e-01, 9.936066088596109e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.963370643981931e-01; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP-7_2</strong>.<br><br>
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
  end SP_minus_7_2;
