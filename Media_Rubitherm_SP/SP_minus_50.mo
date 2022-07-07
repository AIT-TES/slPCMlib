within slPCMlib.Media_Rubitherm_SP;
package SP_minus_50 "Rubitherm SP-50; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP_minus_50";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.161500000000000e+02, 2.251500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.141500000000000e+02, 2.231500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.900000000000000e+05
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
    constant Real[10] data_x   = {-5.700000000000000e+01, -5.600000000000000e+01, -5.500000000000000e+01, -5.400000000000000e+01, -5.300000000000000e+01, -5.200000000000000e+01, -5.100000000000000e+01, -5.000000000000000e+01, -4.900000000000000e+01, -4.800000000000000e+01}; 
    constant Real[10] data_y   = {0.000000000000000e+00, 5.244310786019246e-03, 1.046111877641813e-02, 1.565056820440438e-02, 6.801159962095515e-02, 3.193818433844124e-01, 3.649057028124660e-01, 1.794411463886531e-01, 3.690371002667172e-02, 0.000000000000000e+00}; 
    constant Real[10] m_k      = {5.244310786019246e-03, 5.230559388209067e-03, 2.770141499422345e-03, 1.531991463304607e-02, 1.518656375900040e-01, 1.365715782841609e-01, 0.000000000000000e+00, -1.640009963928972e-01, -8.972057319432653e-02, -3.690371002667172e-02}; 
    constant Real[10] iy_start = {0.000000000000000e+00, 2.614119678608732e-03, 1.064366682387671e-02, 2.261166028336793e-02, 5.295734972047034e-02, 2.472461687330232e-01, 5.995335553465562e-01, 8.843732771031174e-01, 9.859987274113414e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.964999582516468e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 10; 
    constant Real[10] data_x   = {-5.900000000000000e+01, -5.800000000000000e+01, -5.700000000000000e+01, -5.600000000000000e+01, -5.500000000000000e+01, -5.400000000000000e+01, -5.300000000000000e+01, -5.200000000000000e+01, -5.100000000000000e+01, -5.000000000000000e+01}; 
    constant Real[10] data_y   = {0.000000000000000e+00, 1.242236024844720e-02, 3.105590062111801e-02, 8.074534161490683e-02, 9.937888198757763e-02, 3.167701863354037e-01, 2.981366459627329e-01, 8.695652173913043e-02, 7.453416149068323e-02, 0.000000000000000e+00}; 
    constant Real[10] m_k      = {1.242236024844720e-02, 1.552795031055901e-02, 3.416149068322981e-02, 1.554361918552493e-02, 5.369613900454068e-02, 0.000000000000000e+00, -5.590062111801236e-02, -3.473309334936506e-02, -1.350731408030863e-02, -7.453416149068323e-02}; 
    constant Real[10] iy_start = {0.000000000000000e+00, 5.909558067831450e-03, 2.595066803699898e-02, 8.298945391116677e-02, 1.692471334122657e-01, 3.802672147995889e-01, 6.901336073997945e-01, 8.795436937457595e-01, 9.579520439706731e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.928057553956836e-01; 
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.100000000000000e+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP-50</strong>.<br><br>
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
  end SP_minus_50;
