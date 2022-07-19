within slPCMlib.Media_Rubitherm_RT;
package RT_minus_9HC "Rubitherm RT-9HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "RT-9HC";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.591500000000000e+02, 2.651500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.591500000000000e+02, 2.651500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.612610235396359e+05
             "scalar phase transition enthalpy";
      
    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 2.591500000000000e+02
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
    constant Integer pieces  = 7; 
    constant Integer[7] order  = {1, 6, 6, 5, 5, 6, 1}; 
    constant Real[8] breaks = {-1.140000000000000e+02, -1.400000000000000e+01, -1.300000000000000e+01, -1.200000000000000e+01, -1.000000000000000e+01, -9.000000000000000e+00, -8.000000000000000e+00, 9.200000000000000e+01}; 
    constant Real[7,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.084685888209339e-12, 1.816382978456760e-15, 1.020605866502170e-19, 1.807003620809174e-20, 8.689061394103988e-05, 5.918551431572996e-06},  {9.280916645911526e-05, 4.699643793826760e-04, 9.576844108858105e-04, 9.872771680418545e-04, 5.232313411787943e-04, -2.090709358483553e-03, 1.434182172743347e-03},  {2.374439280208044e-03, 5.591636314037881e-03, 7.664542968398806e-03, 1.085675240278843e-02, 1.158241713991123e-02, -1.673110808087716e-03, 0.000000000000000e+00},  {2.628490313839592e-01, 4.033193208461502e-01, 2.169342040959730e-01, 3.659165719856960e-02, -5.148690940965935e-03, -1.048878930856597e-02, 0.000000000000000e+00},  {9.040567332751200e-01, 8.739239903271089e-01, 1.909291369602265e-01, -8.889099965095384e-02, -5.759263748379578e-02, 7.274140988232279e-02, -2.040762746185455e-02},  {1.874760005848174e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer pieces  = 8; 
    constant Integer[8] order  = {1, 6, 5, 6, 5, 6, 5, 1}; 
    constant Real[9] breaks = {-1.140000000000000e+02, -1.400000000000000e+01, -1.300000000000000e+01, -1.200000000000000e+01, -1.100000000000000e+01, -1.000000000000000e+01, -9.000000000000000e+00, -8.000000000000000e+00, 9.200000000000000e+01}; 
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 3.564828337260037e-13, 7.591788370725007e-16, -3.682260235453612e-19, 0.000000000000000e+00, 2.781505089224387e-04, -4.329733905729787e-05},  {2.348531702223824e-04, 1.130968510638413e-03, 2.132045003365702e-03, 1.915558308078430e-03, 7.412924587527253e-04, -1.650006803998685e-04, 0.000000000000000e+00},  {5.989716770657784e-03, 1.328189987478962e-02, 1.067646787611903e-02, 3.230721339090673e-03, -8.371094324661696e-05, -5.400369112857268e-03, 5.122112435661662e-03},  {3.281683824021488e-02, 4.772298492034562e-02, 4.269436164026222e-02, 5.133443515076477e-02, 4.974613002739198e-02, -1.855395123430744e-02, 0.000000000000000e+00},  {2.057607987446720e-01, 3.933297775912043e-01, 3.096349349138340e-01, 6.477944291725823e-02, -4.302362614414525e-02, -4.564501081745728e-02, 2.212525060593518e-02},  {9.069615678113011e-01, 9.393699211423734e-01, 1.212601577151920e-01, -1.212601577151920e-01, 6.063007885759602e-02, -1.212601577151920e-02, 0.000000000000000e+00},  {1.894835552039751e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT-9HC</strong>.<br><br>
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
  end RT-9HC;
