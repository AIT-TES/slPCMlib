within slPCMlib.Media_Rubitherm_RT;
package RT70HC "Rubitherm RT70HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "RT70HC";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.401500000000000e+02, 3.451500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.371500000000000e+02, 3.441500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.123177848153893e+05
             "scalar phase transition enthalpy";
      
    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 3.401500000000000e+02
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
    constant Real[8] breaks = {-3.300000000000000e+01, 6.700000000000000e+01, 6.800000000000000e+01, 6.900000000000000e+01, 7.000000000000000e+01, 7.100000000000000e+01, 7.200000000000000e+01, 1.720000000000000e+02}; 
    constant Real[7,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 2.900998305326960e-14, 3.415868232434442e-17, 2.491560366366319e-18, -5.782411586589357e-19, 1.627145685097248e-03, -4.301720549541999e-04},  {1.196973630172094e-03, 5.554696095786461e-03, 9.818876026659519e-03, 7.668015751888456e-03, 1.683147601173241e-03, -1.372799732043600e-02, 9.460037184894098e-03},  {2.165374897013787e-02, 4.304932231664894e-02, 4.754239345841579e-02, 6.632137665010338e-02, 7.494371877240470e-02, -3.490527426685485e-02, 0.000000000000000e+00},  {2.186052859008558e-01, 4.623467429390474e-01, 3.471160933746056e-01, 1.704350907117371e-02, -9.958265256186953e-02, 2.902541997031388e-02, 0.000000000000000e+00},  {9.745543986941270e-01, 9.545059465058925e-01, 9.100490492004842e-02, -9.103290147316556e-02, 4.554444728969988e-02, -9.125687389810237e-03, 5.599310623419742e-06},  {1.965456707857415e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer pieces  = 9; 
    constant Integer[9] order  = {1, 6, 5, 5, 5, 5, 6, 6, 1}; 
    constant Real[10] breaks = {-3.600000000000000e+01, 6.400000000000000e+01, 6.500000000000000e+01, 6.600000000000000e+01, 6.700000000000000e+01, 6.800000000000000e+01, 6.900000000000000e+01, 7.000000000000000e+01, 7.100000000000000e+01, 1.710000000000000e+02}; 
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 3.772328462862349e-15, 1.020456514258352e-17, 2.380889872394849e-18, 0.000000000000000e+00, 3.117819790861473e-03, -1.446107492332027e-03},  {1.671712298533230e-03, 6.912454000319087e-03, 9.486585523634334e-03, 2.256048061974178e-03, -6.102513430673050e-03, 5.833309605764335e-03, 0.000000000000000e+00},  {2.005759605955212e-02, 3.741026353964082e-02, 3.797274518316192e-02, 3.617909039692532e-02, 2.306403459814862e-02, -1.384261916366290e-02, 0.000000000000000e+00},  {1.408411106137659e-01, 2.449360676710284e-01, 1.464680323262006e-01, -9.990962847109154e-03, -4.614906122016588e-02, 1.681287178328673e-02, 0.000000000000000e+00},  {4.929180583270066e-01, 4.073673578178773e-01, 7.729494296745271e-03, -2.645848989490534e-02, 3.791529769626780e-02, -1.253255393017824e-02, 0.000000000000000e+00},  {9.069391643128134e-01, 4.324492978608366e-01, 3.052027148785370e-02, -1.228384116161820e-04, -2.474747195462338e-02, 1.533922464389481e-01, -7.183362927323197e-02},  {1.426597040460980e+00, 7.300908943385266e-01, 3.380849498162656e-01, -1.862847305268665e-03, -3.352906788583628e-01, 2.687913972782709e-01, -6.724442050219943e-02},  {2.359166335228212e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT70HC</strong>.<br><br>
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
  end RT70HC;
