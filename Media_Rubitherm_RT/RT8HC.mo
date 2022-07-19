within slPCMlib.Media_Rubitherm_RT;
package RT8HC "Rubitherm RT8HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "RT8HC";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.761500000000000e+02, 2.851500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.741500000000000e+02, 2.821500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.409814306188960e+05
             "scalar phase transition enthalpy";
      
    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 2.761500000000000e+02
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
    constant Integer pieces  = 11; 
    constant Integer[11] order  = {1, 6, 6, 5, 5, 5, 5, 5, 5, 6, 1}; 
    constant Real[12] breaks = {-9.700000000000000e+01, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.120000000000000e+02}; 
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.074869249398261e-14, 7.624190254558202e-17, -5.495976928597619e-19, 1.445602896647339e-19, 7.345065200705361e-04, -2.522121947603656e-04},  {4.822943253209951e-04, 2.159259431817529e-03, 3.561882279300066e-03, 2.300821305498053e-03, -1.106503210528033e-04, 2.628565770767117e-03, -9.762240326518820e-04},  {1.004594875899907e-02, 2.302837128062723e-02, 2.144274148737035e-02, 8.619397075920363e-03, -1.611181956995451e-03, 9.109065674642099e-04, 0.000000000000000e+00},  {6.243618321338579e-02, 8.988185049250179e-02, 4.674290664780108e-02, 1.128373492258066e-02, 2.943350880325598e-03, 6.242620351042242e-04, 0.000000000000000e+00},  {2.139122881916991e-01, 2.321135822527485e-01, 1.044968370485394e-01, 2.929975879492529e-02, 6.064661055846719e-03, -5.125974059271472e-03, 0.000000000000000e+00},  {5.807611532844875e-01, 5.276353066614055e-01, 1.775243391756794e-01, 2.298662425597460e-03, -1.956520924051064e-02, 4.172412528504814e-03, 0.000000000000000e+00},  {1.272826664835164e+00, 8.321811979701399e-01, 1.087531962944568e-01, -3.423804925139693e-02, 1.296853402013436e-03, 1.426361944983438e-03, 0.000000000000000e+00},  {2.182246225195360e+00, 9.592926661377150e-01, 2.808378840218013e-02, -1.478701619350882e-02, 8.428663126930623e-03, -2.370494539382951e-03, 0.000000000000000e+00},  {3.160893832129295e+00, 9.929613741723469e-01, 1.058977318940784e-02, -4.777309079615866e-03, -3.423809569984128e-03, 4.172240379872063e-03, -1.162492821958412e-03},  {4.159253608399363e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer pieces  = 10; 
    constant Integer[10] order  = {1, 6, 6, 5, 5, 5, 5, 5, 6, 1}; 
    constant Real[11] breaks = {-9.900000000000000e+01, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.090000000000000e+02}; 
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -6.750283657038463e-13, -3.357231572401635e-15, 1.564675699744506e-18, -2.891205793294678e-19, 8.461559504347133e-04, -3.054253047999075e-04},  {5.407306449564214e-04, 2.398227922723941e-03, 3.880179932345320e-03, 2.353053408348983e-03, -3.505998198250471e-04, 3.844811956091218e-04, -6.596975789902851e-05},  {9.140103526259713e-03, 1.734193616382387e-02, 1.169100682604783e-02, 3.476070927159442e-03, 5.822597897351341e-04, -2.142604876433312e-04, 0.000000000000000e+00},  {4.201711674538266e-02, 5.240989931818635e-02, 2.347017346950397e-02, 3.662505209666666e-03, -4.890426484815219e-04, 3.839418530667863e-04, 0.000000000000000e+00},  {1.214545939473249e-01, 1.103013005576506e-01, 3.536285173828294e-02, 5.545753146408441e-03, 1.430666616852410e-03, 1.008103988016194e-03, 0.000000000000000e+00},  {2.751032699945355e-01, 2.084274498809608e-01, 7.066515075878479e-02, 2.134945949398001e-02, 6.471186556933377e-03, -1.153622441013365e-04, 0.000000000000000e+00},  {5.819011544410930e-01, 4.391140648878200e-01, 1.723870261413123e-01, 4.608058328070017e-02, 5.894375336426695e-03, -1.521255100926709e-02, 0.000000000000000e+00},  {1.230164653078085e+00, 8.696446133113409e-01, 1.938695179092993e-01, -8.246742546626398e-02, -7.016837970990876e-02, 8.087493140780622e-02, -2.228041848860815e-02},  {2.199637492041751e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT8HC</strong>.<br><br>
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
  end RT8HC;
