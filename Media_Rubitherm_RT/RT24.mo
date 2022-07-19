within slPCMlib.Media_Rubitherm_RT;
package RT24 "Rubitherm RT24; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "RT24";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.891500000000000e+02, 2.991500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.891500000000000e+02, 2.981500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 9.780000975944677e+04
             "scalar phase transition enthalpy";
      
    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 2.891500000000000e+02
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
    constant Integer pieces  = 12; 
    constant Integer[12] order  = {1, 6, 5, 5, 5, 6, 5, 5, 5, 5, 6, 1}; 
    constant Real[13] breaks = {-8.400000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 1.260000000000000e+02}; 
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -8.928495384073478e-16, -2.206257160019698e-17, 1.075840157084096e-18, -5.782411586589357e-19, 2.620404937429120e-03, -1.054108711809066e-03},  {1.566296225619140e-03, 6.777372416291734e-03, 1.039241869715523e-02, 5.121875138109886e-03, -2.709605989990386e-03, 6.712891034222948e-04, 0.000000000000000e+00},  {2.181964559060789e-02, 3.544585678208340e-02, 1.621329920576555e-02, 9.963422123712869e-04, 6.468395271210869e-04, -8.287646203543920e-05, 0.000000000000000e+00},  {7.503910685591379e-02, 7.303445762903690e-02, 2.225459838525158e-02, 2.754935700501242e-03, 2.324572169438910e-04, 7.092538033304659e-04, 0.000000000000000e+00},  {1.740248095909779e-01, 1.302845593854718e-01, 3.900668682172331e-02, 1.077730260158146e-02, 3.778726233596220e-03, -2.722816343477104e-03, 3.812679330735790e-04},  {3.555305362229471e-01, 2.444182716491122e-01, 7.250180758937788e-02, 6.289402762666884e-03, -4.116336487685615e-03, 1.994362174095247e-03, 0.000000000000000e+00},  {6.766180439105136e-01, 4.017965600356033e-01, 8.661561869221733e-02, 9.767678552876892e-03, 5.855474382790619e-03, -2.853008913917942e-03, 0.000000000000000e+00},  {1.177800366660084e+00, 6.134876860402386e-01, 1.225214115084122e-01, 4.659486944859950e-03, -8.409570186799090e-03, 1.151567406216715e-03, 0.000000000000000e+00},  {1.911210948373012e+00, 8.446285261755349e-01, 9.755812528436482e-02, -1.746311974016925e-02, -2.651733155715512e-03, 5.166148190232893e-04, 0.000000000000000e+00},  {2.833799361756051e+00, 9.793315589960018e-01, 3.442451531979666e-02, -2.290390417279841e-02, -6.865906059906612e-05, 6.926098500318775e-03, -2.304122229399654e-03},  {3.829204849109371e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer pieces  = 11; 
    constant Integer[11] order  = {1, 6, 5, 5, 5, 5, 6, 5, 5, 6, 1}; 
    constant Real[12] breaks = {-8.400000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 1.250000000000000e+02}; 
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 9.207430603867839e-14, 1.887863436945520e-15, -3.653568527904668e-18, -5.782411586589357e-19, 2.515859465054182e-03, -1.008426409062838e-03},  {1.507433056085302e-03, 6.528738870990856e-03, 1.003219851460114e-02, 4.990066469285048e-03, -2.547098810671665e-03, 5.830776224382072e-04, 0.000000000000000e+00},  {2.109441572272888e-02, 3.429032817755095e-02, 1.555058128280833e-02, 6.324474509804580e-04, 3.682893015193702e-04, 1.564929305039414e-04, 0.000000000000000e+00},  {7.209255486609194e-02, 6.954445495470542e-02, 2.122258874990532e-02, 3.670533962097352e-03, 1.150753954039077e-03, -1.897636720541952e-04, 0.000000000000000e+00},  {1.674911228147849e-01, 1.266554317966954e-01, 3.724107763988993e-02, 6.375913057711707e-03, 2.019355937681010e-04, -3.592772038808242e-04, 0.000000000000000e+00},  {3.376062036989692e-01, 2.192766826052774e-01, 5.398765833682539e-02, 3.590883393975873e-03, -1.594450425636020e-03, 2.482565918740327e-03, -5.797946091459607e-04},  {6.147697489190062e-01, 3.405809096971384e-01, 7.132234601515076e-02, 1.044284869591585e-02, 2.121460030876203e-03, 1.927135741351502e-03, 0.000000000000000e+00},  {1.041164449099439e+00, 5.326756666454562e-01, 1.346510097016707e-01, 3.820004623293569e-02, 1.175713873763371e-02, -1.531313404463731e-02, 0.000000000000000e+00},  {1.743135176372498e+00, 8.870407094749799e-01, 1.666626403799075e-01, -6.790273926290256e-02, -6.480853148555284e-02, 7.221764696731303e-02, -1.975198022340082e-02},  {2.716592922222842e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}}; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT24</strong>.<br><br>
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
  end RT24;
