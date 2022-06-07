within slPCMlib.Media;
package RT21 "Rubitherm RT21, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT21";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {2.871500000000000e+02, 2.981500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {2.861500000000000e+02, 2.961500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 9.876933991753223e+04
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.SIunits.Temp_K            Tref = 273.15+25
             "reference temperature";
    constant Modelica.SIunits.SpecificEnthalpy  href = 0.0
             "reference enthalpy at Tref";

  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  protected
    constant Integer pieces =   data_H.pieces;
    constant Integer order[:] = data_H.order;
    constant Real breaks[:] =   data_H.breaks;
    constant Real coefs[:,:] =  data_H.coefs;
  algorithm
    (Xi, dXi) := BasicUtilities.splineEval(T-273.15,
                 pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
  protected
    constant Integer pieces =   data_C.pieces;
    constant Integer order[:] = data_C.order;
    constant Real breaks[:] =   data_C.breaks;
    constant Real coefs[:,:] =  data_C.coefs;
  algorithm
    (Xi, dXi) := BasicUtilities.splineEval(T-273.15,
                     pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complSolidification;

  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  13;
    constant Integer[13] order =  {1, 6, 5, 5, 5, 5, 6, 5, 5, 5, 5, 6, 1};
    constant Real[14] breaks = {-8.600000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 1.250000000000000e+02};
    constant Real[13,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.776092534681254e-14, -2.023000359407822e-16, 2.268843998587045e-19, 5.782411586589357e-19, 2.588646227674020e-03, -1.045726885566050e-03},  {1.542919342090007e-03, 6.668869824974102e-03, 1.020055899324946e-02, 4.971924565419197e-03, -2.742672145120656e-03, 6.919272232205410e-04, 0.000000000000000e+00},  {2.133352780383265e-02, 3.447470904335843e-02, 1.557957205098860e-02, 9.205082171419817e-04, 7.169639709820485e-04, -1.591042190674741e-04, 0.000000000000000e+00},  {7.286617686723623e-02, 7.046771258535649e-02, 2.105183833763215e-02, 2.197321910395435e-03, -7.855712435532198e-05, 4.581139143745571e-04, 0.000000000000000e+00},  {1.669626064906395e-01, 1.211396960662616e-01, 3.175360046643212e-02, 6.464232556719716e-03, 2.212012447517463e-03, -1.143567559130013e-03, 0.000000000000000e+00},  {3.273885804684404e-01, 2.071698066637094e-01, 5.298269723039598e-02, 3.876606755489421e-03, -3.505825348132603e-03, 5.371172511215415e-03, -1.511673917002537e-03},  {5.917713643641155e-01, 3.285275390524990e-01, 7.461418176518470e-02, 1.333155213506241e-02, 6.749284529064169e-04, -1.367584627779652e-04, 0.000000000000000e+00},  {1.008782807306990e+00, 5.197664804857940e-01, 1.172908242600308e-01, 1.466368131890843e-02, -8.863860983409147e-06, -3.081617450044793e-03, 0.000000000000000e+00},  {1.657413312060695e+00, 7.828956302684322e-01, 1.304125105504078e-01, -1.618794862547314e-02, -1.541695111120737e-02, 5.195533280829803e-03, 0.000000000000000e+00},  {2.544312086423684e+00, 9.594666674522045e-01, 4.130229081504279e-02, -2.590042026200461e-02, 1.056071529294164e-02, -2.067989948142095e-03, 0.000000000000000e+00},  {3.527673349773726e+00, 9.962728997272888e-01, 6.285422305257327e-03, -4.337458571659017e-03, 2.207655522311634e-04, 1.124625129712774e-03, -3.895927467196690e-04},  {4.526850011169837e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 5, 5, 5, 5, 5, 5, 6, 5, 6, 1};
    constant Real[13] breaks = {-8.700000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 1.230000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -4.127035376036931e-14, -8.360091411896434e-16, 3.052657963716548e-18, -5.782411586589357e-19, 1.694475593220665e-03, -7.560275536915509e-04},  {9.384480394870101e-04, 3.936212643910895e-03, 5.604342626832549e-03, 1.824204858375633e-03, -2.868035339269940e-03, 1.196227812050536e-03, 0.000000000000000e+00},  {1.063140064138668e-02, 1.512651017587489e-02, 5.831023286845146e-03, 2.314341621801226e-03, 3.113103720982736e-03, -1.408597366734790e-03, 0.000000000000000e+00},  {3.560778208015589e-02, 3.914100966522904e-02, 1.736669681079741e-02, 6.807828383842728e-04, -3.929883112691213e-03, 1.660708028408967e-03, 0.000000000000000e+00},  {9.052709631028437e-02, 6.850075949325442e-02, 1.243682693389258e-02, 1.568330671709094e-03, 4.373657029353623e-03, -1.355686432556932e-03, 0.000000000000000e+00},  {1.760509840059372e-01, 1.087956013308000e-01, 2.982689679957234e-02, 5.506094463554265e-03, -2.404775133431038e-03, 4.832072510780839e-04, 0.000000000000000e+00},  {3.182580087175109e-01, 1.777646140422734e-01, 3.674860190042974e-02, 7.190664406109518e-04, 1.126112195938122e-05, 1.046386664907303e-03, 0.000000000000000e+00},  {5.345479388876916e-01, 2.586959949773380e-01, 4.943723460309188e-02, 1.122797757752148e-02, 5.243194446495895e-03, -1.094674517641763e-02, 6.826765872392752e-03},  {8.550323611881140e-01, 3.984540440543339e-01, 1.075143703363466e-01, 5.926862104718382e-02, 5.291095665029905e-02, -3.540218671675204e-02, 0.000000000000000e+00},  {1.437778166559525e+00, 8.259215408859718e-01, 2.487641062121712e-01, -8.310941951914039e-02, -1.240999769334612e-01, 1.242128074025110e-01, -3.313093733860627e-02},  {2.396336287268972e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT21</strong>.
  It also contains the phase transition functions for complete melting
  (and solidification), which are modelled by piece-wise splines,
  see
  <blockquote>
  <p>
  Barz, T., Krämer, J., & Emhofer, J. (2020). Identification of Phase
  Fraction–Temperature Curves from Heat Capacity Data for Numerical
  Modeling of Heat Transfer in Commercial Paraffin Waxes.
  Energies, 13(19), 5149.
  <a href>doi.org/10.3390/en13195149</a>.
  </p>
  </blockquote>
  <p>
  </p></html>",
  revisions="<html>
  <ul>
  <li>01-Jun-2022  </ul>
  </html>"));
end RT21;
