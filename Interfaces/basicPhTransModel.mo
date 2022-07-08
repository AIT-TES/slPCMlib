within slPCMlib.Interfaces;
partial model basicPhTransModel
  "Basic phase transition model, computes properties for given T and xi"

    replaceable package PCM =
      slPCMlib.Media_generic.generic_7thOrderSmoothStep;
//             constrainedby slPCMlib.Interfaces.partialPCM     "PCM"
//             annotation (Dialog(group="PCM"),
//                         choicesAllMatching=true);

  Modelica.Units.SI.MassFraction xi(
    min=0,
    max=1,
    start=0.5) "liquid mass phase fraction";
  Modelica.Units.SI.VolumeFraction phi(
    min=0,
    max=1,
    start=0.5) "liquid volume phase fraction";
  Modelica.Units.SI.Density rho(start=1e3) "effective density";
  Modelica.Units.SI.SpecificHeatCapacity cp(start=1e3)
    "effective specific heat capacity";
  Modelica.Units.SI.SpecificHeatCapacity cp_BL(start=1e3)
    "baseline specific heat capacity";
  Modelica.Units.SI.SpecificEnthalpy h(start=1e3) "specific enthalpy";
  Modelica.Units.SI.SpecificEnthalpy h_S(start=1e3)
    "solid specific enthalpy (extrapolated)";
  Modelica.Units.SI.SpecificEnthalpy h_L(start=1e3)
    "liquid specific enthalpy (extrapolated)";
  Modelica.Units.SI.ThermalConductivity lambda(start=0.2)
    "effective thermal conductivity";

    // "Connector for temperature (input signal)"
    inductionAtNode indVar;

//     Modelica.SIunits.Temp_K T(start=273.15+20) = indVar.T
//       "Temperature";
//     Modelica.SIunits.TemperatureSlope der_T(start=0) = indVar.der_T
//       "Time derivative of temperature (= der(T))";

     package enthFcts =
        slPCMlib.BasicUtilities.enthalpyHelpers(redeclare package PCM = PCM);

    Real dxi(unit="1/K");
protected
  parameter Modelica.Units.SI.SpecificEnthalpy h_L_at_Tmax(start=1e3, fixed=
        false) "phase transition enthalpy at T_max, where melting is finished";
  parameter Modelica.Units.SI.SpecificEnthalpy h_BL_at_Tmax(start=1e3, fixed=
        false) "baseline enthalpy at T_max, where melting is finished";
  parameter Modelica.Units.SI.SpecificEnthalpy h_offset(start=1e3, fixed=false)
    "offset für h_L";

initial equation
//   assert(PCM.propData.rangeTmelting[1] >= PCM.propData.rangeTsolidification[1],
//          "PCM.propData.rangeTmelting[1] <= PCM.propData.rangeTsolidification[1].
//           Phase transition function for complete melting should give always smaller values
//           compared with the function for solidification!",
//          AssertionLevel.error);
//    assert(PCM.propData.rangeTmelting[2] >= PCM.propData.rangeTsolidification[2],
//           "PCM.propData.rangeTmelting[2] <= PCM.propData.rangeTsolidification[2].
//            Phase transition function for complete melting should give always smaller values
//            compared with the function for solidification!",
//           AssertionLevel.error);

  h_BL_at_Tmax =Modelica.Math.Nonlinear.quadratureLobatto(
  function BasicUtilities.enthalpyHelpers.spHeatCap_baselineMelting(),
  PCM.propData.Tref,
  PCM.propData.rangeTmelting[2],
  tolerance=100*Modelica.Constants.eps);

  h_L_at_Tmax = h_BL_at_Tmax + PCM.propData.phTrEnth; // pcmData.phTrEnth;
  h_offset    = - enthFcts.enthalpy_liquid(PCM.propData.rangeTmelting[2])
              + h_L_at_Tmax;
equation

    h       = xi*(enthFcts.enthalpy_liquid(indVar.T) + h_offset)
            + (1.0 - xi)*enthFcts.enthalpy_solid(indVar.T);
    h_S     = enthFcts.enthalpy_solid(indVar.T);
    h_L     = enthFcts.enthalpy_liquid(indVar.T) + h_offset;
    cp_BL   =  xi*enthFcts.spHeatCap_liquid(indVar.T)
            + (1.0 - xi)*enthFcts.spHeatCap_solid(indVar.T);
    cp      = cp_BL + dxi*(h_L - h_S);
    phi    = xi/(xi + (1 - xi)*PCM.density_liquid(indVar.T)/PCM.density_solid(indVar.T));
    rho     = phi*PCM.density_liquid(indVar.T) + (1.0 - phi)*PCM.density_solid(indVar.T);
    lambda  = phi*PCM.conductivity_liquid(indVar.T) + (1.0 - phi)*PCM.conductivity_solid(indVar.T);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false),
              graphics={Text(
              lineColor={0,0,255},
              extent={{-150,105},{150,145}},
              textString="%name"), Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={253,212,198},
              radius=25,
              fillColor={235,238,219},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5,
              pattern=LinePattern.None)}),
              Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
          <p>
Assumptions for modeling effective PCM properties:
<ul>

  <li>There are only two phases (two-phase model): a
solid and a liquid phase.</li> 
  <li>Phase transitions are induced by temperature and are
independent of pressure.</li> 
  <li>Phase transitions extend over a temperature range
(non-isothermal phase transitions) and are continuous.</li> 
  <li>Within the phase transition temperature range the
  solid and liquid phases coexist as a homogenous mixture (macroscopic view). 
  The material is then in a semi-solid or semi-liquid state which produces a
mushy zone in the PCM domain.</li> 
  <li>Properties of the mushy state are local
effective (also apparent) mixture properties, which
are defined by a weighting of contributions from
solid and liquid phases. The weighting is based on
the phase change progress, i.e. the mass (or volume) phase fraction.  </li> 
        </ul>
        </p>
      <p>
      Temperature is the input to the model using the 
      &lt; inductionAtNode &gt; connector. 
      For a given temperature input <var> T </var>  
      the liquid mass phase fraction <var> xi </var> is 
      computed, 
      and the following variables are derived:
      </p>  
      <ul>
      <li>liquid volume phase fractions <var> phi </var> </li> 
      <li>effective density <var> rho </var> </li>
      <li>effective thermal conductivity <var> lambda </var> </li>  
      <li>effective specific heat capacity <var> cp </var></li>   
      <li>and specific enthalpy <var> h </var></li>
      <li>baseline heat capacity <var> c_BL </var>, which 
      describes the mixture heat capacity (without the effect of the 
      phase transition enthalpy)</li>
      <li>solid <var> h_S </var> and liquid <var> h_L </var> enthalpies, 
      where the difference <var> h_L - h_S </var> 
      is the temperature-dependent phase transition enthalpy. 
      </ul>
      </p> 
      <p>
      The <a href>slPCMlib.Interfaces.basicPcTransModel</a>  is extended by 
      specific phase transition models which compute <var> xi </var> 
      as function of temperature <var> T </var> and are also 
      contained in 
      <a href>slPCMlib.Interfaces</a>.
      </p>
                </html>",
     revisions="<html>
       <ul>
       <li>2022-06-01; initial version; by Tilman Barz </li>
       </ul>
       </html>"));
end basicPhTransModel;
