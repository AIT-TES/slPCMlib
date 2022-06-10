﻿within slPCMlib.Interfaces;
partial model basicPhTransModel
  "Basic phase transition model, computes properties for given T and xi_m"

    replaceable package PCM =
      slPCMlib.Media.generic_7thOrderSmoothStep
            constrainedby slPCMlib.Interfaces.partialPCM     "PCM"
            annotation (Dialog(group="PCM"),
                        choicesAllMatching=true);

    Modelica.SIunits.MassFraction          xi_m(min=0,max=1,start=0.5)
      "liquid mass phase fraction";
    Modelica.SIunits.VolumeFraction        xi_v(min=0,max=1,start=0.5)
      "liquid volume phase fraction";
    Modelica.SIunits.Density               rho(start=1e3)      "effective density";
    Modelica.SIunits.SpecificHeatCapacity  cp(start=1e3)       "effective specific heat capacity";
    Modelica.SIunits.SpecificHeatCapacity  cp_BL(start=1e3)    "baseline specific heat capacity";
    Modelica.SIunits.SpecificEnthalpy      h(start=1e3)        "specific enthalpy";
    Modelica.SIunits.SpecificEnthalpy      h_S(start=1e3)      "solid specific enthalpy (extrapolated)";
    Modelica.SIunits.SpecificEnthalpy      h_L(start=1e3)      "liquid specific enthalpy (extrapolated)";
    Modelica.SIunits.ThermalConductivity   lambda(start=0.2)   "effective thermal conductivity";

    // "Connector for temperature (input signal)"
    inductionAtNode indVar;

//     Modelica.SIunits.Temp_K T(start=273.15+20) = indVar.T
//       "Temperature";
//     Modelica.SIunits.TemperatureSlope der_T(start=0) = indVar.der_T
//       "Time derivative of temperature (= der(T))";

     package enthFcts =
        slPCMlib.BasicUtilities.enthalpyHelpers(redeclare package PCM = PCM);

    Real dxi_m(unit="1/K");
protected
    parameter Modelica.SIunits.SpecificEnthalpy    h_L_at_Tmax(start=1e3,fixed=false)
    "phase transition enthalpy at T_max, where melting is finished";
    parameter Modelica.SIunits.SpecificEnthalpy    h_BL_at_Tmax(start=1e3,fixed=false)
    "baseline enthalpy at T_max, where melting is finished";
    parameter Modelica.SIunits.SpecificEnthalpy    h_offset(start=1e3,fixed=false)
    "offset für h_L";


initial equation

    h_BL_at_Tmax =Modelica.Math.Nonlinear.quadratureLobatto(
    function BasicUtilities.enthalpyHelpers.spHeatCap_baselineMelting(),
    PCM.propData.Tref,
    PCM.propData.rangeTmelting[2],
    tolerance=100*Modelica.Constants.eps);

    h_L_at_Tmax = h_BL_at_Tmax + PCM.propData.phTrEnth; // pcmData.phTrEnth;
    h_offset    = - enthFcts.enthalpy_liquid(PCM.propData.rangeTmelting[2])
                  + h_L_at_Tmax;
equation

    h       = xi_m*(enthFcts.enthalpy_liquid(indVar.T) + h_offset)
            + (1.0 - xi_m)*enthFcts.enthalpy_solid(indVar.T);
    h_S     = enthFcts.enthalpy_solid(indVar.T);
    h_L     = enthFcts.enthalpy_liquid(indVar.T) + h_offset;
    cp_BL   =  xi_m*enthFcts.spHeatCap_liquid(indVar.T)
            + (1.0 - xi_m)*enthFcts.spHeatCap_solid(indVar.T);
    cp      = cp_BL + dxi_m*(h_L - h_S);
    xi_v    = xi_m/(xi_m + (1 - xi_m)*PCM.density_liquid(indVar.T)/PCM.density_solid(indVar.T));
    rho     = xi_v*PCM.density_liquid(indVar.T) + (1.0 - xi_v)*PCM.density_solid(indVar.T);
    lambda  = xi_v*PCM.conductivity_liquid(indVar.T) + (1.0 - xi_v)*PCM.conductivity_solid(indVar.T);

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
          Two phases are considered, a solid and a liquid phase. 
          It is assumed that these phases co-exist as homogenous mixture 
          (macroscopic view) in an extended phase transition temperature range. 
          Effective (also apparent) mixture properties are computed by a weighting of 
          contributions from solid and liquid phase. 
          The weighting factor is the mass (or volume) liquid phase fraction.  
          </p>
          <p>
          It is assumed that phase transitions are pressure independent,  
          and are induced by temperature. 
          Temperature is input to the model using the 
          <var> inductionAtNode </var> connector. 
          </p>
          <p>
          For a given temperature input signal the following variables 
          are computed: 
          </p>  
      <ul>
          <li>liquid mass <var> xi_m </var> and volume <var> xi_v </var> 
              phase fractions </li>         
          <li>effective density <var> rho </var> </li>
          <li>effective thermal conductivity <var> lambda </var> </li>  
          <li>effective specific heat capacity <var> cp </var></li>   
          <li>and specific enthalpy <var> h </var></li>
          </ul>
          </p>
          <p>
          Furthermore:
      <ul>
          <li>baseline heat capacity <var> c_BL </var>, which 
          describes the mixture heat capacity (without the effect of the 
          phase transition enthalpy)</li>
          <li>solid <var> h_S </var> and liquid <var> h_L </var> enthalpies, 
          where the difference <var> h_L - h_S </var> 
          is the temperature-dependent phase transition enthalpy. 
          </ul>
          </p>      
          </html>",
     revisions="<html>
       <ul>
       <li>2022-06-01; initial version; by Tilman Barz </li>
       </ul>
       </html>"));
end basicPhTransModel;