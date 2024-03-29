within slPCMlib.Interfaces;
partial package partialPCM "Partial PCM model"
  extends Modelica.Icons.MaterialPropertiesPackage;

  replaceable record propData "basic PCM record"

    constant String mediumName = "";

    // --- parameters for phase transition functions ---
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={0,0}
      "temperature range solidification {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={0,0}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={0.0,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={0.0,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=0.0
      "scalar phase transition enthalpy";

    constant Modelica.Units.SI.Temperature Tref=0.0 "reference Temperature";
    constant Modelica.Units.SI.SpecificEnthalpy href=0.0
      "reference enthalpy at Tref";

  end propData;
  // -------------------------------------------------------
//protected
  replaceable partial function phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.MassFraction xi;
    output Real dxi(unit="1/K");
  end phaseFrac_complMelting;

  replaceable partial function phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.MassFraction xi;
    output Real dxi(unit="1/K");
  end phaseFrac_complSolidification;

  replaceable partial function density_solid "Returns solid density"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.Density rho;
  end density_solid;

  replaceable partial function density_liquid "Returns liquid density"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.Density rho;
  end density_liquid;

  replaceable partial function conductivity_solid "Returns solid thermal conductivity"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.ThermalConductivity lambda;
  end conductivity_solid;

  replaceable partial function conductivity_liquid "Returns liquid thermal conductivity"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.ThermalConductivity lambda;
  end conductivity_liquid;
// ------------------------------------------------------------

  annotation (Icon(graphics,
      coordinateSystem(preserveAspectRatio=false)), Diagram(graphics,
      coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
          <p>
          Package containing partial PCM models. 
          </p>
          <p>
          The record &lt; propData &gt; holds information on: 
          <ul>
          <li>phase transition temperature ranges 
          &lt; rangeTmeltingsolidification[2]&gt;, 
          &lt; rangeTsolidification[2]&gt;</li> 
          <li>coefficients of linear models of the 
          solid and liquid specific heat capacities, 
          and the reference temperature: </li> 
          <pre>
          cp_S(T) = cpS_linCoef[1] + cpS_linCoef[2]*(T-Tref) 
          cp_L(T) = cpL_linCoef[1] + cpL_linCoef[2]*(T-Tref)  </pre>
          <li>reference temperature &lt; Tref &gt; and enthalpy &lt; href = h(Tref) &gt;</li> 
          </ul>
          <p>
          The package also contains functions for the computation of 
          liquid mass phase fraction for melting and solidification processes,  
          functions for solid and liquid density and thermal 
          conductivity. 
          The above functions can be arbitrary nonlinear functions. 
          In contrast, the solid and liquid specific heat capacities are 
          defined by the linear coefficients above. 
          The solid, liquid and mixture enthalpy are then derived in 
          <a href>slPCMlib.Interfaces.partialPCM</a>. 
          </p>          
          </html>",
     revisions="<html>
       <ul>
       <li>2022-06-01; initial version; by Tilman Barz </li>
       </ul>
       </html>"));
end partialPCM;
