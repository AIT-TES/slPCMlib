within slPCMlib.Components;
model SingleLayerSlPCMlib
  "Model for single layer heat conductance in a PCM"
  extends Buildings.HeatTransfer.Conduction.BaseClasses.PartialConductor(
   final R=thickness/k/A);
   //final R=if (material.R < Modelica.Constants.eps) then thickness/k/A else material.R/A);
   // if material.R == 0, then the material specifies material.k, and this model specifies x
   // For resistances, material.k need not be specified, and hence we use material.R

  // The value T[:].start is used by the solver when finding initial states
  // that satisfy dT/dt=0, which requires solving a system of nonlinear equations
  // if the convection coefficient is a function of temperature.

  // -> new
  parameter Modelica.Units.SI.Density
    density = ( PCM.density_solid(PCM.propData.rangeTsolidification[1])
                   + PCM.density_liquid(PCM.propData.rangeTmelting[2]))  / 2.0
    "Average (constant) solid/liquid PCM density";

  // -> new
  parameter Modelica.Units.SI.Length thickness=0.03 "Thickness of layer";

  // -> new
  // BuildingsLib assumes a constant k
  parameter Modelica.Units.SI.ThermalConductivity
    k = ( PCM.conductivity_solid(PCM.propData.rangeTsolidification[1])
        + PCM.conductivity_liquid(PCM.propData.rangeTmelting[2]))  / 2.0
    "Average (constant) thermal conductivity";

  Modelica.Units.SI.Temperature T[nSta](start=if stateAtSurface_a then cat(
        1,
        {T_a_start},
        {(T_a_start + (T_b_start - T_a_start)*UA*sum(RNod[k] for k in 1:i - 1))
          for i in 2:nSta}) else {(T_a_start + (T_b_start - T_a_start)*UA*sum(
        RNod[k] for k in 1:i)) for i in 1:nSta}, each nominal=300)
    "Temperature at the states";

  Modelica.Units.SI.HeatFlowRate Q_flow[nSta + 1](each start=0)
    "Heat flow rates to each state";
  Modelica.Units.SI.SpecificInternalEnergy u[nSta](each start=2.7E5, each
      nominal=2.7E5) "Definition of specific internal energy";

  parameter Boolean stateAtSurface_a=true
    "=true, a state will be at the surface a"
    annotation (Dialog(tab="Dynamics"),
                Evaluate=true);
  parameter Boolean stateAtSurface_b=true
    "=true, a state will be at the surface b"
    annotation (Dialog(tab="Dynamics"),
                Evaluate=true);

  // -> new
  replaceable package PCM =
    slPCMlib.Media_generic.generic_7thOrderSmoothStep
    constrainedby slPCMlib.Interfaces.partialPCM
    annotation (Dialog(group="PCM and phase transition model"),
                choicesAllMatching=true);

  // -> new
  replaceable slPCMlib.Interfaces.phTransMod_MeltingCurve_Algebraic phTrModel[
    nSta](redeclare package PCM = PCM) constrainedby
    slPCMlib.Interfaces.basicPhTransModel(redeclare package PCM = PCM)
    annotation (Dialog(group="PCM and phase transition model"),
      choicesAllMatching=true);

  //replaceable parameter Buildings.HeatTransfer.Data.BaseClasses.Material material
    //"Material from Data.Solids, Data.SolidsPCM or Data.Resistances" annotation (
 //    choicesAllMatching=true, Placement(transformation(extent={{60,60},{80,80}})));

  parameter Boolean steadyStateInitial=false
    "=true initializes dT(0)/dt=0, false initializes T(0) at fixed temperature using T_a_start and T_b_start"
        annotation (Dialog(group="Initialization"), Evaluate=true);
  parameter Modelica.Units.SI.Temperature T_a_start=293.15
    "Initial temperature at port_a, used if steadyStateInitial = false"
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
  parameter Modelica.Units.SI.Temperature T_b_start=293.15
    "Initial temperature at port_b, used if steadyStateInitial = false"
    annotation (Dialog(group="Initialization", enable=not steadyStateInitial));

  // -> new
  parameter Integer nSta2=1
  "Number of states in a material (do not overwrite, used to work around Dymola 2017 bug)"
     annotation (Evaluate=true, HideResult=true, Dialog(enable=false, tab="Advanced"));
protected
  final parameter Integer nSta=
    max(nSta2,
        if stateAtSurface_a or stateAtSurface_b then 2 else 1)
    "Number of state variables";
  final parameter Integer nR=nSta+1 "Number of thermal resistances";
  parameter Modelica.Units.SI.ThermalResistance RNod[nR]=if (stateAtSurface_a
       and stateAtSurface_b) then if (nSta == 2) then {(if i == 1 or i == nR
       then 0 else R/(nSta - 1)) for i in 1:nR} else {(if i == 1 or i == nR
       then 0 elseif i == 2 or i == nR - 1 then R/(2*(nSta - 2)) else R/(nSta
       - 2)) for i in 1:nR} elseif (stateAtSurface_a and (not stateAtSurface_b))
       then {(if i == 1 then 0 elseif i == 2 or i == nR then R/(2*(nSta - 1))
       else R/(nSta - 1)) for i in 1:nR} elseif (stateAtSurface_b and (not
      stateAtSurface_a)) then {(if i == nR then 0 elseif i == 1 or i == nR - 1
       then R/(2*(nSta - 1)) else R/(nSta - 1)) for i in 1:nR} else {R/(if i
       == 1 or i == nR then (2*nSta) else nSta) for i in 1:nR}
    "Thermal resistance";

  parameter Modelica.Units.SI.Mass m[nSta]=(A*thickness*density)*(if (
      stateAtSurface_a and stateAtSurface_b) then if (nSta == 2) then {1/(2*(
      nSta - 1)) for i in 1:nSta} elseif (nSta == 3) then {1/(if i == 1 or i
       == nSta then (2*(nSta - 1)) else (nSta - 1)) for i in 1:nSta} else {1/(
      if i == 1 or i == nSta or i == 2 or i == nSta - 1 then (2*(nSta - 2))
       else (nSta - 2)) for i in 1:nSta} elseif (stateAtSurface_a and (not
      stateAtSurface_b)) then {1/(if i == 1 or i == 2 then (2*(nSta - 1)) else
      (nSta - 1)) for i in 1:nSta} elseif (stateAtSurface_b and (not
      stateAtSurface_a)) then {1/(if i == nSta or i == nSta - 1 then (2*(nSta
       - 1)) else (nSta - 1)) for i in 1:nSta} else {1/(nSta) for i in 1:nSta})
    "Mass associated with the temperature state";

  final parameter Real mInv[nSta]={1/m[i] for i in 1:nSta}
    "Inverse of the mass associated with the temperature state";
  //  if material.steadyState then zeros(nSta) else {1/m[i] for i in 1:nSta}

 // final parameter Modelica.Units.SI.HeatCapacity C[nSta]=m*material.c
 //   "Heat capacity associated with the temperature state";
//  final parameter Real CInv[nSta]=
//    if material.steadyState then zeros(nSta) else {1/C[i] for i in 1:nSta}
//    "Inverse of heat capacity associated with the temperature state";

initial equation
  // The initialization is only done for materials that store energy.
  //  if not material.steadyState then
      if steadyStateInitial then
        //if material.phasechange then
        //  der(u) = zeros(nSta);
       // else
          der(T) = zeros(nSta);
      end if;
  //    else

  // if material.phasechange then
  //   (ud, Td, dT_du) = Buildings.HeatTransfer.Conduction.BaseClasses.der_temperature_u(
  //     c =  material.c,
  //     TSol=material.TSol,
  //     TLiq=material.TLiq,
  //     LHea=material.LHea,
  //     ensureMonotonicity=material.ensureMonotonicity);
  // else
 //    ud    = zeros(Buildings.HeatTransfer.Conduction.nSupPCM);
  //   Td    = zeros(Buildings.HeatTransfer.Conduction.nSupPCM);
 //    dT_du = zeros(Buildings.HeatTransfer.Conduction.nSupPCM);
 //  end if;
equation
    port_a.Q_flow = +Q_flow[1];
    port_b.Q_flow = -Q_flow[end];

    port_a.T-T[1]    = if stateAtSurface_a then 0 else Q_flow[1]*RNod[1];
    T[nSta]-port_b.T = if stateAtSurface_b then 0 else Q_flow[end]*RNod[end];

    for i in 1:nSta-1 loop
       // Q_flow[i+1] is heat flowing from (i) to (i+1)
       // because T[1] has Q_flow[1] and Q_flow[2] acting on it.
       T[i]-T[i+1] = Q_flow[i+1]*RNod[i+1];
    end for;

    // Steady-state heat balance

      // Transient heat conduction
      //if material.phasechange then
        // Phase change material
     //   for i in 1:nSta loop
     //     der(u[i]) = (Q_flow[i]-Q_flow[i+1])*mInv[i];
          // Recalculation of temperature based on specific internal energy
      //    T[i]=Buildings.HeatTransfer.Conduction.BaseClasses.temperature_u(
     //               ud=ud,
     //               Td=Td,
      //              dT_du=dT_du,
      //              u=u[i]);
      //  end for;
     // else
        // Regular material
        for i in 1:nSta loop
          //der(T[i]) = (Q_flow[i]-Q_flow[i+1])*CInv[i];
          phTrModel[i].cp*m[i]*der(T[i]) = Q_flow[i]-Q_flow[i+1];
        end for;
        for i in 1:nSta loop
          u[i]=0; // u is not required in this case
        end for;
    //  end if;

    for i in 1:nSta loop
      phTrModel[i].indVar.T = T[i];
      phTrModel[i].indVar.der_T = der(T[i]);

    end for;

    // Example HeatCapacitorPCM
   // T = port.T;
  //der_T = der(T);

  // input temperature signal to the model
  //phTrModel.indVar.T = T;
  //phTrModel.indVar.der_T = der_T;

 // phTrModel.cp*m*der(port.T) = port.Q_flow;

  annotation ( Icon(coordinateSystem(
          preserveAspectRatio=false,extent={{-100,-100},{100,100}}), graphics={
        Text(
          extent={{-100,-80},{6,-98}},
          textColor={0,0,255},
          textString="%thickness"),
        Text(
          extent={{8,-74},{86,-104}},
          textColor={0,0,255},
          textString="%nSta"),
   Rectangle(
    extent={{-60,80},{60,-80}},     fillColor={215,215,215},
   fillPattern=FillPattern.Solid,    lineColor={175,175,175}),
   Line(points={{-92,0},{90,0}},      color = {0, 0, 0}, thickness = 0.5,
   smooth = Smooth.None),
   Line(points={{8,-40},{-6,-40}},        color = {0, 0, 0}, thickness = 0.5,
   smooth = Smooth.None),
   Line(points={{14,-32},{-12,-32}},      color = {0, 0, 0}, thickness = 0.5,
   smooth = Smooth.None),            Line(
          points={{0,0},{0,-32}},
          color={0,0,0},
          thickness=0.5,
          smooth=Smooth.None),       Rectangle(extent={{-40,6},{-20,-6}},
   lineColor = {0, 0, 0}, lineThickness =  0.5, fillColor = {255, 255, 255},
   fillPattern = FillPattern.Solid), Rectangle(extent={{20,6},{40,-6}},
   lineColor = {0, 0, 0}, lineThickness =  0.5, fillColor = {255, 255, 255},
   fillPattern = FillPattern.Solid),
   Line(points={{66,-40},{52,-40}},       color = {0, 0, 0}, thickness = 0.5,
   smooth = Smooth.None,
   visible=stateAtSurface_b),
   Line(points={{72,-32},{46,-32}},       color = {0, 0, 0}, thickness = 0.5,
   smooth = Smooth.None,
   visible=stateAtSurface_b),            Line(points={{59,0},{59,-32}},
   color = {0, 0, 0}, thickness = 0.5, smooth = Smooth.None,
   visible=stateAtSurface_b),
   Line(points={{-59,0},{-59,-32}},
   color = {0, 0, 0}, thickness = 0.5, smooth = Smooth.None,
   visible=stateAtSurface_a),
   Line(points={{-46,-32},{-72,-32}},     color = {0, 0, 0}, thickness = 0.5,
   smooth = Smooth.None,
   visible=stateAtSurface_a),
   Line(points={{-52,-40},{-66,-40}},     color = {0, 0, 0}, thickness = 0.5,
   smooth = Smooth.None,
   visible=stateAtSurface_a)}),
defaultComponentName="lay",
    Documentation(info="<html>
    <p>
    This is a model of a heat conductor for a single layer of a PCM 
    that computes transient heat conduction.  
    </p> <p>    
    It is a modification of 
        <a href=\"modelica://Buildings.HeatTransfer.Conduction.SingleLayer\">
        Buildings.HeatTransfer.Conduction.SingleLayer</a>
        and can be used with materials and phase transition models 
        contained in <strong>slPCMlib</strong>. 
        </p>
        <p> The 1D heat conduction problem with phase change is formulated adopting the 
        apparent heat capacity method. This method is suitable for real PCM 
        where the phase change occurs over a temperature range of several K. 
        Note that the <a href>SingleLayer</a> 
        in the Buildings Library uses the enthalpy method which can also handle 
        ideal PCM showing a phase change at a single temperature or 
        over a small (<span>&#60;</span> 1K) phase change temperature range. 
        </p><p>
        For details see the documentation in: 
        <a href=\"modelica://Buildings.HeatTransfer.Conduction\">
        Buildings.HeatTransfer.Conduction</a> and
        <a href=\"modelica://slPCMlib\">
        slPCMlib</a>. 
        </p>
</html>",
revisions="<html>
<ul>
<li>2022-06-01; initial version; by Tilman Barz </li>
</ul>
</html>"));
end SingleLayerSlPCMlib;
