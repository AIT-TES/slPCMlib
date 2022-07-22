within slPCMlib.Components;
model PCMlayer_1D_1port_1symmetry
  "PCM layer, 1D model, planar geometry, symmetry boundary condition"

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port
  annotation (Placement(transformation(
        origin={0,-100},
        extent={{-10,-10},{10,10}},
        rotation=90), iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-90,0})));

  parameter Modelica.Units.SI.Length width=0.003
    "Width of the PCM layer (direction of heat transfer into the PCM layer)";

  parameter Modelica.Units.SI.Length height=0.10 "Height of the PCM layer"
    annotation (Dialog(tab="General", group="Parameters"));

  parameter Modelica.Units.SI.Length length=0.10 "Length of the PCM layer"
    annotation (Dialog(tab="General", group="Parameters"));

  parameter Modelica.Units.SI.Area htrfArea=height*length "Heat transfer area"
    annotation (Dialog(group="Parameters", enable=false));

  replaceable package PCM =
    slPCMlib.Media_generic.generic_7thOrderSmoothStep
    constrainedby slPCMlib.Interfaces.partialPCM
    annotation (Dialog(group="PCM and phase transition model"),
                choicesAllMatching=true);

  replaceable slPCMlib.Interfaces.phTransModMeltingCurve
    phTrModel_j[n_FD + 1](redeclare package PCM = PCM)
    constrainedby slPCMlib.Interfaces.basicPhTransModel(redeclare package PCM = PCM)
    annotation(Dialog(group="PCM and phase transition model"),
               choicesAllMatching=true);

  parameter Modelica.Units.SI.Density
    densitySLPCM = ( PCM.density_solid(PCM.propData.rangeTsolidification[1])
                   + PCM.density_liquid(PCM.propData.rangeTmelting[2]))  / 2.0
    "Average (constant) solid/liquid PCM density"
    annotation (Dialog(group="PCM and phase transition model"));

  parameter Modelica.Units.SI.Temperature initT(start=273.15 + 20, fixed=true)
    "initial temperatures inside the PCM layer (homogenous T field assumed)"
    annotation (Dialog(group="Initial PCM state"), choicesAllMatching=true);

  Modelica.Units.SI.Temperature T_j[n_FD](start=ones(n_FD)*(initT), fixed=true);

  parameter Integer n_FD(min=1,max=9)=6
    "Number of internal nodes (into the PCM)"
    annotation(Dialog(tab = "General", group = "Discretization"));

  parameter Modelica.Units.SI.Mass mass=width*height*length*densitySLPCM
    "Mass of the PCM element";

  Modelica.Units.SI.MassFraction stateOfCharge(start=0.0)
    "SoC defined by (liquid mass) fraction of stored/absorbed latent heat";

  Modelica.Units.SI.Energy storedEnergy(start=0.0)
    "Total stored/absorbed energy includes sensible and latent heat";



// ---------------------------------------------------------------------------
protected
  parameter Modelica.Units.SI.Length width_i=2*width/(2*n_FD + 1)
    "Thickness of a discrete cell";
//     % Calculation of position of discrete cells
//     z_FD = [z0:dz:zEnd zEnd+dz/2]';

  Real lambda_jp12_BC;
  Real T_jm1[n_FD], T_jp1[n_FD];
  Real lambda_j[n_FD], lambda_jm1[n_FD], lambda_jp1[n_FD];
  Real lambda_jm12[n_FD], lambda_jp12[n_FD];
  Real D_j[n_FD];
  final constant Real eps = Modelica.Constants.small;

// ---------------------------------------------------------------------------
equation

  // input temperature signal to the model
  // (now all properties are automatically updated)
  phTrModel_j[1:n_FD].indVar.T     = T_j[1:n_FD];
  phTrModel_j[1:n_FD].indVar.der_T = der(T_j[1:n_FD]);

  // in addition for BC
  phTrModel_j[n_FD+1].indVar.T =  port.T; // BC
  phTrModel_j[n_FD+1].indVar.der_T =  der(port.T); // BC

  //  der_T_j[:]  = der(T_j[:]); // for plotting only

    // --- compute ODE rhs for each node ---
    for j in 1:n_FD loop

        // --- evaluate properties at node with index j ---
        lambda_j[j] = phTrModel_j[j].lambda;

        // --- use temperatures with index j-1 & j+1 ---
        if (j==1) then
            T_jm1[j]      = port.T; // 1st BC
            T_jp1[j]      = T_j[j+1];
            lambda_jm1[j] = phTrModel_j[n_FD+1].lambda;
            lambda_jp1[j] = phTrModel_j[j+1].lambda;

        elseif (j==n_FD) then
            T_jm1[j]      = T_j[j-1];
            T_jp1[j]      = T_j[j]; // 2nd BC
            lambda_jm1[j] = phTrModel_j[j-1].lambda;
            lambda_jp1[j] = phTrModel_j[j].lambda; // 2nd BC

        else
            T_jm1[j]      = T_j[j-1];
            T_jp1[j]      = T_j[j+1];
            lambda_jm1[j] = phTrModel_j[j-1].lambda;
            lambda_jp1[j] = phTrModel_j[j+1].lambda;

        end if;

        // use harmonic mean
        lambda_jm12[j] = 2 * lambda_j[j] * lambda_jm1[j] / max(lambda_j[j] + lambda_jm1[j],eps);
        lambda_jp12[j] = 2 * lambda_j[j] * lambda_jp1[j] / max(lambda_j[j] + lambda_jp1[j],eps);

        // --- Finite differences for planar geometry ---
        D_j[j] = (lambda_jp12[j] * (T_jp1[j]-T_j[j]) / width_i - lambda_jm12[j] * (T_j[j]-T_jm1[j]) / width_i) / width_i;

        // densitySLPCM =!= constant /-/ instead of modelPhaseTrans_j[j].rho

        // discrete conduction equations
        der(T_j[j]) =  D_j[j] /  (densitySLPCM * phTrModel_j[j].cp);

    end for;

    // --- boundary condition at the port - connect T and Q bound ---
    // thermal conductivity coefficient is calculated using the harmonic mean:
    // NOTE: [n_FD+1] position is used for BC
    lambda_jp12_BC = 2 * phTrModel_j[1].lambda * phTrModel_j[n_FD+1].lambda
                      / max((phTrModel_j[1].lambda + phTrModel_j[n_FD+1].lambda),eps);
    // "Heat flow rate (positive if flowing from outside into the component)"
    port.Q_flow = lambda_jp12_BC * (port.T - T_j[1]) / width_i * htrfArea;


algorithm

    stateOfCharge :=(phTrModel_j[n_FD + 1].xi +phTrModel_j[1].xi)/2.0;
    storedEnergy  :=(phTrModel_j[n_FD + 1].h  +phTrModel_j[1].h)/2.0*mass;
    for j in 1:n_FD-1 loop
      stateOfCharge :=stateOfCharge +(phTrModel_j[j].xi+phTrModel_j[j+1].xi)/2.0;
      storedEnergy  :=storedEnergy  +(phTrModel_j[j].h+phTrModel_j[j+1].h)/2.0*mass;
    end for;
    stateOfCharge :=(stateOfCharge +phTrModel_j[n_FD].xi/2.0)/(n_FD+1);
    storedEnergy  :=(storedEnergy  +phTrModel_j[n_FD].h/2.0*mass)/(n_FD+1);

//     Modelica.Utilities.Streams.writeRealMatrix(
//     fileName="C:/temp/test.mat",
//     matrixName="testMatrix",
//     matrix=[T_j],
//     append=false,
//     format="7");    //"modelica://slPCMlib/out.mat",




      annotation(Icon(
      coordinateSystem(preserveAspectRatio=true),
      graphics={
          Rectangle(
          extent={{80,100},{100,-100}},
          lineColor={0,0,0},
          lineThickness=0.2,
          fillColor={80,180,180},
          fillPattern=fillPattern.Solid),
          Rectangle(
          extent={{-80,100},{80,-100}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid),
          Line(
          points={{100,60},{80,40}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{100,80},{80,60}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{100,40},{80,20}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{100,20},{80,0}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{100,0},{80,-20}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{100,-20},{80,-40}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{100,-40},{80,-60}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{100,-60},{80,-80}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{100,100},{80,80}},
          color={0,0,0},
          thickness=0.2),
        Text(
          extent={{-80,-100},{80,-76}},
          lineColor={0,0,0},
          textStyle={TextStyle.Bold},
          textString="heat conduction"),
        Line(
          points={{-80,-50},{72,-50}},
          color={1,1,1}),
        Polygon(
          points={{22,-50},{14,-60},{14,-40},{22,-50}},
          lineColor={1,1,1},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{80,-50},{72,-54},{72,-46},{80,-50}},
          lineColor={1,1,1},
          lineThickness=0.5),
        Polygon(
          points={{-22,-50},{-36,-64},{-36,-36},{-22,-50}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-54,-50},{-74,-66},{-74,-34},{-54,-50}},
          lineColor={1,1,1},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
          Line(
          points={{100,-80},{80,-100}},
          color={0,0,0},
          thickness=0.2),
        Text(
          extent={{-80,46},{80,-32}},
          lineColor={28,108,200},
          fillColor={225,225,225},
          fillPattern=FillPattern.Solid,
          textStyle={TextStyle.Bold},
          fontName="Adobe Thai",
          textString="PCM"),
        Line(
          points={{-80,60},{72,60}},
          color={1,1,1}),
        Polygon(
          points={{22,60},{14,50},{14,70},{22,60}},
          lineColor={1,1,1},
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{80,60},{72,56},{72,64},{80,60}},
          lineColor={1,1,1},
          lineThickness=0.5),
        Polygon(
          points={{-22,60},{-36,46},{-36,74},{-22,60}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-54,60},{-74,44},{-74,76},{-54,60}},
          lineColor={1,1,1},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(preserveAspectRatio=true)),
      Documentation(info="<html>
    <p>
    This is a 1D heat conduction model in a planar PCM layer 
    (e.g. a wall element). 
    </p>
    <p> 
    The PCM thermal conductivity and apparent (effective) heat capacity are 
    functions of temperature. They are computed by the selected 
    (replaceable) model &lt; phTransMod_XXX_ &gt; 
    for a selected PCM &lt; slPCMlib.Media._XXX_ &gt;.  
    PCM density is assumed <strong>constant</strong>. This means that volume 
    changes due to varying temperature are neglected. 
    The constant density needs to be defined as a parameter 
    &lt; densitySLPCM &gt;. 
    </p>
          <blockquote>          
          <var>lambda = f( T, optional{sgn(dT/dt)} )  <br> 
          cp_app = f( T, optional{sgn(dT/dt)} ) <br>
          rho = constant </var> 
          </blockquote>     
          The discretized 1D modeling domain <code> (z direction --->) </code> is: 
    </p>
    <pre>  
    zStart            zEnd
    ^                 ^      
    |----*----*----*--|--  symmetry condition (Neumann boundary condition) 
         1    2    3       n_FD = 3;  (internal nodes) 
    1    2    3    4       z_FD = 4  
    </pre>
    The model uses a heat port 
    <a href>Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a</a>.
    </p>
    <p>
    Select a &lt;PCM&gt; from <a href>slPCMlib.Media</a>, 
    and a phase transition model 
    &lt;phTrModel&gt; from <a href>slPCMlib.Interfaces</a>. 
    Changes in PCM properties are induced by temperature, 
    use: 
    &lt; phTrModel.indVar.T = T; &gt;. 
    </html>",
  revisions="<html>
       <ul>
       <li>2022-06-01; initial version; by Tilman Barz </li>
       </ul>
       </html>"));
end PCMlayer_1D_1port_1symmetry;
