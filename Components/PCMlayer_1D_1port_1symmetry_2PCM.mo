within slPCMlib.Components;
model PCMlayer_1D_1port_1symmetry_2PCM
  "PCM layer (planar geometry) 1-dimensional model, 1 port & 1 symmetry boundary condition"

extends Modelica.Icons.ObsoleteModel;
  //
  //  heat conduction equation, second order PDE with d^2x/dz^2
  //
  //  lambda = f(T); rho = constant; cp_app = f(T, optional{direction(T)})
  //
  //      zStart            zEnd
  //      ^                 ^
  //      |----*----*----*--|--  symmetry condition = Neumann boundary condition
  //           1    2    3       n_FD = 3;  (internal points)
  //      1    2    3    4       z_FD = 4,
  //

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port
  annotation (Placement(transformation(
        origin={0,-100},
        extent={{-10,-10},{10,10}},
        rotation=90), iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-70,0})));

  parameter Modelica.Units.SI.Length width=0.003
    "Width of the PCM layer (direction of heat transfer into the PCM layer)"
    annotation (Dialog(tab="General", group="Parameters"));

  parameter Modelica.Units.SI.Length height=0.10 "Height of the PCM layer"
    annotation (Dialog(tab="General", group="Parameters"));

  parameter Modelica.Units.SI.Length length=0.10 "Length of the PCM layer"
    annotation (Dialog(tab="General", group="Parameters"));

  parameter Modelica.Units.SI.Area htrfArea=height*length "Heat transfer area";

  parameter Integer n_FD(min=1)=6
      "Number of internal nodes (into the PCM)"
      annotation(Dialog(tab = "General", group = "Discretization"));

  // --- PCM and phase transition model ---
  replaceable package PCM1 =
      slPCMlib.Media_generic.generic_7thOrderSmoothStep
            constrainedby slPCMlib.Interfaces.partialPCM
            annotation (Dialog(tab = "General", group="PCM"),
                        choicesAllMatching=true);

  replaceable package PCM2 =
      slPCMlib.Media_generic.generic_GumbelMinimum
            constrainedby slPCMlib.Interfaces.partialPCM
            annotation (Dialog(tab = "General", group="PCM"),
                        choicesAllMatching=true);

  // (n_FD+1 is for evaluating the properties at the left boundary/ port)
  PCM1.ModelPhaseTrans model1PhaseTrans_j[n_FD+1]
       "vector of phase transition models for each discrete node"
                    annotation (Dialog(tab = "General", group="PCM"),
                          choicesAllMatching=true);
  PCM2.ModelPhaseTrans model2PhaseTrans_j[n_FD+1]
       "vector of phase transition models for each discrete node"
                    annotation (Dialog(tab = "General", group="PCM"),
                          choicesAllMatching=true);



  parameter Modelica.Units.SI.Density densitySLPCM=(PCM1.rho_liquid + PCM1.rho_solid)
      /2.0 "Average (constant) solid/liquid PCM density"
    annotation (Dialog(tab="General", group="PCM"));

  // --- port and internal temperatures ---
  Modelica.Units.SI.Temperature T_j[n_FD](start=ones(n_FD)*(PCM1.rangeTsolidification[
        1] - (PCM1.rangeTmelting[2] - PCM1.rangeTsolidification[1])/2),
      displayUnit="degC") "PCM temperatures inside the layer";
                                                       //*(47.5+273.15),
//             *(   PCM1.rangeTsolidification[1]
//                - (PCM1.rangeTmelting[2] - PCM1.rangeTsolidification[1])/2),
//      annotation(Dialog(tab = "General", group = "PCM"));

//   Modelica.SIunits.TemperatureSlope der_T_j[n_FD](start=zeros(n_FD))
//     "Time derivative of temperature (= der(T))";

//     % Calculation of position of discrete cells
//     z_FD = [z0:dz:zEnd zEnd+dz/2]';
    // % Thickness of a discrete cell in radial direction
    // dz = 2*(zEnd-z0) / (2*n_FD +1);

// ---------------------------------------------------------------------------

protected
  parameter Modelica.Units.SI.Length width_i=2*width/(2*n_FD + 1)
    "Thickness of a discrete cell";

  Real lambda_jp12_BC;
  Real T_jm1[n_FD], T_jp1[n_FD];
  Real lambda_j[n_FD], lambda_jm1[n_FD], lambda_jp1[n_FD];
  Real lambda_jm12[n_FD], lambda_jp12[n_FD];
  Real D_j[n_FD];

// ---------------------------------------------------------------------------
equation

    // --- connect temperatures to nodes ---
    // (now all properties are automatically updated)
//    modelPhaseTrans_j[1:n_FD].indVar.T =  Modelica.SIunits.Conversions.to_degC(T_j[1:n_FD]);
    model1PhaseTrans_j[1:n_FD].indVar.T =  (T_j[1:n_FD]);
    // in addition for BC
    model1PhaseTrans_j[n_FD+1].indVar.T =  (port.T); // BC

    model2PhaseTrans_j[1:n_FD].indVar.T =  (T_j[1:n_FD]);
    // in addition for BC
    model2PhaseTrans_j[n_FD+1].indVar.T =  (port.T); // BC

  //  der_T_j[:]  = der(T_j[:]); // for plotting only

    // --- compute ODE rhs for each node ---
    for j in 1:n_FD loop

        // --- evaluate properties at node with index j ---
        lambda_j[j] = model1PhaseTrans_j[j].lambda;

        // --- use temperatures with index j-1 & j+1 ---
        if (j==1) then
            T_jm1[j]      = port.T; // 1st BC
            T_jp1[j]      = T_j[j+1];
            lambda_jm1[j] = model1PhaseTrans_j[n_FD+1].lambda;
            lambda_jp1[j] = model1PhaseTrans_j[j+1].lambda;

        elseif (j==n_FD) then
            T_jm1[j]      = T_j[j-1];
            T_jp1[j]      = T_j[j]; // 2nd BC
            lambda_jm1[j] = model1PhaseTrans_j[j-1].lambda;
            lambda_jp1[j] = model1PhaseTrans_j[j].lambda; // 2nd BC

        else
            T_jm1[j]      = T_j[j-1];
            T_jp1[j]      = T_j[j+1];
            lambda_jm1[j] = model1PhaseTrans_j[j-1].lambda;
            lambda_jp1[j] = model1PhaseTrans_j[j+1].lambda;

        end if;

        // densitySLPCM =!= constant /-/ instead of modelPhaseTrans_j[j].rho

        lambda_jm12[j] = 2 * lambda_j[j] * lambda_jm1[j] / (lambda_j[j] + lambda_jm1[j]);
        lambda_jp12[j] = 2 * lambda_j[j] * lambda_jp1[j] / (lambda_j[j] + lambda_jp1[j]);

        // --- Finite differences for planar geometry ---
        D_j[j] = (lambda_jp12[j] * (T_jp1[j]-T_j[j]) / width_i - lambda_jm12[j] * (T_j[j]-T_jm1[j]) / width_i) / width_i;

        // discrete conduction equations
        der(T_j[j]) =  D_j[j] /  (densitySLPCM * (model1PhaseTrans_j[j].cp + model2PhaseTrans_j[j].cp));

    end for;

    // --- boundary condition at the port - connect T and Q bound ---
    // thermal conductivity coefficient is calculated using the harmonic mean:
    lambda_jp12_BC = 2 * model1PhaseTrans_j[1].lambda * model1PhaseTrans_j[n_FD+1].lambda
                      / (model1PhaseTrans_j[1].lambda + model1PhaseTrans_j[n_FD+1].lambda);
    // "Heat flow rate (positive if flowing from outside into the component)"
    port.Q_flow = lambda_jp12_BC * (port.T - T_j[1]) / width_i * htrfArea;

//     Modelica.Utilities.Streams.writeRealMatrix(
//     fileName="C:/temp/test.mat",
//     matrixName="testMatrix",
//     matrix=[T_j],
//     append=false,
//     format="7");                 //"modelica://slPCMlib/out.mat",
                           //

      annotation(Icon(coordinateSystem(extent={{-60,-120},{40,120}}), graphics={
          Rectangle(
          extent={{20,120},{40,-120}},
          lineColor={0,0,0},
          lineThickness=0.2,
          fillColor={80,180,180},
          fillPattern=fillPattern.Solid),
          Rectangle(
          extent={{-60,120},{20,-120}},
          lineColor={28,108,200},
          radius=0,
          fillColor={200,200,200},
          fillPattern=FillPattern.Solid),
          Text(
          extent={{-67,56},{67,-56}},
          lineColor={162,29,33},
          fontName="Adobe Thai",
          textStyle={TextStyle.Bold},
          origin={-10,21},
          rotation=-90,
          textString="2 PCM"),
          Line(
          points={{40,60},{20,40}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,80},{20,60}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,40},{20,20}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,20},{20,0}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,0},{20,-20}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,-20},{20,-40}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,-40},{20,-60}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,-60},{20,-80}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,-100},{20,-120}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,-80},{20,-100}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,100},{20,80}},
          color={0,0,0},
          thickness=0.2),
          Line(
          points={{40,120},{20,100}},
          color={0,0,0},
          thickness=0.2),
        Text(
          extent={{-60,-112},{6,-88}},
          lineColor={0,0,0},
          textStyle={TextStyle.Bold},
          textString="1D heat 
conduction"),
        Line(
            points={{-60,-100},{20,-100},{2,-100}},
            color={1,1,1},
            arrow=yes,
            thickness=0.3),
        Line(
            points={{10,-94},{20,-100}},
            color={1,1,1},
            thickness=0.3),
        Line(
            points={{10,-106},{20,-100}},
            color={1,1,1},
            thickness=0.3)}),   Diagram(coordinateSystem(extent={{-60,-120},{40,
            120}})));
end PCMlayer_1D_1port_1symmetry_2PCM;
