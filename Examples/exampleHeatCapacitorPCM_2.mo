within slPCMlib.Examples;
model exampleHeatCapacitorPCM_2 "Example using PCM heat capacitor"
  extends Modelica.Icons.Example;


  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=0.6)
    annotation (Placement(transformation(extent={{24,-46},{44,-26}})));

  Modelica.Blocks.Tables.CombiTable1D combiTable1D(
    tableOnFile=true,
    tableName="tab1",
    fileName=ModelicaServices.ExternalReferences.loadResource("modelica://slPCMlib/inputT.txt"),
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint)
    annotation (Placement(transformation(extent={{-62,-46},{-40,-24}})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-18,-46},{2,-26}})));
  Modelica.Blocks.Sources.Clock clock
    annotation (Placement(transformation(extent={{-94,-44},{-74,-24}})));
  Components.HeatCapacitorPCM heatCapacitorPCM(m=0.01, T(start=302.25, fixed=true))
    annotation (Placement(transformation(extent={{56,24},{76,44}})));
equation
  connect(combiTable1D.y[1], prescribedTemperature.T) annotation (Line(points={{
          -38.9,-35},{-38.9,-35.5},{-20,-35.5},{-20,-36}}, color={0,0,127}));
  connect(prescribedTemperature.port, thermalConductor.port_a)
    annotation (Line(points={{2,-36},{24,-36}}, color={191,0,0}));
  connect(clock.y, combiTable1D.u[1]) annotation (Line(points={{-73,-34},{-66,-34},
          {-66,-35},{-64.2,-35}},      color={0,0,127}));
  connect(thermalConductor.port_b, heatCapacitorPCM.port) annotation (Line(
      points={{44,-36},{56,-36},{56,24},{66,24}},
      color={191,0,0},
      smooth=Smooth.Bezier));
annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
      coordinateSystem(preserveAspectRatio=false)));
end exampleHeatCapacitorPCM_2;
