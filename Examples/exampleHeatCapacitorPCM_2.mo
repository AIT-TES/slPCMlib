within slPCMlib.Examples;
model exampleHeatCapacitorPCM_2 "Example using PCM heat capacitor"
  extends Modelica.Icons.Example;

  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=5)
    annotation (Placement(transformation(extent={{18,-46},{40,-24}})));

  Modelica.Blocks.Tables.CombiTable1Dv combiTable1D(
    tableOnFile=true,
    tableName="tab1",
    fileName=ModelicaServices.ExternalReferences.loadResource(
        "modelica://slPCMlib/Resources/inputT.tabOnFile"),
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint)
    annotation (Placement(transformation(extent={{-50,-46},{-28,-24}})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-16,-46},{6,-24}})));
  Modelica.Blocks.Sources.ContinuousClock clock
    annotation (Placement(transformation(extent={{-80,-44},{-62,-26}})));
  Components.HeatCapacitorSlPCMlib heatCapacitorPCM(
    T(fixed=true),
    redeclare package PCM = slPCMlib.Media_PLUSS_HS.PLUSS_savE_HS36,
    redeclare slPCMlib.Interfaces.phTransModCurveScaleHysteresisDifferentiated
      phTrModel) annotation (Placement(transformation(extent={{60,8},{80,28}})));
equation
  connect(combiTable1D.y[1], prescribedTemperature.T) annotation (Line(points={{-26.9,
          -35},{-18.2,-35}},                               color={0,0,127}));
  connect(prescribedTemperature.port, thermalConductor.port_a)
    annotation (Line(points={{6,-35},{18,-35}}, color={191,0,0}));
  connect(clock.y, combiTable1D.u[1]) annotation (Line(points={{-61.1,-35},{
          -52.2,-35}},                 color={0,0,127}));
  connect(thermalConductor.port_b, heatCapacitorPCM.port) annotation (Line(
        points={{40,-35},{56,-35},{56,8},{70,8}}, color={191,0,0}));
annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
      coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=2000,
      Tolerance=1e-05,
      __Dymola_Algorithm="Dassl"));
end exampleHeatCapacitorPCM_2;
