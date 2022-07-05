within slPCMlib.Examples;
model heatConduction1D_withVariableInput
  "Example considering 1D heat conduction problem"
  extends Modelica.Icons.Example;


  Modelica.Blocks.Tables.CombiTable1Dv combiTable1D(
    tableOnFile=true,
    tableName="tab1",
    fileName=ModelicaServices.ExternalReferences.loadResource(
        "modelica://slPCMlib/inputT.txt"),
    columns={2},
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint)
    annotation (Placement(transformation(extent={{-62,-46},{-40,-24}})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
    prescribedTemperature
    annotation (Placement(transformation(extent={{-18,-46},{2,-26}})));
  Modelica.Blocks.Sources.ContinuousClock clock
    annotation (Placement(transformation(extent={{-126,-44},{-106,-24}})));
  Components.PCMlayer_1D_1port_1symmetry_2PCM pCMlayer_1D_1port_1symmetry_2PCM(
      redeclare package PCM1 = MediaMyFits.RT35HC_myData_comp1, redeclare
      package PCM2 = MediaMyFits.RT35HC_myData_comp2)
    annotation (Placement(transformation(extent={{30,18},{56,78}})));
equation
  connect(combiTable1D.y[1], prescribedTemperature.T) annotation (Line(points={{
          -38.9,-35},{-38.9,-35.5},{-20,-35.5},{-20,-36}}, color={0,0,127}));
  connect(clock.y, combiTable1D.u[1]) annotation (Line(points={{-105,-34},{-66,
          -34},{-66,-35},{-64.2,-35}}, color={0,0,127}));
  connect(prescribedTemperature.port, pCMlayer_1D_1port_1symmetry_2PCM.port)
    annotation (Line(points={{2,-36},{16,-36},{16,48},{27.4,48}},
                                                                color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=6000,
      __Dymola_NumberOfIntervals=10000,
      Tolerance=1e-07,
      __Dymola_Algorithm="Dassl"));
end heatConduction1D_withVariableInput;
