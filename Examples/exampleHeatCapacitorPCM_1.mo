within slPCMlib.Examples;
model exampleHeatCapacitorPCM_1 "Example using PCM heat capacitor"
  extends Modelica.Icons.Example;

  Components.SineT sineT(ampl=3.5)
    annotation (Placement(transformation(extent={{-66,-20},{-46,0}})));
  Components.HeatCapacitorPCM heatCapacitorPCM(redeclare package PCM =
        Media.generic_7thOrderSmoothStepHysteresis)
    annotation (Placement(transformation(extent={{12,44},{32,64}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=20)
    annotation (Placement(transformation(extent={{-24,-20},{-4,0}})));
equation
  connect(thermalConductor.port_b, heatCapacitorPCM.port) annotation (Line(
        points={{-4,-10},{10,-10},{10,44},{22,44}}, color={191,0,0}));
  connect(sineT.port, thermalConductor.port_a)
    annotation (Line(points={{-45,-10},{-24,-10}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=6000,
      __Dymola_NumberOfIntervals=10000,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end exampleHeatCapacitorPCM_1;
