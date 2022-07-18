within slPCMlib.Examples;
model exampleHeatCapacitorPCM_3 "Example using PCM heat capacitor"
  extends Modelica.Icons.Example;

    // redeclare package PCM = Media.RT64HC,
    // redeclare package PCM = Media.RT60,

  Components.HeatCapacitorSlPCMlib heatCapacitorPCM2(T(start=273.15 + 50),
      redeclare
      slPCMlib.Interfaces.phTransModCurveScaleHysteresisDifferentiated
      phTrModel(redeclare package PCM = Media_Rubitherm_RT.generic_uniform))
    annotation (Placement(transformation(extent={{12,-10},{32,10}})));
  Components.HeatCapacitorSlPCMlib heatCapacitorPCM3(T(start=273.15 + 50),
      redeclare
      slPCMlib.Interfaces.phTransModCurveScaleHysteresisDifferentiated
      phTrModel(redeclare package PCM = Media_Rubitherm_RT.RT64HC))
    annotation (Placement(transformation(extent={{12,-34},{32,-14}})));
  Components.linearT linearT(tValues={0,1,2,3,4,5,6,7,8,9,10,11,12}, TValues={
        331.15,339.15,331.15,337.15,334.15,337.15,334.15,337.15,334.15,337.15,
        334.15,337.15,337.15})
    annotation (Placement(transformation(extent={{-30,-20},{-10,0}})));
equation
  connect(linearT.port, heatCapacitorPCM2.port)
    annotation (Line(points={{-9,-10},{22,-10}}, color={191,0,0}));
  connect(linearT.port, heatCapacitorPCM3.port) annotation (Line(points={{-9,
          -10},{4,-10},{4,-34},{22,-34}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false), graphics={
        Text(
          extent={{28,8},{66,-10}},
          lineColor={0,0,0},
          fontSize=14,
          textString="generic 
uniform"),
        Text(
          extent={{28,-16},{66,-34}},
          lineColor={0,0,0},
          fontSize=14,
          textString="RT64HC")}),
    experiment(
      StopTime=12,
      __Dymola_NumberOfIntervals=10000,
      Tolerance=1e-08,
      __Dymola_Algorithm="Dassl"),
    __Dymola_Commands(file="doSimulationHeatCapacitorPCM1.mos"
        "doSimulationHeatCapacitorPCM1", file=
          "doSimulationHeatCapacitorPCM1.mos" "doSimulationHeatCapacitorPCM1"));
end exampleHeatCapacitorPCM_3;
