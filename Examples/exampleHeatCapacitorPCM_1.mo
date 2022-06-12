within slPCMlib.Examples;
model exampleHeatCapacitorPCM_1 "Example using PCM heat capacitor"
  extends Modelica.Icons.Example;

  Components.SineT sineT(
    ampl=2.5,
    rateT=1/2000,
    startT=273.15 + 26)
    annotation (Placement(transformation(extent={{-22,-20},{-2,0}})));



    // redeclare package PCM = Media.RT64HC,
    // redeclare package PCM = Media.RT60,
  Components.HeatCapacitorPCM heatCapacitorPCM1(
    T(start=273.15 + 50),
    redeclare slPCMlib.Interfaces.phTransModMeltingCurve phTrModel(redeclare
        package PCM = Media.generic_GumbelMinimum))
    annotation (Placement(transformation(extent={{12,14},{32,34}})));

  Components.HeatCapacitorPCM heatCapacitorPCM2(
    T(start=273.15 + 50),
    redeclare slPCMlib.Interfaces.phTransModCurveTrackHysteresis phTrModel(redeclare
        package  PCM = Media.generic_GumbelMinimum))
    annotation (Placement(transformation(extent={{12,-10},{32,10}})));
  Components.HeatCapacitorPCM heatCapacitorPCM3(
    T(start=273.15 + 50),
    redeclare slPCMlib.Interfaces.phTransModCurveScaleHysteresisDifferentiated phTrModel(redeclare
        package PCM =  Media.generic_GumbelMinimum))
    annotation (Placement(transformation(extent={{12,-34},{32,-14}})));
  Components.HeatCapacitorPCM heatCapacitorPCM4(
    T(start=273.15 + 50),
    redeclare slPCMlib.Interfaces.phTransModCurveSwitchHysteresisDifferentiated phTrModel(redeclare
        package  PCM =
          Media.generic_GumbelMinimum))
    annotation (Placement(transformation(extent={{12,-58},{32,-38}})));
equation
  connect(sineT.port, heatCapacitorPCM1.port) annotation (Line(points={{-1,-10},
          {6,-10},{6,14},{22,14}}, color={191,0,0}));
  connect(sineT.port, heatCapacitorPCM2.port)
    annotation (Line(points={{-1,-10},{22,-10}}, color={191,0,0}));
  connect(sineT.port, heatCapacitorPCM3.port) annotation (Line(points={{-1,-10},
          {6,-10},{6,-34},{22,-34}}, color={191,0,0}));
  connect(sineT.port, heatCapacitorPCM4.port) annotation (Line(points={{-1,-10},
          {6,-10},{6,-58},{22,-58}}, color={191,0,0}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false), graphics={
        Text(
          extent={{26,8},{64,-10}},
          lineColor={0,0,0},
          fontSize=14,
          textString="curve track
model"),Text(
          extent={{26,32},{64,14}},
          lineColor={0,0,0},
          fontSize=14,
          textString="melting curve
model"),Text(
          extent={{26,-16},{64,-34}},
          lineColor={0,0,0},
          fontSize=14,
          textString="curve scale
model"),Text(
          extent={{26,-40},{64,-58}},
          lineColor={0,0,0},
          fontSize=14,
          textString="curve switch
model")}),
    experiment(
      StopTime=5000,
      __Dymola_NumberOfIntervals=100000,
      Tolerance=1e-08,
      __Dymola_Algorithm="Dassl"),
    __Dymola_Commands(file="doSimulationHeatCapacitorPCM1.mos"
        "doSimulationHeatCapacitorPCM1", file=
          "doSimulationHeatCapacitorPCM1.mos" "doSimulationHeatCapacitorPCM1"));
end exampleHeatCapacitorPCM_1;
