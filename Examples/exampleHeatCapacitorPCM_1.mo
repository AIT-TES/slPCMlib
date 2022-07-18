within slPCMlib.Examples;
model exampleHeatCapacitorPCM_1 "Example using PCM heat capacitor"
  extends Modelica.Icons.Example;

  Components.SineT sineT(
    ampl=4.5,
    rateT=1/2000,
    startT=273.15 + 17)
    annotation (Placement(transformation(extent={{-22,-20},{-2,0}})));

    // redeclare package PCM = Media.RT64HC,
    // redeclare package PCM = Media.RT60,
  Components.HeatCapacitorSlPCMlib heatCapacitorPCM1(
    T(start=273.15 + 50),
    redeclare package PCM = slPCMlib.Media_Rubitherm_SP.SP29Eu,
    redeclare slPCMlib.Interfaces.phTransModMeltingCurve phTrModel)
    annotation (Placement(transformation(extent={{12,14},{32,34}})));

  Components.HeatCapacitorSlPCMlib heatCapacitorPCM2(
    T(start=273.15 + 50),
    redeclare package PCM = slPCMlib.Media_Rubitherm_SP.SP25E2,
    redeclare slPCMlib.Interfaces.phTransModCurveTrackHysteresis phTrModel)
    annotation (Placement(transformation(extent={{12,-10},{32,10}})));
  Components.HeatCapacitorSlPCMlib heatCapacitorPCM3(
    T(start=273.15 + 50),
    redeclare package PCM = slPCMlib.Media_Rubitherm_SP.SP24E,
    redeclare slPCMlib.Interfaces.phTransModCurveSwitchHysteresisDifferentiated
      phTrModel)
    annotation (Placement(transformation(extent={{12,-34},{32,-14}})));
  Components.HeatCapacitorSlPCMlib heatCapacitorPCM4(
    T(start=273.15 + 50),
    redeclare package PCM = slPCMlib.Media_Rubitherm_SP.SP26E,
    redeclare slPCMlib.Interfaces.phTransModCurveSwitchHysteresisDifferentiated
      phTrModel)
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
          fontSize=12,
          textString="curve track
model"),Text(
          extent={{26,32},{64,14}},
          lineColor={0,0,0},
          fontSize=12,
          textString="melting curve
model"),Text(
          extent={{26,-40},{64,-58}},
          lineColor={0,0,0},
          fontSize=12,
          textString="curve scale
model"),Text(
          extent={{26,-16},{64,-34}},
          lineColor={0,0,0},
          fontSize=12,
          textString="curve switch
model")}),
    experiment(
      StopTime=5000,
      __Dymola_NumberOfIntervals=100000,
      Tolerance=1e-08,
      __Dymola_Algorithm="Dassl"));
end exampleHeatCapacitorPCM_1;
