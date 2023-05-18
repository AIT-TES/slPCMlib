within slPCMlib.Examples;
model exampleHeatCapacitorPCM_3 "Example using PCM heat capacitor"
  extends Modelica.Icons.Example;



  Components.linearT linearT(tValues={0,1,2,3,4,5,6,7,8,9,10,11,12}, TValues={331.15,
        339.15,331.15,337.15,334.15,337.15,334.15,337.15,334.15,337.15,334.15,337.15,
        337.15})
    annotation (Placement(transformation(extent={{-30,-20},{-10,0}})));
  Components.HeatCapacitorSlPCMlib heatCapacitorSlPCMlib(redeclare package PCM =
        slPCMlib.Media_generic.generic_uniform, redeclare
      slPCMlib.Interfaces.phTransModCurveScaleHysteresisDifferentiated
      phTrModel)
    annotation (Placement(transformation(extent={{12,-10},{32,10}})));
  Components.HeatCapacitorSlPCMlib heatCapacitorSlPCMlib1(redeclare package PCM
      = slPCMlib.Media_Rubitherm_RT.Rubitherm_RT64HC,
                                            redeclare
      slPCMlib.Interfaces.phTransModCurveScaleHysteresisDifferentiated
      phTrModel)
    annotation (Placement(transformation(extent={{12,-34},{32,-14}})));
equation
  connect(linearT.port, heatCapacitorSlPCMlib.port)
    annotation (Line(points={{-9,-10},{22,-10}}, color={191,0,0}));
  connect(linearT.port, heatCapacitorSlPCMlib1.port) annotation (Line(points={{
          -9,-10},{4,-10},{4,-34},{22,-34}}, color={191,0,0}));
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
      Tolerance=1e-05,
      __Dymola_Algorithm="Dassl"));
end exampleHeatCapacitorPCM_3;
