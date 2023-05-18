within slPCMlib.Examples;
model heatConduction1D
      "Example considering 1D heat conduction problem"
  extends Modelica.Icons.Example;

  Components.PCMlayer_1D_1port_1symmetry pCMlayer_1D_1port_1symmetry(
    width=0.003,
    redeclare package PCM = slPCMlib.Media_Rubitherm_SP.Rubitherm_SP24E,
    redeclare slPCMlib.Interfaces.phTransModCurveScaleHysteresisDifferentiated
      phTrModel_j,
    densitySLPCM=800,
    initT=297.15,
    n_FD=8)  annotation (Placement(transformation(extent={{30,16},{66,52}})));
  Components.SineT sineT(
    ampl=5,
    timePeriod=240,
    rateT=1/180,
    startT=22 + 273.15)
    annotation (Placement(transformation(extent={{-30,24},{-10,44}})));
equation
  connect(sineT.port, pCMlayer_1D_1port_1symmetry.port) annotation (Line(points={{-9,34},
          {31.8,34}},                            color={191,0,0}));

   //  Advanced.Define.DAEsolver = true;

//annotation (experiment(__Dymola_Algorithm="Radau"));

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=2000,
      __Dymola_NumberOfIntervals=10000,
      Tolerance=1e-09,
      __Dymola_Algorithm="Cvode"),
      __Dymola_experimentFlags="Advanced.Define.DAEsolver = true");
end heatConduction1D;
