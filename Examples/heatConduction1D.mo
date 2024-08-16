within slPCMlib.Examples;
model heatConduction1D
      "Example considering 1D heat conduction problem"
  extends Modelica.Icons.Example;

  Components.PCMlayer_1D_1port_1symmetry pCMlayer_1D_1port_1symmetry(n_FD=1)
    annotation (Placement(transformation(extent={{48,24},{68,44}})));
  Components.SineT sineT(
    ampl=5,
    timePeriod=240,
    rateT=1/180,
    startT=22 + 273.15)
    annotation (Placement(transformation(extent={{-30,24},{-10,44}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(Q_flow=2000)
    annotation (Placement(transformation(extent={{-28,-10},{-8,10}})));
  Components.PCMlayer_1D_1port_1symmetry pCMlayer_1D_1port_1symmetry1(
      redeclare slPCMlib.Interfaces.phTransMod_MeltingCurve_Differentiated
      phTrModel_profile, n_FD=1)
    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
  Components.PCMlayer_1D_2ports pCMlayer_1D_2ports(n_FD=1)
    annotation (Placement(transformation(extent={{16,24},{36,44}})));
equation

   //  Advanced.Define.DAEsolver = true;

//annotation (experiment(__Dymola_Algorithm="Radau"));

  connect(fixedHeatFlow.port, pCMlayer_1D_1port_1symmetry1.port)
    annotation (Line(points={{-8,0},{41,0}},      color={191,0,0}));
  connect(sineT.port, pCMlayer_1D_2ports.portA)
    annotation (Line(points={{-9,34},{17,34}}, color={191,0,0}));
  connect(pCMlayer_1D_2ports.portB, pCMlayer_1D_1port_1symmetry.port)
    annotation (Line(points={{35,34},{49,34}}, color={191,0,0}));
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
