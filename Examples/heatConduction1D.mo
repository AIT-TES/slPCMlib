within slPCMlib.Examples;
model heatConduction1D
      "Example considering 1D heat conduction problem"
  extends Modelica.Icons.Example;

  Components.SineT sineT(
    ampl=5,
    timePeriod=240,
    rateT=1/180,
    startT=22 + 273.15)
    annotation (Placement(transformation(extent={{-30,24},{-10,44}})));
  Components.PCMlayer_1D_1port_1symmetry pCMlayer_1D_1port_1symmetry
    annotation (Placement(transformation(extent={{14,56},{34,76}})));
  Components.PCMlayer_1D_2ports pCMlayer_1D_2ports
    annotation (Placement(transformation(extent={{16,24},{36,44}})));
equation

   //  Advanced.Define.DAEsolver = true;

//annotation (experiment(__Dymola_Algorithm="Radau"));

  connect(sineT.port, pCMlayer_1D_1port_1symmetry.port)
    annotation (Line(points={{-9,34},{6,34},{6,66},{15,66}}, color={191,0,0}));
  connect(sineT.port, pCMlayer_1D_2ports.portA)
    annotation (Line(points={{-9,34},{17,34}}, color={191,0,0}));
  connect(sineT.port, pCMlayer_1D_2ports.portB) annotation (Line(points={{-9,34},
          {6,34},{6,50},{38,50},{38,34},{35,34}}, color={191,0,0}));
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
