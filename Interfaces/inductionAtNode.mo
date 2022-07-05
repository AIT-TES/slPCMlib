within slPCMlib.Interfaces;
connector inductionAtNode "Connector for temperature (input signal)"

  input Modelica.Units.SI.Temperature T(start=273.15 + 20) "temperature";
  input Modelica.Units.SI.TemperatureSlope der_T(start=0)
    "Time derivative of temperature (= der(T))";

annotation(Documentation(info="<html>
  It is assumed that phase transitions are 
  <strong>induced by temperature</strong>. 
  Therefore, all properties are computed as function of temperature 
  and temperature (and temperature slope) is the input to the model. 
  </html>"));
end inductionAtNode;
