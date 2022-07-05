within slPCMlib.BasicUtilities;
function cubicHermiteSplineEval
  "Evaluate piecewise cubic hermite splines"
  extends Modelica.Icons.Function;
  input Real T;
  input Integer len_x;
  input Real data_x[len_x];
  input Real data_y[len_x];
  input Real m_k[len_x];
  input Real iy_start[len_x];
  input Real iy_scaler;
  output Real Xi, dXi;

protected
  Integer intNum "interval number";
  Real delta, t;
  Real ih00, ih01, ih10, ih11;
  Real h00,   h01,  h10,  h11;

algorithm

//   assert(  (breaks[1] <= T and T <= breaks[pieces+1]),
//        ("Medium model outside feasible range! Problem with T/°C = " + String(T)),
//        AssertionLevel.error);                          // if not true then

//  Modelica.Utilities.Streams.print(" ---> enter splineEval => " + " T = " + String(T));


  // find the interval
  intNum := 0;
  if T < data_x[1] then
    Xi  :=0.0;
    dXi :=0.0;
  elseif T > data_x[len_x] then
    Xi  :=1.0;
    dXi :=0.0;
  else
//     Xi  :=1.0;
//     dXi :=0.0;
    while T > data_x[intNum+1] loop
      intNum := intNum + 1;
  assert(  intNum < len_x,
       ("intNum outside range! Problem with T/°C = " + String(T)),
       AssertionLevel.error);                          // if not true then
//    Modelica.Utilities.Streams.print(" - loop = " + String[intNum] + " | T = " + String(T)
//                                        + " | breaks = " + String(breaks[intNum]));
    end while;
//   Modelica.Utilities.Streams.print(" - int found: breaks[" + String[intNum] + "]=" + String(breaks[intNum]) + " <= " + String(T));

    // eval spline in element
    delta := data_x[intNum + 1] - data_x[intNum];
    t     := (T - data_x[intNum])/delta;

    ih00 := 2.0/4.0*t*t*t*t - 3.0/3.0*t*t*t + 1.0*t;
    ih10 := 1.0/4.0*t*t*t*t - 2.0/3.0*t*t*t + 1.0/2.0*t*t;
    ih01 := -2.0/4.0*t*t*t*t + 3.0/3.0*t*t*t;
    ih11 := 1.0/4.0*t*t*t*t - 1.0/3.0*t*t*t;

    h00 := 2.0*t*t*t - 3.0*t*t + 1.0;
    h10 := t*t*t - 2.0*t*t + t;
    h01 := -2.0*t*t*t + 3.0*t*t;
    h11 := t*t*t - t*t;


    Xi  :=(  ih00*data_y[intNum]
           + ih10*delta*m_k[intNum]
           + ih01*data_y[intNum+1]
           + ih11*delta*m_k[intNum+1])
              * delta * iy_scaler + iy_start[intNum];

    dXi :=(  h00*data_y[intNum]
          + h10*delta*m_k[intNum]
          + h01*data_y[intNum+1]
          + h11*delta*m_k[intNum+1])  * iy_scaler;

  end if;

  annotation (Icon(graphics,
      coordinateSystem(preserveAspectRatio=false)), Diagram(graphics,
      coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
          <p>
          Function for the evaluation of piecewise cubic hermite splines 
          including their integrals and derivatives. 
          See e.g.  
          </p>
          <blockquote>          
          <p>
          <a href>https://en.wikipedia.org/wiki/Cubic_Hermite_spline</a>. 
          </p>          
          </blockquote>          
          <p> 
          The result (piecewise polynom) is one time continuously 
          differentiable (C1 smooth). 
          It is used for modeling the phase transition function xi .           
          </p>
          </html>",
    revisions="<html>
        <ul>
        <li>2022-06-01; initial version; by Tilman Barz </li>
        </ul>
        </html>"));
end cubicHermiteSplineEval;
