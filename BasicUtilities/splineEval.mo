within slPCMlib.BasicUtilities;
function splineEval "Evaluate spline"
  extends Modelica.Icons.Function;
  input Real T;
  input Integer pieces; // num intvals
  input Integer order[pieces];
  input Real breaks[pieces+1];
  input Real coefs[:,:];
  output Real Xi, dXi;
  // Integer n=scalar(size(breaks)) - 1 "Max value";
protected
  Integer intNum "interval number";
  Integer nO;
  Real dT;
algorithm

  assert(  (breaks[1] <= T and T <= breaks[pieces+1]),
       ("Medium model outside feasible range! Problem with T/°C = " + String(T)),
       AssertionLevel.error);                          // if not true then
  //      Modelica.Utilities.Streams.print(" ---> enter splineEval => " + " T = " + String(T));

  // find the interval
  intNum := 0;
  while T >= breaks[intNum+1] loop
    intNum := intNum + 1;
    //           Modelica.Utilities.Streams.print(" - loop = " + String(intNum) + " | T = " + String(T)
    //                                             + " | breaks = " + String(breaks[intNum]));
  end while;
  //          Modelica.Utilities.Streams.print(" - int found: breaks[" + String(intNum) + "]=" + String(breaks[intNum]) + " <= " + String(T));

  nO :=order[intNum]; // 0 order is integral
  dT :=T - breaks[intNum];
  //           Modelica.Utilities.Streams.print(" - nO = " + String(nO));
  //           Modelica.Utilities.Streams.print(" - dT = " + String(dT));

//      iXi :=  coefs[intNum,1]
//           +  coefs[intNum,2]*dT
//           +  coefs[intNum,3]*dT*dT
//           +  coefs[intNum,4]*dT*dT*dT
//           +  coefs[intNum,5]*dT*dT*dT*dT
//           +  coefs[intNum,6]*dT*dT*dT*dT*dT
//           +  coefs[intNum,7]*dT*dT*dT*dT*dT*dT;
                                // value is integral
      Xi :=      coefs[intNum,2]
          +  2.0*coefs[intNum,3]*dT
          +  3.0*coefs[intNum,4]*dT*dT
          +  4.0*coefs[intNum,5]*dT*dT*dT
          +  5.0*coefs[intNum,6]*dT*dT*dT*dT
          +  6.0*coefs[intNum,7]*dT*dT*dT*dT*dT;
                                 // 1st deriv is value
     dXi :=  2.0*coefs[intNum,3]
          +  6.0*coefs[intNum,4]*dT
          + 12.0*coefs[intNum,5]*dT*dT
          + 20.0*coefs[intNum,6]*dT*dT*dT
          + 30.0*coefs[intNum,7]*dT*dT*dT*dT;
                                // 2nd deriv is 1st deriv

  // Modelica.Utilities.Streams.print("intNum = " + String(intNum) + " | nO = " + String(nO)
  //                                   + " | breaks = " + String(breaks[intNum]));


  annotation (Icon(graphics,
      coordinateSystem(preserveAspectRatio=false)), Diagram(graphics,
      coordinateSystem(preserveAspectRatio=false)),
          Documentation(info="<html>
          <p>
          Function for the evaluation of splines described in  
          </p>
          <blockquote>          
          <p>
          Barz, T., Krämer, J., & Emhofer, J. (2020). Identification of Phase 
          Fraction–Temperature Curves from Heat Capacity Data for Numerical 
          Modeling of Heat Transfer in Commercial Paraffin Waxes. 
          Energies, 13(19), 5149.
          <a href>doi.org/10.3390/en13195149</a>. 
          </p>          
          </blockquote>          
          <p> 
          </p>
          </html>",
    revisions="<html>
        <ul>
        <li>2022-06-01; initial version; by Tilman Barz </li>
        </ul>
        </html>"));
end splineEval;
