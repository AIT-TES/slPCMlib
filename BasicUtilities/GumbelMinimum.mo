within slPCMlib.BasicUtilities;
package GumbelMinimum

  function PDF "Gumbel Minimum PDF"
    extends Modelica.Icons.Function;
      input  Real x     "Random number";
      input  Real mu    "Maximum, location parameter";
      input  Real beta  "Shape parameter, smaller is sharper, beta>0";
      output Real phi   "Density";
  protected
      Real tmp1, tmp2, tmp3;
  algorithm
      tmp1 :=(x - mu)/beta;
      tmp2 :=Modelica.Math.exp(tmp1);
      tmp3 :=Modelica.Math.exp(-tmp2);
      phi :=1.0/beta*tmp2 .* tmp3;
  end PDF;

  function CDF "Gumbel Minimum CDF"
    extends Modelica.Icons.Function;
      input  Real x       "Random number";
      input  Real mu      "Maximum, location parameter";
      input  Real beta    "Shape parameter, smaller is sharper, beta>0";
      output Real P       "Probability";
  protected
      Real tmp1, tmp2, tmp3;
  algorithm
      tmp1 :=(x - mu)/beta;
      tmp2 :=Modelica.Math.exp(tmp1);
      tmp3 :=Modelica.Math.exp(-tmp2);
      P    :=(1.0 - tmp3);
  end CDF;

  function invCDF
    "Gumbel Minimum inverse CDF"
    extends Modelica.Icons.Function;
      input  Real P       "Probability";
      input  Real mu      "Maximum, location parameter";
      input  Real beta    "Shape parameter, smaller is sharper, beta>0";
      output Real x       "Random number";
  algorithm
      x :=Modelica.Math.log(-Modelica.Math.log(1.0 - P))*beta + mu;
  end invCDF;
annotation (Documentation(info="<html>
      <p>
      This package provides the PDF, CDF and inverse CDF of the 
      Gumbel extreme value type I (Minimum) distribution, see
      <blockquote>          
      NIST/SEMATECH. e-Handbook of Statistical Methods. 
      <a href>http://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm</a>
      </blockquote> 
      The functions are asymmetric (steep on the right decline). The PDF reads <br>
      </p>          
      <pre>
      <var>y_GumbMin = 1/beta * exp( (x-mue)/beta) .* exp(-exp((x-mue)/beta))    </var> 
      <var>mue</var>  : maximum, location parameter
      <var>beta</var>  : shape parameter, smaller is sharper, <var>beta>0</var> 
      </pre>
      <p>
      The functions are asymmetric (PDF is steep on the right decline). 
      </html>",
    revisions="<html>
        <ul>
        <li>2022-06-01; initial version; by Tilman Barz </li>
        </ul>
        </html>"));
end GumbelMinimum;
