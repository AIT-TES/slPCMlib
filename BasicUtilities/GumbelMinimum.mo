within slPCMlib.BasicUtilities;
package GumbelMinimum
  "Library of Gumbel Minimum distribution functions"
   extends Modelica.Icons.Package;

  function density "Density of Gumbel Minimum distribution"
    extends Modelica.Math.Distributions.Interfaces.partialDensity;
      input  Real mu = 0     "Maximum, location parameter" annotation(Dialog);
      input  Real beta = 1   "Shape parameter, smaller is sharper, beta>0" annotation(Dialog);
  protected
      Real tmp1, tmp2, tmp3;
  algorithm
      tmp1 :=(u - mu)/beta;
      tmp2 :=Modelica.Math.exp(tmp1);
      tmp3 :=Modelica.Math.exp(-tmp2);
      y :=1.0/beta*tmp2 .* tmp3;
  end density;

  function cumulative
    "Cumulative distribution function of Gumbel Minimum distribution"
    extends Modelica.Math.Distributions.Interfaces.partialCumulative;
      input  Real mu =   0   "Maximum, location parameter" annotation(Dialog);
      input  Real beta = 1   "Shape parameter, smaller is sharper, beta>0" annotation(Dialog);
  protected
      Real tmp1, tmp2, tmp3;
  algorithm
      tmp1 :=(u - mu)/beta;
      tmp2 :=Modelica.Math.exp(tmp1);
      tmp3 :=Modelica.Math.exp(-tmp2);
      y    :=(1.0 - tmp3);
  end cumulative;

  function quantile "Quantile of Gumbel Minimum distribution"
    extends Modelica.Math.Distributions.Interfaces.partialQuantile;
      input  Real mu =   0   "Maximum, location parameter" annotation(Dialog);
      input  Real beta = 1   "Shape parameter, smaller is sharper, beta>0" annotation(Dialog);
  algorithm
      y :=Modelica.Math.log(-Modelica.Math.log(1.0 - u))*beta + mu;
  end quantile;
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
