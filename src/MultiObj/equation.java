package MultiObj;


import IMRT_DAO.TreatmentPlan;

import java.util.ArrayList;
/**
 * This class is thinked by the 2d version of the IMRT problem.
 * To extend to other dimensions is necessary add other surface function.
 */
public class equation {
	double Pa,Pb,Pc;//Coeficent of the parabola formula
	double m,b;//number linear equation
	double x1, y1, x2, y2, x3, y3;//Variable that represent point a b and c where the under number 1 represente the x value and 2 the y value
	
	

	public equation(TreatmentPlan up,TreatmentPlan down) {
		x1=up.scores[1];
		y1=up.scores[2];
		x2=down.scores[1];
		y2=down.scores[2];	
		
		getLinearEquation();
	}
	
	public equation(TreatmentPlan up,TreatmentPlan down,TreatmentPlan middle) {
		x1=up.scores[0];
		y1=up.scores[1];
		x2=down.scores[0];
		y2=down.scores[1];
		x3=middle.scores[0];
		y3=middle.scores[1];
		
		
		
		parabolaFunctionGeneration();
	}
	
	/**
	 * This function calulate the coeficient of the prabola function
	 * y = c + b*x + a*x^2
	 * @param x1 x value of the point 1
	 * @param y1 y value of the point 1
	 * @param x2 x value of the point 2
	 * @param y2 y value of the point 2
	 * @param x3 x value of the point 3
	 * @param y4 y value of the point c3
	 * @return
	 */
	private void parabolaFunctionGeneration() {
		
		double denom = (x1 - x2)*(x1 - x3)*(x2 - x3);
		Pa = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2))/ denom;
		Pb = (Math.pow(x3, 2) * (y1 - y2) + Math.pow(x2, 2) * (y3 - y1) + Math.pow(x1, 2) * (y2 - y3))/ denom;
		Pc = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3)/ denom	;
		
	}
	/**
    * Given two points {(x1,y1),(x2,y2)} returns {m,b} from the equation y = mx + b for a line
    * which passes through both points.
    * 
    * @param x1 First x coordinate
    * @param y1 First y coordinate
    * @param x2 Second x coordinate
    * @param y2 Second y coordinate
    * @return The slope and y intercept of the line passing through the points provided
    */
   public double[] getLinearEquation() {
       m = (y2 - y1) / (x2 - x1);
       b = -(m * x1) + y1;
       double[] linearEquation=new double[] { m, b };
       
       return linearEquation;
   }
	
	
	
	public double getPointParabolaFunction() {
		System.out.println("Puntos: "+x1+" "+y1);
		System.out.println("Puntos: "+x2+" "+y2);
		System.out.println("Puntos: "+x3+" "+y3);
		System.out.print(Pa+"x^2 + "+Pb+"x + "+Pc);
		return 0;
	}
	
	/**	
	 * This function use the lagrange interpolation formula
	 * For the moment only work for three point (and thinking for 2d dimension)
	 * @param a1 x value of the point a
	 * @param a2 y value of the point a
	 * @param b1 x value of the point b
	 * @param b2 y value of the point b
	 * @param c1 x value of the point c
	 * @param c2 y value of the point c
	 * @return
	 */
	public double lagrangeInterpolationFunction(double x) {
		double y=0,p1,p2,p3;
		p1=y1*(((x-x2)*(x-x3))/((x1-x2)*(x1-x3)));
		p2=y2*(((x-x1)*(x-x3))/((x2-x1)*(x2-x3)));
		p3=y3*(((x-x1)*(x-x2))/((x3-x1)*(x3-x2)));
		y=p1+p2+p3;
		return y;
		
	
	}
	
	public double[] intersect(double lineM, double lineB) {
	    double x = (lineB - this.b) / (this.m - lineM);
	    double y = this.m * x + this.b;
	    return new double[] {x,y};
	}
	

	
	
	
	
	
	
	
	

}




