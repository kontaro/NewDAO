package MultiObj;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import IMRT_DAO.TreatmentPlan;

public class AttachedSurface {
	ArrayList<surface> surfaceList=new ArrayList<surface>();
	ArrayList<TreatmentPlan> treatmentPlanFront=new ArrayList<TreatmentPlan>();
	ArrayList<Integer> treatmentPlanIndexSurface=new ArrayList<Integer>();
	HashMap<Integer, surface> currentAttachedSurface=new HashMap<Integer, surface>();
	
	ArrayList<TreatmentPlan> intersectionListAttached=new ArrayList<TreatmentPlan>();
	
	public AttachedSurface(ArrayList<surface> initialSurfaces) {
		// TODO Auto-generated constructor stub
		surfaceList.addAll(initialSurfaces);
		mergeSurfaces2();
	}
	public AttachedSurface(surface initialSurface) {
		// TODO Auto-generated constructor stub
		surfaceList.add(initialSurface);
	}
	
	private void mergeInitialSurfaces(ArrayList<surface> initialSurfaces) {
		surfaceList.add(initialSurfaces.get(0));
		if(initialSurfaces.size()>=1) {
			for(int i=1;i<initialSurfaces.size();i++) {
				mergeSurfaces(initialSurfaces.get(i));
			}
		}
	}
	private boolean mergeSurfaces(surface newSurface) {
		/* If any point of the surface list is dominated a merge is required */
		ArrayList<Integer> dominatedSurfacesList=dominatedPointList(newSurface);
		if(dominatedSurfacesList.size()>0) {
			if(isDominated(newSurface))
			for(int i=0;i<surfaceList.size();i++) {
				mergeTwoSurfaces(surfaceList.get(dominatedSurfacesList.get(i)),newSurface);
				
			}
			return true;
		}
		else {
			/* If all the points of the new surface is non dominated the surface is add*/
			
			if(true){
				surfaceList.add(newSurface);
				return true;
			}
			else {
				return false;
			}
		}
		
	}
	
	
	private boolean isDominated(surface newSurface) {
		// TODO Auto-generated method stub
		return false;
	}
	private boolean mergeSurfaces2() {
		/* If any point of the surface list is dominated a merge is required */

		ArrayList<TreatmentPlan> intersectionList =mergeTwoSurfaces(surfaceList.get(0),surfaceList.get(1));
		intersectionListAttached.addAll(intersectionList);
		return false;
		
		
		
	}
	
	
	/**genera el tratamiento que se deberia encontrar con el porcentaje obtenido el valor Y del punto perperdicular
	 *  al punto interseccion dado con los puntos extremos de una surface
	 * @param intersectionPoint
	 * @param motherSurface
	 */
	public TreatmentPlan generateIntersectionTreatmentPlan(double[] intersectionPoint, surface motherSurface) {
		double[] extremePointUp=motherSurface.getExtreamPointTreatment(0);
		double[] extremePointdown=motherSurface.getExtreamPointTreatment(1);
		
		double[] linearEquation=getLinearEquation(extremePointUp,extremePointdown);
		double[] perpendicularPoint=encontrarPuntoPerpendicular(intersectionPoint[0], intersectionPoint[1], linearEquation[0], linearEquation[1]);
		double percent=calcularPorcentajeEnLinea(perpendicularPoint, extremePointUp, extremePointdown);
		
		return motherSurface.generatePointWithExtreme((1-percent));
	}
	
	
    private double calcularPorcentajeEnLinea(double[] punto, double[] puntoA, double[] puntoB) {
        // Calcular la longitud total de la línea
        double longitudTotal = calcularDistanciaEntrePuntos(puntoA, puntoB);

        // Calcular la longitud desde el puntoA hasta el punto evaluado
        double longitudDesdeA = calcularDistanciaEntrePuntos(puntoA, punto);

        // Calcular el porcentaje
        return (longitudDesdeA / longitudTotal);
    }

    // Función para calcular la distancia entre dos puntos
    private double calcularDistanciaEntrePuntos(double[] punto1, double[] punto2) {
        return Math.sqrt(Math.pow(punto2[0] - punto1[0], 2) + Math.pow(punto2[1] - punto1[1], 2));
    }
	
   public double[] getLinearEquation(double[] p1,double[] p2) {
       double m = (p2[1] - p1[1]) / (p2[0] - p1[0]);
       double b = -(m * p1[0]) + p1[1];
       double[] linearEquation=new double[] { m, b };
       
       return linearEquation;
   }
    // Función para encontrar el punto que forma un ángulo de 90 grados con la recta desde el punto A
    private double[] encontrarPuntoPerpendicular(double xA, double yA, double m, double b) {
        // Calcular la pendiente de la línea perpendicular
        double mPerpendicular = -1 / m;

        // Calcular las coordenadas del punto perpendicular
        
        double bPerpendicular=yA - (mPerpendicular * xA);
        
        double[] puntoPerpendicular = intersectionPoint(m, b, mPerpendicular, bPerpendicular);

 
        return intersectionPoint(m, b, mPerpendicular, bPerpendicular);
    }
	private double calcPercentWithRange(double value, double initialRange, double finalRange) {
            return ((value - initialRange) / (finalRange - initialRange)) ;
        
    }
	private ArrayList<TreatmentPlan> mergeTwoSurfaces(surface s1, surface s2) {
		ArrayList<TreatmentPlan> newExtremeArrayList = new ArrayList<TreatmentPlan>();
		equation eq21=s2.shape.get(0);//superior segment surface 2
		equation eq22=s2.shape.get(1);
		equation eq11=s1.shape.get(0);//superior segment surface 1
		equation eq12=s1.shape.get(1);
		double[] intersectionEq21Eq11=intersectionPoint(eq21,eq11);
		double[] intersectionEq21Eq12=intersectionPoint(eq21,eq12);
		
		
		double[] intersectionEq22Eq11=intersectionPoint(eq22,eq11);
		double[] intersectionEq22Eq12=intersectionPoint(eq22,eq12);
		
		
		if(validPoint(intersectionEq21Eq12, eq21)&&validPoint(intersectionEq21Eq11, eq21)) {
			TreatmentPlan interSectionPointEq21Eq11=generateIntersectionTreatmentPlan(intersectionEq21Eq11, s2);
			surface newSurface= new surface(s2.extremeSolution.get(0),interSectionPointEq21Eq11);
			//Aqui falta un if si es que el punto de s1 no es dominado se tiene que crear una superficie mas
			//la cual cubre desde el el up point de s1 hasta el up point de s2
			s1.changeUpExtreme(interSectionPointEq21Eq11);
			TreatmentPlan interSectionPointEq21Eq12=generateIntersectionTreatmentPlan(intersectionEq21Eq12, s2);
			s1.changeDownExtreme(interSectionPointEq21Eq12);
			s2.changeUpExtreme(interSectionPointEq21Eq12);
			newExtremeArrayList.add(interSectionPointEq21Eq12);
			return newExtremeArrayList;
			
		}
		if(validPoint(intersectionEq22Eq12, eq22)&&validPoint(intersectionEq22Eq11, eq22)) {
			TreatmentPlan interSectionPointEq21Eq11=generateIntersectionTreatmentPlan(intersectionEq21Eq11, s2);
			surface newSurface= new surface(s2.extremeSolution.get(0),interSectionPointEq21Eq11);
			//Aqui falta un if si es que el punto de s1 no es dominado se tiene que crear una superficie mas
			//la cual cubre desde el el up point de s1 hasta el up point de s2
			s1.changeUpExtreme(interSectionPointEq21Eq11);
			TreatmentPlan interSectionPointEq21Eq12=generateIntersectionTreatmentPlan(intersectionEq21Eq12, s2);
			s1.changeDownExtreme(interSectionPointEq21Eq12);
			s2.changeUpExtreme(interSectionPointEq21Eq12);		
			newExtremeArrayList.add(interSectionPointEq21Eq12);
			return newExtremeArrayList;
			
		}
		if(validPoint(intersectionEq21Eq12, eq21)&&!validPoint(intersectionEq21Eq11, eq21)) {
			TreatmentPlan interSectionPointEq21Eq12EndPointS2=generateIntersectionTreatmentPlan(intersectionEq21Eq12, s2);
			TreatmentPlan interSectionPointEq21Eq12EndPointS1=generateIntersectionTreatmentPlan(intersectionEq21Eq12, s1);
			s2.changeUpExtreme(interSectionPointEq21Eq12EndPointS2);
			s1.changeDownExtreme(interSectionPointEq21Eq12EndPointS1);
			newExtremeArrayList.add(interSectionPointEq21Eq12EndPointS1);
			newExtremeArrayList.add(interSectionPointEq21Eq12EndPointS2);
			return newExtremeArrayList;
		}
		
		return null;
		
	}
	
	/**
	 * return if a point is in the range of the limit of given equation
	 * @param point
	 * @param range
	 * @return
	 */
	private boolean validPoint(double[] point,equation range) {
		return (point[0] >= Math.min(range.x1, range.x2) && point[1] <= Math.max(range.y1, range.y2));
	
	}
	/**
	 * Return the x of the point that intersect two lines. If the lines are parallel return -1 because the objective function can not be negative
	 * y=MX+B
	 * @return
	 */
	private double[] intersectionPoint(equation e1, equation e2) {
        return intersectionPoint(e1.m,e1.b,e2.m,e2.b);
    }	
	/**
	 * Return the x of the point that intersect two lines. If the lines are parallel return -1 because the objective function can not be negative
	 * y=MX+B
	 * @return
	 */
	private double[] intersectionPoint(double m1, double b1, double m2, double b2) {
        if (m1 == m2) {
            return new double[] {-1,-1}; // Las rectas son paralelas y no tienen intersección
        } else {
            double x = (b2 - b1) / (m1 - m2);
            double y = m1 * x + b1;
            return new double[] {x,y};
        }
    }
	/**
	 * This function return wich surface is dominated by the newsurface
	 * If the newSurface is dominated for any surface the function return null value
	 * @param newSurface
	 * @return
	 */
	private ArrayList<Integer> dominatedPointList(surface newSurface) {
		// TODO Auto-generated method stub
		return null;
	}
	private ArrayList<Integer> dominatedPoint(surface newSurface) {
		// TODO Auto-generated method stub
		return null;
	}

	public void printSurfaces() {
		
	}
	public void printSurface(int index) {
		
	}
	
	
}
