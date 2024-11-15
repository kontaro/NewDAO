package ApertureBase;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import IMRT_DAO.Aperture;
import IMRT_DAO.DDM;
import IMRT_DAO.Organs;
import IMRT_DAO.TreatmentPlan;
import gurobi.GRBException;

public class sequentialGenerator {
	public Vector<int[][]> apertureShapes;
	public  double [] fmo;
	int[] beamletAngles;
	public sequentialGenerator(int totalbmlt,DDM M, Organs[] o, double[] dd,int[] beamletAngles,int numAperture,int rounded, int numAngles, TreatmentPlan solution) {
		
		this.beamletAngles=beamletAngles;
		
		
		
		try {
			FMO_Solver newModel=new FMO_Solver(totalbmlt,M,o,beamletAngles,dd,solution.weights);

			System.out.println("OBJ: "+newModel.objVal);
			fmo=newModel.newIntensity;
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		MLC obt=new MLC();
		roundBMLTS(rounded);
		
		setInitialApertures( solution, fmo,numAngles);
		
	    for(int j=0;j<apertureShapes.size();j++) {
	    	List<Aperture> list_aperture = new ArrayList<Aperture>();
	    	obt.getAperture(apertureShapes.get(j));
	    	//setInitialSolution(solution,numAperture,j,obt,list_aperture);
	    	setAllSolution(solution,numAperture,j,obt,list_aperture);
	    	solution.addAperture(list_aperture);
	    }
	   
		
	}
	
	
	
	private void setInitialSolution(TreatmentPlan solution, int numAperture, int angle, MLC obt, List<Aperture> list_aperture) {
		Vector<double[]> v =indexSQFperBeam(1, 1, angle, obt);
		orderVector(v);
		
		//obtenermos el factor SQF y agregamos los numAperture mejores al tratamiento en el angulo angle
		for(int i=0;i<numAperture;i++) {
			list_aperture.add(obt.decompositionAperture.get((int) v.get(i)[0]));
			
		}
		// TODO Auto-generated method stub
		
	}
	private void setAllSolution(TreatmentPlan solution, int numAperture, int angle, MLC obt, List<Aperture> list_aperture) {
		
		
		//obtenermos el factor SQF y agregamos los numAperture mejores al tratamiento en el angulo angle
		for(int i=0;i<obt.decompositionAperture.size();i++) {
			list_aperture.add(obt.decompositionAperture.get(i));
			
		}
		// TODO Auto-generated method stub
		
	}

	 public double segmentQualityFactor(float segmentWeight,float intensityWeight,MLC obt,int aperturePosition, int angle,int segmentTotal, double weightMax) {
    	 double SQF=0; 
    	 int segmentith=obt.decompositionAperture.get(aperturePosition).getTotalApertureOpen(); 
    	 double weightith=obt.decompositionAperture.get(aperturePosition).intensity; 
    	 SQF=(segmentWeight*(segmentith/segmentTotal))+(intensityWeight*(weightith/weightMax));
    	 return SQF;
    }
	 
	  public Vector<double[]> indexSQFperBeam(float segmentWeight,float intensityWeight, int beamPosition,MLC obt) {
		    Vector<double[]> v = new Vector<double[]>();
		    int segmentTotal=0;
		    segmentTotal=beamletAngles[beamPosition];
		    double SQF=0;   
		   	double weightMax=obt.maxIntensity;
		   
		   	for(int j=0; j<obt.decompositionAperture.size(); j++){
		   		double[] SQFBeam=new double[2];
		   		SQF=segmentQualityFactor(segmentWeight, intensityWeight,obt, beamPosition, j, segmentTotal, weightMax);
				//list_apertures.get(j).intensity=newIntensity[i][j];
		   		SQFBeam[0]=j;
		   		SQFBeam[1]=SQF;
		   		v.add(SQFBeam);
		   		
				
			}
	  
	   	 return v;
	   }

	public void roundBMLTS(int mod) {
		double [] x;
		x=new double[fmo.length];
		for(int i =0;i<fmo.length;i++) {
			double aux=fmo[i]/mod;
			aux=Math.round(aux);
			x[i]=aux*mod;
		}
		fmo=x;
		
	}
	public void orderVector(Vector<double[]> v) {
		
		int n=v.size();
		double aux[]=new double[2];
		for(int i=0;i<n;i++) {
			for(int j=1; j < (n-i); j++){ 
				if(v.get(j-1)[1]>v.get(j)[1]) {
					aux [0]= v.get(j-1)[0];  
					aux [1]= v.get(j-1)[1];
					v.get(j-1)[0] = v.get(j)[0];  
					v.get(j-1) [1]= v.get(j)[1]; 
					v.get(j)[0] = aux[0]; 
					v.get(j)[1] = aux[1] ;  
				}
			}
		}
	}
	
	
	public void setInitialApertures(TreatmentPlan solution,double[] solutionVector, int numAngles){
		int sizeX, sizeY, x, y;
		apertureShapes=new Vector<int[][]>();
		
		int beamlets=0;

		for(int i=0;i<numAngles;i++){
			sizeX  = solution.selAngles[i].maxAbsX;
			sizeY = solution.selAngles[i].maxAbsY;
			int [][] apertureShape = new int[sizeX][sizeY];
			for(x=0; x<sizeX; x++){
				for(y=0; y<sizeY; y++){
					apertureShape[x][y] = 0;
				}
			}
			
			double[][] beamletsCoord = solution.selAngles[i].beamletsCoord;
			for(int j=0;j<solution.selAngles[i].beamlets;j++){
				int xcoord=(int)beamletsCoord[j][1]-1;
				int ycoord=(int)beamletsCoord[j][2]-1;
				apertureShape[xcoord][ycoord] =(int)solutionVector[beamlets] ;
				beamlets=beamlets+1;	
			}
			apertureShapes.add(apertureShape);
			

		}
	}


}
