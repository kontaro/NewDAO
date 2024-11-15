package SingleObj;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;
import java.util.Vector;

import IMRT_DAO.Aperture;
import IMRT_DAO.Beam;
import IMRT_DAO.DDM;
import IMRT_DAO.Organs;
import IMRT_DAO.TreatmentPlan;


/**
 * Clse que permite movimientos de vecindad multihilos
 * @author Mauricio Moyano
 *TODO
 *-Intercambiar los nombres de las variables state y direction
 */
public class TreatmentNeighborhood extends Thread {
	
	TreatmentPlan actualSolution;
	TreatmentPlan actual;
	TreatmentPlan bestSolution;
	ArrayList<TreatmentPlan> neigborhoodTreatment;
	public static Vector<int[][]> initial_aperture_shape;
	public int init_intensity;               // # DAO : Initial intensity for apertures
    public int max_intensity;                // # DAO : Maximum intensity for apertures
    public int max_delta;     
    public int max_iter;              		// # DAO : Cantidad m�xima de iteraciones del algoritmo
    public int max_time;  
    public int seed;                         // # DAO : semilla para partir la busqueda local
    public int step_intensity;	
    public int[] selAngles;
	public int numOrgans;
    public int numAngles; 
    public Organs[] o;
    public int [] dd;
    public DDM M;
    public int leaf; // #DAO: hoja en la que esta posicionado
    public boolean state; // #DAO: Direccion del movimiento Izquierda true  / Derecha False
    public int direction; // #DAO: Estado el movimiento 0 cierra/ 1 abre
    public boolean isBetter;
    public int called; // #number of time's called the run function
    public double scores[];//
    public double singleObjectiveValue;	
    public boolean solver;
    public  int [] zmax,zmin;
    public int idchangedBeamlet;
    public double changedBeamlet;
    public boolean refresh;
    public int angleChoose;
    public int apertureChoose;
    public boolean feasible;
    public int pChanged;
    
    public TreatmentNeighborhood(TreatmentPlan actual,int[] selAngles,Vector<int[][]> initial_aperture_shape,int angles,int numOrgans,DDM M,Organs[] o,int [] dd) {
    	actualSolution=actual;
    	this.singleObjectiveValue=actual.singleObjectiveValue;
    	TreatmentNeighborhood.initial_aperture_shape=initial_aperture_shape;
    	this.M=M;
    	this.o=o;
    	this.numAngles=angles;
    	this.numOrgans= numOrgans;
    	this.dd=dd;
    	neigborhoodTreatment=new ArrayList<TreatmentPlan>();
    	isBetter=false;
    	called=0;
    	refresh=false;
    	solver=true;
    	this.selAngles =selAngles ;
    	this.zmax=zmax;
    	this.zmax=zmin;
    	setPriority(MAX_PRIORITY);

    }
    public TreatmentNeighborhood(TreatmentPlan actual,int[] selAngles,int init, int max_i, int delta, int iter, int time, int seed, int step,Vector<int[][]> initial_aperture_shape,int angles,int numOrgans,DDM M,Organs[] o,int [] dd,int leaf,boolean state,int direction) {
    	this.actualSolution=new TreatmentPlan(init, max_i, delta, iter, time, seed, step, angles,numOrgans);
    	this.actualSolution.updateSolDAO(actual);
    	this.singleObjectiveValue=actual.singleObjectiveValue;
    	TreatmentNeighborhood.initial_aperture_shape=initial_aperture_shape;
    	this.M=M;
    	this.o=o;
    	this.numAngles=angles;
    	this.numOrgans= numOrgans;
    	this.dd=dd;
    	neigborhoodTreatment=new ArrayList<TreatmentPlan>();
    	this.leaf=leaf;
    	this.state=state;
    	this.direction=direction;
    	isBetter=false;
    	called=0;
    	solver=true;
    	refresh=false;
    	this.selAngles =selAngles ;
    	setPriority(MAX_PRIORITY);
    }
    public TreatmentNeighborhood(TreatmentPlan actual,int[] selAngles,int init, int max_i, int delta, int iter, int time, int seed, int step,Vector<int[][]> initial_aperture_shape,int angles,int numOrgans,DDM M,Organs[] o,int [] dd,int leaf,boolean state,int direction,int angleChoose, int apertureChoose) {
    	this.actualSolution=new TreatmentPlan(init, max_i, delta, iter, time, seed, step, angles,numOrgans);
    	this.actualSolution.updateSolDAO(actual);
    	this.actualSolution.evaluationMethod=actual.evaluationMethod;
    	this.actualSolution.evaluationFunction=actual.evaluationFunction;
    	this.singleObjectiveValue=actual.singleObjectiveValue;
    	TreatmentNeighborhood.initial_aperture_shape=initial_aperture_shape;
    	this.M=M;
    	this.o=o;
    	this.numAngles=angles;
    	this.numOrgans= numOrgans;
    	this.dd=dd;
    	neigborhoodTreatment=new ArrayList<TreatmentPlan>();
    	this.leaf=leaf;
    	this.state=state;
    	this.direction=direction;
    	isBetter=false;
    	called=0;
    	solver=true;
    	refresh=false;
    	this.selAngles =selAngles ;
    	setPriority(MAX_PRIORITY);
    	this.angleChoose=angleChoose;
    	this.apertureChoose=apertureChoose;

    }
    public TreatmentNeighborhood(TreatmentPlan actual,int[] selAngles,int init, int max_i, int delta, int iter, int time, int seed, int step,Vector<int[][]> initial_aperture_shape,int angles,int numOrgans,DDM M,Organs[] o,int [] dd,int angleChoose) {
    	this.actualSolution=new TreatmentPlan(init, max_i, delta, iter, time, seed, step, angles,numOrgans);
    	this.actualSolution.updateSolDAO(actual);
    	this.singleObjectiveValue=actual.singleObjectiveValue;
    	TreatmentNeighborhood.initial_aperture_shape=initial_aperture_shape;
    	this.M=M;
    	this.o=o;
    	this.numAngles=angles;
    	this.numOrgans= numOrgans;
    	this.dd=dd;
    	neigborhoodTreatment=new ArrayList<TreatmentPlan>();
    	isBetter=false;
    	called=0;
    	solver=true;
    	refresh=false;
    	this.selAngles =selAngles ;
    	setPriority(MAX_PRIORITY);
    	this.angleChoose=angleChoose;
  

    }
    public void id() {
    	System.out.println("id Solver: "+leaf+state+direction);
    }
    public int[] getApertureChanged() {
    	return actualSolution.apertures.get(angleChoose).get(apertureChoose).aperture[leaf];
    }
    public int[] getApertureChanged2() {
    	return actual.apertures.get(angleChoose).get(apertureChoose).aperture[leaf];
    }
    public String idstring() {
    	String id="";
    	id=id+angleChoose+apertureChoose+leaf+state+direction;
    	return id;
    }
    
 
    
    
    public void generateNeighborhood() {
		if(!neigborhoodTreatment.isEmpty()) actual=neigborhoodTreatment.get(0); //significa que en la anterior generacion no mejoro
		neigborhoodTreatment.clear();
		int positionChanged=0;
		idchangedBeamlet=-1;
		
		if(actualSolution.apertures.size()>5) {
			System.out.println("error 11");
		}
		
		int[][] list_shape = initial_aperture_shape.get(angleChoose);
		//int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		//for(int i=0;i<sizeAperture;i++){
			
		//int[][] aperture = listApertures.get(apertureChoose).aperture;
		TreatmentPlan neigbor=new TreatmentPlan(numAngles,numOrgans);
		neigbor.updateSolDAO(actualSolution);
		neigbor.evaluationMethod=actualSolution.evaluationMethod;
		neigbor.evaluationFunction=actualSolution.evaluationFunction;
		
		positionChanged=neigbor.apertures.get(angleChoose).get(apertureChoose).changeRow2(leaf, direction, state,list_shape[leaf][0],list_shape[leaf][1]);
		pChanged=positionChanged;
		//neigbor.changeAperture(aperture, angleChoose, apertureChoose);
		if(pChanged==-1) {
			pChanged=positionChanged;
		}
		
		//int[][] algo= neigbor.beamletCoord();
		if(isFeasible(neigbor) && positionChanged>0 ){
			positionChanged=positionChanged;
			neigborhoodTreatment.add(neigbor);
			idchangedBeamlet=neigbor.searchBeamletID(angleChoose,leaf,positionChanged);
			//Se guarda el valor de la intensidad anterior para luego evaluar el cambio
			changedBeamlet=neigbor.intensity[idchangedBeamlet];
			
			
			
			//se actualiza el valor de la intensidad 
			if(direction==1) {
				
				//double t=neigbor.intensity[idchangedBeamlet]+neigbor.apertures.get(angleChoose).get(apertureChoose).intensity;
				neigbor.intensity[idchangedBeamlet]=neigbor.intensity[idchangedBeamlet]+neigbor.apertures.get(angleChoose).get(apertureChoose).intensity;
				
			}
			else
				neigbor.intensity[idchangedBeamlet]=neigbor.intensity[idchangedBeamlet]-neigbor.apertures.get(angleChoose).get(apertureChoose).intensity;

		}

	}
    
    public void generateNeighborhood1() {



		boolean flag=true;
		bestSolution=new TreatmentPlan(numAngles,numOrgans);
		bestSolution.updateSolDAO(actualSolution);
		bestSolution.evaluationMethod=actualSolution.evaluationMethod;
		bestSolution.evaluationFunction=actualSolution.evaluationFunction;
		int[][] list_shape = initial_aperture_shape.get(angleChoose);
		int positionChanged=0;
		int beamletMo;
		double[] oldIntensity;
		while(flag) {
			flag=false;
			positionChanged= actualSolution.apertures.get(angleChoose).get(apertureChoose).changeRow2(leaf, direction, state,list_shape[leaf][0],list_shape[leaf][1]);

			beamletMo=actualSolution.searchBeamletID(angleChoose,leaf,positionChanged);	
			oldIntensity= this.bestSolution.intensity.clone();
			//intensityChangedBeamlet=bestSolution.intensity[beamletMo];
			//se hace el movimiento
			//se revisa si es factible
			
			if(isFeasible(actualSolution)&&positionChanged!=-3) {
				//long startTime = System.currentTimeMillis();
				//double value = actualSolution.evaluateSolution();
				//double value = actualSolution.evaluateSolution(beamletMo,intensityChangedBeamlet);
				
	
				double value = actualSolution.deltaEvaluateSolution(beamletMo,oldIntensity);
				//System.out.println("Tiempo_evaluate_mauriC: " + totalTime);
				//long totalTime = System.currentTimeMillis() - startTime;
		        
		        //System.out.println("Tiempo value: " + totalTime);
				if(value < bestSolution.singleObjectiveValue) {
					flag=true;
					bestSolution.updateSolDAO(actualSolution);
					//System.out.println("se hace update");
				}
			}else {
				//System.out.println("no es factible en intensificacion");
			}
		}
		

	}
    
    
    
    
	public void generateNeighborhood(int angleChoose,int apertureChoose) {
		if(!neigborhoodTreatment.isEmpty()) actual=neigborhoodTreatment.get(0); //significa que en la anterior generacion no mejoro
		neigborhoodTreatment.clear();
		this.angleChoose=angleChoose;
		this.apertureChoose=apertureChoose;
		int positionChanged=0;
		idchangedBeamlet=-1;
		int[][] list_shape = initial_aperture_shape.get(angleChoose);
		//int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		//for(int i=0;i<sizeAperture;i++){
			
		//int[][] aperture = listApertures.get(apertureChoose).aperture;
		TreatmentPlan neigbor=new TreatmentPlan(numAngles,numOrgans);
		neigbor.updateSolDAO(actualSolution);
		
		positionChanged=neigbor.apertures.get(angleChoose).get(apertureChoose).changeRow2(leaf, direction, state,list_shape[leaf][0],list_shape[leaf][1]);
		if(positionChanged==0 && angleChoose!=0) {
			pChanged=positionChanged;
		}
		
		if(positionChanged==-1 ||positionChanged==-2) {
			pChanged=positionChanged;
		}
		//neigbor.changeAperture(aperture, angleChoose, apertureChoose);
		//int[][] algo= neigbor.beamletCoord();
		if(isFeasible(neigbor) && positionChanged>0 ){
			positionChanged=positionChanged;
			neigborhoodTreatment.add(neigbor);
			idchangedBeamlet=neigbor.searchBeamletID(angleChoose,leaf,positionChanged);
			//Se guarda el valor de la intensidad anterior para luego evaluar el cambio
			changedBeamlet=neigbor.intensity[idchangedBeamlet];
			
			
			
			//se actualiza el valor de la intensidad 
			if(direction==1) {
				
				//double t=neigbor.intensity[idchangedBeamlet]+neigbor.apertures.get(angleChoose).get(apertureChoose).intensity;
				neigbor.intensity[idchangedBeamlet]=neigbor.intensity[idchangedBeamlet]+neigbor.apertures.get(angleChoose).get(apertureChoose).intensity;
				
			}
			else
				neigbor.intensity[idchangedBeamlet]=neigbor.intensity[idchangedBeamlet]-neigbor.apertures.get(angleChoose).get(apertureChoose).intensity;
			
			

			
			//positionChanged=positionChanged;
			//double[] w = new double[]{1,1,1}; 
			//neigbor.aperturesToIntensityMatrixPerStation();
			//neigbor.intensityMatrixStationToIntensityMatrixPlan();
			//deltaEvaluateSolutionpoof(neigbor, M, o, w);
			//System.out.println(neigbor.singleObjectiveValue);
			//evaluateSolution(neigbor, M, o, w);
			//System.out.println(neigbor.singleObjectiveValue);
			

		}

		
	
	}

	
	@Override
	/**
	 * TODO:
	 * Cambiar el w de pesos que esta declarado dentro de la funcion
	 */
	public void run() {
		generateNeighborhood1();
	}
	

	public static boolean isFeasible(TreatmentPlan sol){
		Vector<List<Aperture>> stations = sol.apertures;
		boolean feasible=true;
		//recorrido de angulos
		for(int i=0; i<stations.size(); i++){
			
			List<Aperture> list_apertures = stations.get(i);
			int[][] limit_shape = initial_aperture_shape.get(i);
			//recorrido de aperturas por angulo
			for(int j=0; j<list_apertures.size(); j++){
				//int intensity = list_apertures.get(j).intensity;
				int[][] matrix = list_apertures.get(j).aperture;
				
				for(int x=0;x<matrix.length; x++){
					if((limit_shape[x][0]>matrix[x][0]&& matrix[x][1]!=-1) || (limit_shape[x][1]<matrix[x][1]&&matrix[x][0]!=-1) ){
						//System.out.println(limit_shape[x][0]+ " " + matrix[x][0]+ " "+limit_shape[x][1]+ " " + matrix[x][1]);
						feasible=false;
						return feasible; 
					}
					if(matrix[x][0]>matrix[x][1]) {
						//System.out.println(limit_shape[x][0]+ " " + matrix[x][0]+ " "+limit_shape[x][1]+ " " + matrix[x][1]);
						feasible=false;
						return feasible; 
					}

					
				}	
			}
				
		}
		return feasible;
	}

	
}
