package SingleObj;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import IMRT_DAO.Aperture;
import IMRT_DAO.DDM;
import IMRT_DAO.Organs;
import IMRT_DAO.TreatmentPlan;
import gurobi.GRBException;

public class Algorithm {
	public int numOrgans;
    public int numAngles;                   //#Beams in TreatmentPlan
    public Organs[] o;
    public String pathFile = "";
	public int num_apertures = 5;            // # DAO : Number of apertures per station
	public int init_intensity;               // # DAO : Initial intensity for apertures
    public int max_intensity;                // # DAO : Maximum intensity for apertures
    public int max_delta;                    // # DAO : Delta máximo para la variación de intesidad por apertura
    public int max_iter;              		// # DAO : Cantidad máxima de iteraciones del algoritmo
    public int max_time;                     // # DAO : tiempo máximo de ejecución del algoritmo
    public int seed;                         // # DAO : semilla para partir la busqueda local
    public int step_intensity;	
	public DDM M;
	public int nThreads=12;
    public int[][] beamletSet;
    public int[] initialBACs;
    public Vector<int[][]> initial_aperture_shape;  // # DAO : Conjunto de formas de las aperturas según restricciones (rango sup y rango inf)
    public Vector<int[][]> apertureShapes;
	 int[] beamletAngles;
	 int[] dd;
	int[][] apertureCombination;   // numero de combinaciones que hay en una apertura
    
	public int iterations; //# DAO: Numero de iteraciones que a hecho el algoritmo
	
	
	public Algorithm() {
		
	}
	
	public Algorithm(int[][] beamletSet,int[] initialBACs,DDM M, Organs[] o, double[] dd, int[] beamletAngles, int init_intensity, int max_intensity,
			int max_delta, int max_iter, int max_time, int seed, int step_intensity, int numAngles, int numOrgans) {
	
		this.M=M;
		this.o=o;
		
		this.beamletAngles=beamletAngles;
		this.init_intensity=init_intensity;
		this.max_intensity=max_intensity;
		this.max_delta=max_delta;
		this.max_iter=max_iter;
		this.max_time=max_time; 
		this.seed=seed; 
		this.step_intensity=step_intensity; 
		this.numAngles=numAngles; 
		this.numOrgans=numOrgans;
		this.initialBACs=initialBACs;
		this.beamletSet=beamletSet;

		this.dd=doubletoint(dd);
		apertureCombination=null;
		combination();
	}

	
	//------------Algoritmos de Optimizacion Principales
	
	
	//----------Algoritmos de optimisacion de formas
	/**
	 * Funcion que llama a local search con solver
	 * @param M es la ddm
	 * @param solution es la solucion inicial
	 * @param moviment
	 * @return
	 * TODO:
	 * dejar las variables w y selAngles como parametros de la funcion
	 */
	public TreatmentPlan LSmath(DDM M,TreatmentPlan solution) throws IOException, GRBException {
	       
		setInitialApertureShape(solution); //guardamos la figura
	    
	    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
	    solution.evaluateSolution();
	   
	    System.out.print("Initial: "+solution.singleObjectiveValue+" ");

	    solution.updateSolDAO(LocalSearch(solution,M,2));

	    solution.aperturesToIntensityMatrixPerStation();

	    solution.intensityMatrixStationToIntensityMatrixPlan();

	    return solution;
	}
	
	/**
	 * Funcion que llama a local search con solver
	 * @param M es la ddm
	 * @param solution es la solucion inicial
	 * @param moviment
	 * @return
	 * TODO:
	 * dejar las variables w y selAngles como parametros de la funcion
	 */
	public TreatmentPlan LSPen(DDM M,TreatmentPlan solution) throws IOException, GRBException {
	       
		setInitialApertureShape(solution); //guardamos la figura
	    
	    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
	    solution.evaluateSolution();
	   
	    //System.out.print("Initial: "+solution.singleObjectiveValue+" ");

	    solution.updateSolDAO(LocalSearchPen(solution,M,2));

	    solution.aperturesToIntensityMatrixPerStation();

	    solution.intensityMatrixStationToIntensityMatrixPlan();

	    
	    //solution.print();
	    return solution;
	}
	
	
	
	
	
	/**
	 * Funcion que llama a local search con solver
	 * @param M es la ddm
	 * @param solution es la solucion inicial
	 * @param moviment
	 * @return
	 * TODO:
	 * dejar las variables w y selAngles como parametros de la funcion
	 */
	public TreatmentPlan rVNS(DDM M,TreatmentPlan solution) throws IOException, GRBException {
	       
		setInitialApertureShape(solution); //guardamos la figura
	    
	    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
	    solution.evaluateSolution();
	   
	    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
	    //solution.updateSolDAO(LocalSearchWithSolver(solution,M));
	    //solution.updateSolDAO(LocalSearchTwoStep(solution,M));
	    solution.updateSolDAO(reducedVariableNeighbourhoodSearch(solution,M));
	    //System.out.println("value:"+solution.singleObjectiveValue);
	    solution.aperturesToIntensityMatrixPerStation();

	    solution.intensityMatrixStationToIntensityMatrixPlan();
		evaluateSolution(solution, M, o,solution.weights);
	    return solution;
	}
	public TreatmentPlan VND(DDM M,TreatmentPlan solution) throws IOException, GRBException {
	       
		setInitialApertureShape(solution); //guardamos la figura
	    
	    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
	    solution.evaluateSolution();
	   
	    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
	    //solution.updateSolDAO(LocalSearchWithSolver(solution,M));
	    //solution.updateSolDAO(LocalSearchTwoStep(solution,M));
	    solution.updateSolDAO(VariableNeigborhoodDescent2(solution,M));
	    //System.out.println("value:"+solution.singleObjectiveValue);
	    solution.aperturesToIntensityMatrixPerStation();

	    solution.intensityMatrixStationToIntensityMatrixPlan();
		evaluateSolution(solution, M, o,solution.weights);
	    return solution;
	}

	/**
	 * Funcion que llama a local search con solver
	 * @param M es la ddm
	 * @param solution es la solucion inicial
	 * @param moviment
	 * @return
	 * TODO:
	 * dejar las variables w y selAngles como parametros de la funcion
	 */
	public TreatmentPlan GA(DDM M,TreatmentPlan solution) throws IOException, GRBException {
	       
		setInitialApertureShape(solution); //guardamos la figura
	    //initialPopulation();
	    //TreatmentNeighborhood
		//solution.roundIntensity();
	    //solution.aperturesToIntensityMatrixPerStation();
	    //solution.intensityMatrixStationToIntensityMatrixPlan();
	    solution.beamletCoord();
	    evaluateSolution(solution, M, o,solution.weights);
	    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
	    //solution.updateSolDAO(LocalSearchWithSolver(solution,M));
	    //solution.updateSolDAO(LocalSearchTwoStep(solution,M));
	    solution.updateSolDAO(LocalSearch(solution,M,2));
	    //System.out.println("value:"+solution.singleObjectiveValue);
	    solution.aperturesToIntensityMatrixPerStation();

	    solution.intensityMatrixStationToIntensityMatrixPlan();
		evaluateSolution(solution, M, o,solution.weights);
		//solution.printIntensity();
		//solution.solverIntensity2();
		//solution.printIntensity();
		//System.out.println("final value Revision:"+solution.singleObjectiveValue);
	    //printAperture(solution);
	    //printFirstSolution(solution);
	    
	    //solution.print();
	    return solution;
	}
	//---------algoritmos de optimizacion de intensidades
	
	
	
	
	//---------Funciones de apoyo
	
	
	/**
	 * Es una busqueda local que al generar un vecino lo prueba lo evalua con el solver
	 * @param sol
	 * @param M
	 * @param conditionCase
	 * @param conditionValue
	 * @param moviment
	 * @return
	 * @throws IOException
	 */
	public  TreatmentPlan LocalSearchWithSolver(TreatmentPlan sol,DDM M) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		ArrayList<TreatmentPlan> neigborhood = new ArrayList<TreatmentPlan>();
		//TreatmentNeighborhood neigborhoodnew;
		int mov=1;
		iterations=0;
	
		boolean condition=true;
		boolean update=false;

		
	
		while(condition) {
			//BestSolution.print();
			iterations=iterations+1;
			update=false;
			//---------------Se genera la vecindad
			neigborhood.clear();
			//neigborhoodnew=new TreatmentNeighborhood(BestSolution,initial_aperture_shape,numAngles,numOrgans,M,o,dd);
			//movTest(neigborhood,BestSolution,M);
			//generateNeighborhoodPoolwithoutBest(neigborhood,BestSolution,M);
			//generateNeighborhoodPoolMix(neigborhood,BestSolution,M);
			
			
			
			
			
			//movAperture(neigborhood,BestSolution,M);
			//movIntensification(neigborhood,BestSolution,M);
			
			
			//movOnlyLeaf(neigborhood,BestSolution,M);
			movLeafRandom(neigborhood,BestSolution,M,20);
			//movCrossover(neigborhood,BestSolution,M);
			
			
			
			//generateNeighborhoodPoolSolver(neigborhood,BestSolution,M);
			//generateNeighborhoodPoolMixPrueba(neigborhood,BestSolution,M);
			//neigborhoodnew.generateNeighborhood();
			//neigborhood=neigborhoodnew.neigborhoodTreatment;
			//neigborhoodnew.start();
			//System.out.println("entre ");
			//System.out.println("tamaño"+neigborhood.size());
			

			if(neigborhood.size()>0) {
	
				actualSolution.updateSolDAO(bestNeigborPool(neigborhood,M));
				
				//Cuando se obtiene el mejor resultado se optimizan la intensidad
				
				
				//System.out.println(actualSolution.singleObjectiveValue);
				//System.out.println("entre ");
				if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
					BestSolution.updateSolDAO(actualSolution);
					update=true;
					//printIntensitiesAperture(BestSolution);
					//printIntensitiesAperture(BestSolution);
					//BestSolution.print();
					//System.out.println(BestSolution.singleObjectiveValue);
					
					
					
				}
				if(!update) {
					condition=false;
				}
				
			}else {
				condition=false;
			}
			//System.out.println("iterationTime: "+(System.currentTimeMillis()-timeIteration));

		}
		
		return BestSolution;
	}
	
	
	public  TreatmentPlan LocalSearchTwoStep(TreatmentPlan sol,DDM M) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		ArrayList<TreatmentPlan> neigborhood = new ArrayList<TreatmentPlan>();
		//TreatmentNeighborhood neigborhoodnew;
		int mov=1;
		iterations=0;
	
		boolean condition=true;
		boolean update=false;

		int flag=0;
	
		while(condition) {
			//BestSolution.print();
			iterations=iterations+1;
			update=false;
			//---------------Se genera la vecindad
			neigborhood.clear();
			//neigborhoodnew=new TreatmentNeighborhood(BestSolution,initial_aperture_shape,numAngles,numOrgans,M,o,dd);
			//movTest(neigborhood,BestSolution,M);
			//generateNeighborhoodPoolwithoutBest(neigborhood,BestSolution,M);
			//generateNeighborhoodPoolMix(neigborhood,BestSolution,M);
			
			
			
			
			
			//movAperture(neigborhood,BestSolution,M);
			//movIntensification(neigborhood,BestSolution,M);
			
			
			//movOnlyLeaf(neigborhood,BestSolution,M);
			movLeafRandom(neigborhood,BestSolution,M,20);
			//movCrossover(neigborhood,BestSolution,M);
			
			
			
			//generateNeighborhoodPoolSolver(neigborhood,BestSolution,M);
			//generateNeighborhoodPoolMixPrueba(neigborhood,BestSolution,M);
			//neigborhoodnew.generateNeighborhood();
			//neigborhood=neigborhoodnew.neigborhoodTreatment;
			//neigborhoodnew.start();
			//System.out.println("entre ");
			//System.out.println("tamaño"+neigborhood.size());
			

			if(neigborhood.size()>0) {
	
				actualSolution.updateSolDAO(bestNeigborPool(neigborhood,M));
				
				//Cuando se obtiene el mejor resultado se optimizan la intensidad
				
				
				//System.out.println(actualSolution.singleObjectiveValue);
				//System.out.println("entre ");
				if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
					BestSolution.updateSolDAO(actualSolution);
					update=true;
					flag=0;
					//printIntensitiesAperture(BestSolution);
					//printIntensitiesAperture(BestSolution);
					//BestSolution.print();
					System.out.println(BestSolution.singleObjectiveValue);
					
					
					
				}else {
					System.out.println(BestSolution.singleObjectiveValue);
				}
				if(!update) {
					condition=false;
					actualSolution.updateSolDAO(BestSolution);
					System.out.println();
					System.out.println(actualSolution.singleObjectiveValue);
					
					evaluateSolution(actualSolution, M, o, actualSolution.weights);
					System.out.println(actualSolution.singleObjectiveValue);
					solverIntensity(BestSolution);
					System.out.println(BestSolution.singleObjectiveValue);
					BestSolution.aperturesToIntensityMatrixPerStation();
					BestSolution.intensityMatrixStationToIntensityMatrixPlan();
					evaluateSolution(BestSolution, M, o, BestSolution.weights);
					
					if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
						BestSolution.updateSolDAO(actualSolution);
						update=true;
					}
					
				}
				
			}else {
				if(flag==0) {
					//Aqui deberia ir un algoritmo de reparacion el cual perturbe la solcuion
					solverIntensity(BestSolution);
					flag=1;
				}
				else
					condition=false;
				
			}
			//System.out.println("iterationTime: "+(System.currentTimeMillis()-timeIteration));

		}
		
		return BestSolution;
	}
	
	
	public  TreatmentPlan LocalSearchTest(TreatmentPlan sol,DDM M) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		ArrayList<TreatmentPlan> neigborhood = new ArrayList<TreatmentPlan>();
		//TreatmentNeighborhood neigborhoodnew;
		int mov=1;
		iterations=0;
	
		boolean condition=true;
		boolean update=false;

		int flag=0;
	
		while(condition) {
			//BestSolution.print();
			iterations=iterations+1;
			update=false;
			//---------------Se genera la vecindad
			neigborhood.clear();
			//neigborhoodnew=new TreatmentNeighborhood(BestSolution,initial_aperture_shape,numAngles,numOrgans,M,o,dd);
			//movTest(neigborhood,BestSolution,M);
			//generateNeighborhoodPoolwithoutBest(neigborhood,BestSolution,M);
			//generateNeighborhoodPoolMix(neigborhood,BestSolution,M);
			
			
			
			
			
			//movAperture(neigborhood,BestSolution,M);
			//movIntensification(neigborhood,BestSolution,M);
			
			
			//movOnlyLeaf(neigborhood,BestSolution,M);
			//movLeafRandom(neigborhood,BestSolution,M,24);
			//movCrossover(neigborhood,BestSolution,M);
			
			
			
			//generateNeighborhoodPoolSolver(neigborhood,BestSolution,M);
			//generateNeighborhoodPoolMixPrueba(neigborhood,BestSolution,M);
			//neigborhoodnew.generateNeighborhood();
			//neigborhood=neigborhoodnew.neigborhoodTreatment;
			//neigborhoodnew.start();
			//System.out.println("entre ");
			//System.out.println("tamaño"+neigborhood.size());
			
			
			
			//Pruebas algoritmos finales
			movLeafsAperture(neigborhood,BestSolution,M);
			//movLeafsRandom(neigborhood,BestSolution,M,24);
			

			if(neigborhood.size()>0 ) {
	
				actualSolution.updateSolDAO(bestNeigborPool(neigborhood,M));
				
				//Cuando se obtiene el mejor resultado se optimizan la intensidad
				
				
				//System.out.println(actualSolution.singleObjectiveValue);
				//System.out.println("entre ");
				if((actualSolution.singleObjectiveValue-BestSolution.singleObjectiveValue) < (-0.01) ){
					
					BestSolution.updateSolDAO(actualSolution);
					update=true;
					flag=0;
					//printIntensitiesAperture(BestSolution);
					//printIntensitiesAperture(BestSolution);
					//BestSolution.print();
					//System.out.println(BestSolution.singleObjectiveValue);
					
					
					
				}else {
					//System.out.println(BestSolution.singleObjectiveValue);
				}
				if(!update) {
					condition=false;
					actualSolution.updateSolDAO(BestSolution);
					//System.out.println();
					//System.out.println(actualSolution.singleObjectiveValue);
					actualSolution.aperturesToIntensityMatrixPerStation();
					actualSolution.intensityMatrixStationToIntensityMatrixPlan();
					evaluateSolution(actualSolution, M, o, actualSolution.weights);
					//System.out.println(actualSolution.singleObjectiveValue);
					solverIntensity(actualSolution);
					//System.out.println(actualSolution.singleObjectiveValue);
					actualSolution.aperturesToIntensityMatrixPerStation();
					actualSolution.intensityMatrixStationToIntensityMatrixPlan();
					evaluateSolution(actualSolution, M, o, actualSolution.weights);
					
					
					
					
					
					if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
						BestSolution.updateSolDAO(actualSolution);
						update=true;
						condition=true;
					}
					
				}
				
			}else {
				if(flag==0) {
					//Aqui deberia ir un algoritmo de reparacion el cual perturbe la solcuion
					
					solverIntensity(BestSolution);
					flag=1;
					
					
				}
				else
					condition=false;
				
			}
			//System.out.println("iterationTime: "+(System.currentTimeMillis()-timeIteration));

		}
		
		return BestSolution;
	}
	
	public int [][]getStateAperture(){
		int[][] stateAperture=new int[numAngles][num_apertures];
		for(int i=0;i<numAngles;i++) {
			for(int j=0;j<num_apertures;j++) {
				stateAperture[i][j]=0;
			}
		}
		return stateAperture;
		
	}
	public int [][]getLowAperture(int delta){
		int[][] stateAperture=new int[numAngles][num_apertures];
		for(int i=0;i<numAngles;i++) {
			for(int j=0;j<num_apertures;j++) {
				stateAperture[i][j]=0;
			}
		}
		return stateAperture;
		
	}
	
	public  TreatmentPlan LocalSearch(TreatmentPlan sol,DDM M, int mov) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		ArrayList<TreatmentPlan> neigborhood = new ArrayList<TreatmentPlan>();
		int[] apertureChanged =new int[2];
		int[][] stateAperture=getStateAperture();

		iterations=0;
	
		boolean condition=true;
		boolean update=false;
		int flag=0;

		BestSolution.getListLowIntensity(1);
		while(condition && flag<=1) {
			iterations=iterations+1;
			update=false;
			//---------------Se genera la vecindad
			neigborhood.clear();


			switch(mov) {
			
			case 1:
				movAperture2(neigborhood,BestSolution,M);
				//apertureChanged=movLeafsAperture(neigborhood,BestSolution,M);
				break;
				
			case 2:
				movLeafsRandom(neigborhood,BestSolution,M,24);
				break;
				
			case 3:
				movLeafsAperture(neigborhood,BestSolution,M);
				movLeafsRandom(neigborhood,BestSolution,M,24);
			
			default:
				apertureChanged=movLeafsAperture(neigborhood,BestSolution,M);
				break;
			
			}
			

			if(neigborhood.size()>0 ) {
	
				actualSolution.updateSolDAO(bestNeigborPool(neigborhood,M));
				//Si al optimizar la apertura producio mejora hace esto
				if((actualSolution.singleObjectiveValue-BestSolution.singleObjectiveValue) < (-0.01) ){
					
					BestSolution.updateSolDAO(actualSolution);

					update=true;
					
					
				}else {
					
				//En el caso contrario se actualiza la solucion con el solver con delta
					if(!update) {
						condition=false;
						actualSolution.updateSolDAO(BestSolution);
						actualSolution.solverIntensity(1);
				
						
						
						
						//si al usar el solver la solucion mejora se vuelve a optimizar las aperturas
						if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
							BestSolution.updateSolDAO(actualSolution);
							if(mov==1)BestSolution.getListLowIntensity(1);
							update=true;
							condition=true;
						}
						else {
							flag+=1;
							condition=true;
						}
					}					
				}	
			}
			//System.out.println(iterations+" "+BestSolution.singleObjectiveValue);
		}
		
		BestSolution.solverIntensity(0);
		BestSolution.quadraticPenalization();
		
		return BestSolution;
	}
	/**
	 * Este es exactamente el mismo LS que arriba, pero se le cambio el final para que evalue con el modelo penalizado a organos
	 * Este codigo es pura flojera,el uso de la evaluacion deberia estar ligado a un parametro de treatment
	 * Y quizas valga la pena generar una clase con funciones objetivos
	 * @param sol
	 * @param M
	 * @param mov
	 * @return
	 * @throws IOException
	 */
	public  TreatmentPlan LocalSearchPen(TreatmentPlan sol,DDM M, int mov) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		BestSolution.evaluationMethod=sol.evaluationMethod;
		BestSolution.evaluationFunction=sol.evaluationFunction;
		ArrayList<TreatmentPlan> neigborhood = new ArrayList<TreatmentPlan>();
		int[] apertureChanged =new int[2];
		int[][] stateAperture=getStateAperture();

		iterations=0;
	
		boolean condition=true;
		boolean update=false;
		int flag=0;
		long startTime = System.currentTimeMillis();
		BestSolution.getListLowIntensity(1);
		while(condition && flag<=1 ) {

			iterations=iterations+1;
			update=false;
			//---------------Se genera la vecindad
			neigborhood.clear();


			switch(mov) {
			
			case 1:
				movAperture2(neigborhood,BestSolution,M);
				//apertureChanged=movLeafsAperture(neigborhood,BestSolution,M);
				break;
				
			case 2:
				movLeafsRandom(neigborhood,BestSolution,M,24);
				break;
				
			case 3:
				movLeafsAperture(neigborhood,BestSolution,M);
				movLeafsRandom(neigborhood,BestSolution,M,24);
			
			default:
				apertureChanged=movLeafsAperture(neigborhood,BestSolution,M);
				break;
			
			}

			if(neigborhood.size()>0 ) {
	
				actualSolution.updateSolDAO(bestNeigborPoolPen(neigborhood,M));
				//Si al optimizar la apertura producio mejora hace esto
				if(((actualSolution.singleObjectiveValue-BestSolution.singleObjectiveValue) < (-0.01))){
					
					BestSolution.updateSolDAO(actualSolution);
					//if(mov==1) {
					//	BestSolution.deleteLowIntensityApertureNotChanged(apertureChanged);
					//	stateAperture[apertureChanged[0]][apertureChanged[1]]+=1;
					//}
					update=true;

				}else {
					
				//En el caso contrario se actualiza la solucion con el solver con delta
					if(!update) {
						condition=false;
						actualSolution.updateSolDAO(BestSolution);
						actualSolution.solverIntensityPen(1);

						//si al usar el solver la solucion mejora se vuelve a optimizar las aperturas
						if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
							BestSolution.updateSolDAO(actualSolution);
							//if(mov==1)BestSolution.getListLowIntensity(1);
							update=true;
							condition=true;
						}
						//sino se termina el algoritm
						else {
							//System.out.println(flag);
							flag+=1;
							condition=true;
						}
					}					
				}	
			}
			//System.out.println(iterations+" "+BestSolution.singleObjectiveValue);
		}
		
		BestSolution.solverIntensityPen(0);
		BestSolution.quadraticPenalizationOARstPTV();
		
		return BestSolution;
	}
	
	public  TreatmentPlan reducedVariableNeighbourhoodSearch(TreatmentPlan sol,DDM M) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		ArrayList<TreatmentPlan> neigborhood = new ArrayList<TreatmentPlan>();
		int[] apertureChanged =new int[2];
		int[][] stateAperture=getStateAperture();

		iterations=0;
	
		boolean condition=true;
		boolean update=false;
		int flag=0;

		BestSolution.getListLowIntensity(1);
		while(condition && flag<=1) {
			iterations=iterations+1;
			update=false;
			//---------------Se genera la vecindad
			neigborhood.clear();


			movAperture2(neigborhood,BestSolution,M);
			//apertureChanged=movLeafsAperture(neigborhood,BestSolution,M);
			movLeafsRandom(neigborhood,BestSolution,M,24);

			

			if(neigborhood.size()>0 ) {
	
				actualSolution.updateSolDAO(bestNeigborPool(neigborhood,M));
				//Si al optimizar la apertura producio mejora hace esto
				if((actualSolution.singleObjectiveValue-BestSolution.singleObjectiveValue) < (-0.01) ){
					
					BestSolution.updateSolDAO(actualSolution);
					update=true;
					
					
				}else {
					
				//En el caso contrario se actualiza la solucion con el solver con delta
					if(!update) {
						condition=false;
						actualSolution.updateSolDAO(BestSolution);
						actualSolution.solverIntensity(1);
				
						
						
						
						//si al usar el solver la solucion mejora se vuelve a optimizar las aperturas
						if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
							BestSolution.updateSolDAO(actualSolution);
							update=true;
							condition=true;
						}
						else {
							flag+=1;
							condition=true;
						}
					}					
				}	
			}
		}
		
		BestSolution.solverIntensity(0);
		BestSolution.quadraticPenalization();
		
		return BestSolution;
	}
	
	
	
	public  TreatmentPlan VariableNeigborhoodDescent1(TreatmentPlan sol,DDM M) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		ArrayList<TreatmentPlan> neigborhood = new ArrayList<TreatmentPlan>();
		int[] apertureChanged =new int[2];
		int[][] stateAperture=getStateAperture();

		iterations=0;
		int mov=1;
		boolean condition=true;
		boolean update=false;

		BestSolution.getListLowIntensity(1);
		while(condition) {
			iterations=iterations+1;
			update=false;
			//---------------Se genera la vecindad
			neigborhood.clear();

			
			switch(mov) {
			
			case 1:
				movAperture2(neigborhood,BestSolution,M);
				break;
				
			case 2:
				movLeafsRandom(neigborhood,BestSolution,M,24);
				break;
			
			default:
				apertureChanged=movLeafsAperture(neigborhood,BestSolution,M);
				break;
			
			}
			

			if(neigborhood.size()>0 ) {
	
				actualSolution.updateSolDAO(bestNeigborPool(neigborhood,M));
				//Si al optimizar la apertura producio mejora hace esto
				if((actualSolution.singleObjectiveValue-BestSolution.singleObjectiveValue) < (-0.01) ){
					
					BestSolution.updateSolDAO(actualSolution);
					mov=1;
					update=true;
					
					
				}else {
					
				//En el caso contrario se actualiza la solucion con el solver con delta
					if(!update) {
						condition=false;
						actualSolution.updateSolDAO(BestSolution);
						actualSolution.solverIntensity();
				
						
						
						
						//si al usar el solver la solucion mejora se vuelve a optimizar las aperturas
						if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
							BestSolution.updateSolDAO(actualSolution);
							update=true;
							condition=true;
						}
						//si al usar el solver la solucion no mejora se cambia el movimiento
						else {
							mov+=1;
							condition=true;
							if (mov>2) {
								condition=false;
							}
						}
					}					
				}	
			}
		}
		
		BestSolution.solverIntensity(0);
		
		return BestSolution;
	}
	public  TreatmentPlan VariableNeigborhoodDescent2(TreatmentPlan sol,DDM M) throws IOException {
		TreatmentPlan BestSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		TreatmentPlan actualSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		BestSolution.updateSolDAO(sol);
		ArrayList<TreatmentPlan> neigborhood = new ArrayList<TreatmentPlan>();
		int[] apertureChanged =new int[2];
		int[][] stateAperture=getStateAperture();

		iterations=0;
		int mov=1;
		boolean condition=true;
		boolean update=false;
		int flag=0;
		int interarFlagInformation=0;
		BestSolution.getListLowIntensity(1);
		while(condition) {
			iterations=iterations+1;
			update=false;
			//---------------Se genera la vecindad
			neigborhood.clear();

			
			switch(mov) {
			
			case 2:
				movAperture2(neigborhood,BestSolution,M);
				break;
				
			case 1:
				movLeafsRandom(neigborhood,BestSolution,M,24);
				break;
			
			default:
				apertureChanged=movLeafsAperture(neigborhood,BestSolution,M);
				break;
			
			}
			

			if(neigborhood.size()>0 ) {
	
				actualSolution.updateSolDAO(bestNeigborPool(neigborhood,M));
				//Si al optimizar la apertura producio mejora hace esto
				if((actualSolution.singleObjectiveValue-BestSolution.singleObjectiveValue) < (-0.01) ){
					
					BestSolution.updateSolDAO(actualSolution);

					update=true;
					mov=1;
					
					
				}else {
					
				//En el caso contrario se actualiza la solucion con el solver con delta
					if(!update) {
						condition=false;
						actualSolution.updateSolDAO(BestSolution);
						actualSolution.solverIntensity();
				
						
						
						
						//si al usar el solver la solucion mejora se vuelve a optimizar las aperturas
						if(actualSolution.singleObjectiveValue < BestSolution.singleObjectiveValue){
							BestSolution.updateSolDAO(actualSolution);
							update=true;
							condition=true;
						}
						else {
							mov+=1;
							if (mov>2) {
								condition=false;
							}
						}
					}					
				}	
			}
		}
		
		BestSolution.solverIntensity(0);
		
		return BestSolution;
	}
	
	private void generateNeigborhood(ArrayList<TreatmentPlan> neigborhood, int k, TreatmentPlan BestSolution, DDM m2) {
		// TODO Auto-generated method stub
		neigborhood.clear();
		
		if(k==1) {
			movLeafsAperture(neigborhood,BestSolution,M);
		}else {
			
			movLeafRandom(neigborhood,BestSolution,M,24);
		}
	}

	public void fixSolution(TreatmentPlan solution){
		int numberOfApertureRemove=0;
		for(int i=0;i<numAngles;i++){
			
			//Eliminar aperturas con malas intensidades
			numberOfApertureRemove=solution.removeBadAperturePerBeam(i, 1);
			//agregar nuebas aperturas//
			generateRandomAperture(solution,i,numberOfApertureRemove);
			
		}
		solverIntensity(solution);
		
	}
	
	public TreatmentPlan fixIntensity(TreatmentPlan solution) {
		TreatmentPlan newSolution=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		newSolution.updateSolDAO(solution);
		newSolution.setBadIntensityPerAperture(1, 1);
		newSolution.evaluateSolution();
		return newSolution;
	}
	public  void generateRandomAperture(TreatmentPlan sol,int anglePosition, int newApertures){
		int [][] aper_shape, new_aper;

		Random r = new Random();
		int random_intensity,range;
		

		aper_shape = initial_aperture_shape.get(anglePosition); //initial shape for angle
		for(int z=0;z<newApertures;z++){
			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			
			
			//Solución inicial completamente Random
			for(int j=0; j<aper_shape.length; j++){

				
				range=aper_shape[j][1]-aper_shape[j][0];
				new_aper[j][0] = ((int) r.nextInt(range+1))+ aper_shape[j][0];
				range=aper_shape[j][1]-new_aper[j][0];
				new_aper[j][1] = (int) r.nextInt(range+1) + new_aper[j][0];
			}
			
			Aperture aper = new Aperture(new_aper, random_intensity);
			sol.apertures.get(anglePosition).add(aper);
			
		}

		
		
	}
	public void solverIntensity(TreatmentPlan solution) {
		
		 solution.aperturesToIntensityMatrixPerStationwithoutIntensty();
		 solution.ApertureMatrixStationToIntensityMatrixPlan();
		

		    //time=System.currentTimeMillis();
		    
	    Gurobi_Solver newModel;
		try {
			newModel = new Gurobi_Solver(solution,M,o,beamletAngles,dd,solution.weights);
			solution.setIntensity(newModel.newIntensity);
		    solution.setOF(newModel.objVal);
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
	}
	
	public  TreatmentPlan bestNeigborPool(ArrayList<TreatmentPlan> neighborhood,DDM M) {
		 
		TreatmentPlan best=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		best.singleObjectiveValue=9999999;

		for(TreatmentPlan actual:neighborhood) {	
			actual.aperturesToIntensityMatrixPerStation();
			actual.intensityMatrixStationToIntensityMatrixPlan();
			evaluateSolution(actual, M, o,actual.weights);
			if(actual.singleObjectiveValue<best.singleObjectiveValue) {
				best.updateSolDAO(actual);
				//best.aperturesToIntensityMatrixPerStation();
				//best.intensityMatrixStationToIntensityMatrixPlan();
				//evaluateSolution(best, M, o, best.weights);
				
			}
		}
		return best;
	}
	public  TreatmentPlan bestNeigborPoolPen(ArrayList<TreatmentPlan> neighborhood,DDM M) {
		 
		TreatmentPlan best=new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
		best.singleObjectiveValue=9999999;
		neighbourhoodEvaluationSolver(neighborhood);
		//neighbourhoodEvaluation(neighborhood);
		for(TreatmentPlan actual:neighborhood) {

			//actual.aperturesToIntensityMatrixPerStation();
			//actual.intensityMatrixStationToIntensityMatrixPlan();
			//evaluateSolutionORPenstPTV(actual, M, o);

			if(actual.singleObjectiveValue<best.singleObjectiveValue) {
				best.updateSolDAO(actual);
				
			}
		}
		return best;
	}
	//---------------- FUNCIONES GENERACION DE MOVIMIENTO
	
	/**
	 * Genera la vecindad de una apertura
	 * @param BestSolution
	 * @param angleChoose
	 * @param apertureChoose
	 * @param sizeAperture
	 * @return Lista de TreatmentNeighborhood
	 */
	public ArrayList<TreatmentNeighborhood> generateNeighborhoodAperture(TreatmentPlan BestSolution, int angleChoose, int apertureChoose){
		ArrayList<TreatmentNeighborhood> pool =new ArrayList<TreatmentNeighborhood>(); 
		int sizeAperture = BestSolution.apertures.elementAt(angleChoose).get(apertureChoose).aperture.length;
		
		for(int i=0;i<sizeAperture;i++){
			TreatmentNeighborhood neigborhoodnew1=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,true,0);
			neigborhoodnew1.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew1.solver=false;
			TreatmentNeighborhood neigborhoodnew2=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,true,1);
			neigborhoodnew2.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew2.solver=false;
			TreatmentNeighborhood neigborhoodnew3=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,false,0);
			neigborhoodnew3.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew3.solver=false;
			TreatmentNeighborhood neigborhoodnew4=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,false,1);
			neigborhoodnew4.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew4.solver=false;
			pool.add(neigborhoodnew1);
			pool.add(neigborhoodnew2);
			pool.add(neigborhoodnew3);
			pool.add(neigborhoodnew4);
			
		}
		return pool;
	}
	/**
	 * Genera la vecindad de una apertura
	 * @param BestSolution
	 * @param angleChoose
	 * @param apertureChoose
	 * @param sizeAperture
	 * @return Lista de TreatmentNeighborhood
	 */
	public ArrayList<TreatmentNeighborhood> generateNeighborhoodAperture(TreatmentPlan BestSolution, int angleChoose, int apertureChoose, int sizeAperture){
		ArrayList<TreatmentNeighborhood> pool =new ArrayList<TreatmentNeighborhood>(); 
		
		for(int i=0;i<sizeAperture;i++){
			TreatmentNeighborhood neigborhoodnew1=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,true,0);
			neigborhoodnew1.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew1.solver=false;
			TreatmentNeighborhood neigborhoodnew2=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,true,1);
			neigborhoodnew2.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew2.solver=false;
			TreatmentNeighborhood neigborhoodnew3=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,false,0);
			neigborhoodnew3.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew3.solver=false;
			TreatmentNeighborhood neigborhoodnew4=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,false,1);
			neigborhoodnew4.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew4.solver=false;
			pool.add(neigborhoodnew1);
			pool.add(neigborhoodnew2);
			pool.add(neigborhoodnew3);
			pool.add(neigborhoodnew4);
			
		}
		
		return pool;
	}
	/**
	 * Genera la vecindad de una apertura
	 * @param BestSolution
	 * @param angleChoose
	 * @param apertureChoose
	 * @param sizeAperture
	 * @return Lista de TreatmentNeighborhood
	 */
	public ArrayList<TreatmentNeighborhood> generateNeighborhoodAperture(TreatmentPlan BestSolution, int angleChoose, int apertureChoose, int sizeAperture, boolean solver){
		ArrayList<TreatmentNeighborhood> pool =new ArrayList<TreatmentNeighborhood>(); 
		
		for(int i=0;i<sizeAperture;i++){
			TreatmentNeighborhood neigborhoodnew1=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,true,0);
			neigborhoodnew1.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew1.solver=solver;
			TreatmentNeighborhood neigborhoodnew2=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,true,1);
			neigborhoodnew2.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew2.solver=solver;
			TreatmentNeighborhood neigborhoodnew3=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,false,0);
			neigborhoodnew3.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew3.solver=solver;
			TreatmentNeighborhood neigborhoodnew4=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,false,1);
			neigborhoodnew4.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew4.solver=solver;
			pool.add(neigborhoodnew1);
			pool.add(neigborhoodnew2);
			pool.add(neigborhoodnew3);
			pool.add(neigborhoodnew4);
			
		}
		
		return pool;
	}

	
	/**
	 * This fuction generate new aperture swapping two father apertures.
	 * @param BestSolution
	 * @param angleChoose
	 * @param i
	 * @param j
	 * @return
	 */
	private ArrayList<int[][]> CrossoverAperture(TreatmentPlan BestSolution,int angleChoose, int i, int j) {

		return BestSolution.singlePointCrossover(angleChoose, i, j);
	}

	public void combination() {
		
		if(apertureCombination==null) {
			int total=(num_apertures*(num_apertures-1))/2;
			int[][] comb=new int[2][total];
			
			int iter=0;
			for(int i=0;i<num_apertures;i++) {
	
				for(int j=i+1;j<num_apertures;j++) {
					comb[0][iter]=i;
					comb[1][iter]=j;
					iter++;
				}
			}
			
			apertureCombination=comb;
		}
		
		
		
	}

    
/**
 * Genera los vecinos de un conjunto de hojas seleccionadas al azar en toda la solucion
 * @param BestSolution solucion de partida
 * @param nNeighbor cantidad de hoja que se van a seleccionar
 * @return
 */

	public ArrayList<TreatmentNeighborhood> generateNeighborhoodLeafs(TreatmentPlan BestSolution, int nNeighbor){
		ArrayList<TreatmentNeighborhood> pool =new ArrayList<TreatmentNeighborhood>(); 
		Random r = new Random(System.currentTimeMillis());
		int numAperture=BestSolution.apertures.get(0).size();
		int TotalLeafs=0;
		int leafsAperture=0;
		int[]leafShape=new int[numAngles];
		
		//Se obtiene cuantas hojas hay en total por cada angulo
		for(int i =0;i<numAngles;++i) {
			leafShape[i]=initial_aperture_shape.get(i).length;
			leafsAperture=leafsAperture+leafShape[i];
		}
		
		TotalLeafs=leafsAperture*numAperture;//por ahora sirve, pero deberia ser numAperture
		int []listNumber=RandomNumbers(nNeighbor,TotalLeafs);//genera una lista
		
		for(int j=0;j<listNumber.length;j++){
			int[] position =getposition(listNumber[j],TotalLeafs, leafsAperture,leafShape);
			int angleChoose=position[0];
			int apertureChoose=position[1];
			int leaf=position[2];//hoja
			
			
			TreatmentNeighborhood neigborhoodnew1=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,leaf,true,0,angleChoose,apertureChoose);
			neigborhoodnew1.generateNeighborhood();
			neigborhoodnew1.solver=false;
			TreatmentNeighborhood neigborhoodnew2=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,leaf,true,1,angleChoose,apertureChoose);
			neigborhoodnew2.generateNeighborhood();
			neigborhoodnew2.solver=false;
			TreatmentNeighborhood neigborhoodnew3=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,leaf,false,0,angleChoose,apertureChoose);
			neigborhoodnew3.generateNeighborhood();
			neigborhoodnew3.solver=false;
			TreatmentNeighborhood neigborhoodnew4=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,leaf,false,1,angleChoose,apertureChoose);
			neigborhoodnew4.generateNeighborhood();
			neigborhoodnew4.solver=false;
			
		
			
			
			pool.add(neigborhoodnew1);
			pool.add(neigborhoodnew2);
			pool.add(neigborhoodnew3);
			pool.add(neigborhoodnew4);
			
		}
		
		return pool;
	}
	
	
	/**
	 * Genera los vecinos de un conjunto de hojas seleccionadas al azar en toda la solucion
	 * @param BestSolution solucion de partida
	 * @param nNeighbor cantidad de hoja que se van a seleccionar
	 * @return
	 */

		public ArrayList<TreatmentNeighborhood> generateNeighborhoodLeafRandom(TreatmentPlan BestSolution, int nNeighbor){
			ArrayList<TreatmentNeighborhood> pool =new ArrayList<TreatmentNeighborhood>(); 
			num_apertures=BestSolution.apertures.get(0).size();
			int[]leafShape=new int[numAngles];
			
			//Se obtiene cuantas hojas hay en total por cada angulo
			for(int i =0;i<numAngles;++i) {
				leafShape[i]=initial_aperture_shape.get(i).length;
				
			}

			HashMap<String, int[]>mapPosition=RandomNumbersMap(nNeighbor,leafShape);//genera una lista
			ArrayList< int[]> list = new ArrayList< int[]>(mapPosition.values());
			for(int j=0;j<list.size();j++){
				int[] position =list.get(j);
				//int angleChoose=position[0];
				//int apertureChoose=position[1];
				//int leaf=position[2];//hoja
				//int state=position[3]
				if(position[3]==1) {//Si es hoja derecha usa esta
					TreatmentNeighborhood neigborhoodnew1=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,position[2],true,0,position[0],position[1]);
					neigborhoodnew1.generateNeighborhood();
					neigborhoodnew1.solver=false;
					TreatmentNeighborhood neigborhoodnew2=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,position[2],true,1,position[0],position[1]);
					neigborhoodnew2.generateNeighborhood();
					neigborhoodnew2.solver=false;
					pool.add(neigborhoodnew1);
					pool.add(neigborhoodnew2);
				
				}
				
				else {//Si es hoja izquierda usa esta
					TreatmentNeighborhood neigborhoodnew3=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,position[2],false,0,position[0],position[1]);
					neigborhoodnew3.generateNeighborhood();
					neigborhoodnew3.solver=false;
					TreatmentNeighborhood neigborhoodnew4=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,position[2],false,1,position[0],position[1]);
					neigborhoodnew4.generateNeighborhood();
					neigborhoodnew4.solver=false;

					pool.add(neigborhoodnew3);
					pool.add(neigborhoodnew4);
				}
				
			}
			
			return pool;
		}
	
	
	/**
	 * (Funcion deprecada)Parte de una solucion inicial con intensidades optimzas y despues calcula al final las intensidades
	 * @param initial
	 * @param BestSolution
	 * @param M
	 */
	public  void generateNeighborhoodPoolMix(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M) {
		// creamos una copia de 
		long Vecindad=System.currentTimeMillis();
		

		Random r = new Random(System.currentTimeMillis());
		
		
		
		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();
		
		int angleSize=BestSolution.apertures.size()-1;
		int angleChoose= r.nextInt(angleSize+1);
		double delta=1;
		int apertureChoose= getAperture(BestSolution, angleChoose, delta);
		List<Aperture> listApertures = BestSolution.apertures.elementAt(angleChoose);
		int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();
		int size=0;
		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);
		TreatmentNeighborhood[] bestForLeaf=new TreatmentNeighborhood[sizeAperture];
		//System.out.println("Angulo:"+angleChoose+" Apertura:"+apertureChoose);
		ArrayList<TreatmentNeighborhood> pool =	generateNeighborhoodAperture(BestSolution, angleChoose, apertureChoose, sizeAperture);
		

		//System.out.println("CreacionVecindad: "+(System.currentTimeMillis()-Vecindad));
		//long iteracion=System.currentTimeMillis();
		//System.out.println(pool.size());
		for(int i=0;i<pool.size();i++) {
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		for(int i=0;i<pool.size();i++) {
			
			if(pool.get(i).isBetter==true) {
				finalBest.add(pool.get(i));
				
			}
		}
		
		size=finalBest.size();
		pool.clear();
		boolean notChange=false;
		//System.out.println("IteracionInicial: "+(System.currentTimeMillis()-iteracion));
		while(size>0 && notChange != true) {
			//iteracion=System.currentTimeMillis();
			//initial.clear();

			for(int i=0;i<pool.size();i++) {
				if(pool.get(i).isBetter==true) {
					finalBest.add(pool.get(i));
				}
			}
			for(int i=0;i<finalBest.size();i++) {
				if(bestForLeaf[finalBest.get(i).leaf]==null) {
					bestForLeaf[finalBest.get(i).leaf]=finalBest.get(i);
					//finalBest.get(i).id();
					//printIntensitiesAperture(finalBest.get(i).neigborhoodTreatment.get(0));
					
				}else {
					if(!finalBest.get(i).neigborhoodTreatment.isEmpty()) {
						if(bestForLeaf[finalBest.get(i).leaf].singleObjectiveValue>finalBest.get(i).neigborhoodTreatment.get(0).singleObjectiveValue) {
							//System.out.print("entre matriz");
							bestForLeaf[finalBest.get(i).leaf]=finalBest.get(i);
							//finalBest.get(i).id();
							//printIntensitiesAperture(finalBest.get(i).neigborhoodTreatment.get(0));
						}
					}
				}
				
			}
			pool.clear();
			

			 //System.out.println("Entre bucle "+finalBest.size());
			for(int i=0;i<finalBest.size();i++) {
				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
				finalBest.get(i).generateNeighborhood(angleChoose,apertureChoose);
				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
			}
			//System.out.println("antes de llamado"+finalBest.size());
			List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
			ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());

			for(int i=0;i<finalBest.size();i++) {
				ncalls.add(Executors.callable(finalBest.get(i)));
				
			}
			
			try {
				poolThread.invokeAll(ncalls);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			poolThread.shutdown();
			//System.out.println("sali pool bucle");
			notChange=true;
			for(int i=0;i<finalBest.size();i++) {
				if(finalBest.get(i).isBetter==true) {
					
					pool.add(finalBest.get(i));
					notChange=false;
				}
			}
			
			//System.out.println("Sali bucle "+finalBest.size());
			//System.out.println("Sali bucle "+pool.size());
			size=pool.size();
			finalBest.clear();
			//System.out.println("Iteracion: "+(System.currentTimeMillis()-iteracion));
		}
		
		//Se obtiene la mejor intensidad para los mejores resultados
		
		long Solver=System.currentTimeMillis();
		
		finalBest.clear();
		for(int i=0;i<bestForLeaf.length;i++) {
			if(bestForLeaf[i]!=null) {
				//bestForLeaf[i].solver=true;
				bestForLeaf[i].refresh=true;
				//System.out.println(bestForLeaf[i].singleObjectiveValue);
				finalBest.add(bestForLeaf[i]);
			}
		}

		finalBest.add(joinTreatmentNeighborhood(angleChoose,apertureChoose,BestSolution,bestForLeaf,M));
		
		List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
		ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());
		//System.out.println("Antes");
		for(int i=0;i<finalBest.size();i++) {
			finalBest.get(i).refresh=true;
			//finalBest.get(i).actualSolution.aperturesToIntensityMatrixPerStationwithoutIntensty();
			//finalBest.get(i).actualSolution.ApertureMatrixStationToIntensityMatrixPlan();
			//evaluateSolution(finalBest.get(i).actualSolution, M, o, w);
			//System.out.println(finalBest.get(i).singleObjectiveValue);
			ncalls.add(Executors.callable(finalBest.get(i)));
			
		}
		
		try {
			poolThread.invokeAll(ncalls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		poolThread.shutdown();
		//System.out.println("Despues");
		
		for(int i=0;i<finalBest.size();i++) {
			//System.out.println("Tratamiento["+i +"]:"+finalBest.get(i).actualSolution.singleObjectiveValue);
			//finalBest.get(i).actualSolution.aperturesToIntensityMatrixPerStationwithoutIntensty();
			//finalBest.get(i).actualSolution.ApertureMatrixStationToIntensityMatrixPlan();
			//evaluateSolution(finalBest.get(i).actualSolution, M, o, w);
			//System.out.println(finalBest.get(i).actualSolution.singleObjectiveValue);
			initial.add(finalBest.get(i).actualSolution);
			
		}
		//System.out.println("sali bucle "+initial.size());
		//devuelve los mejores

		
	}
	

	public ArrayList<TreatmentNeighborhood> generateTest(TreatmentPlan BestSolution, int angleChoose, int apertureChoose, int sizeAperture){
ArrayList<TreatmentNeighborhood> pool =new ArrayList<TreatmentNeighborhood>(); 
		
		for(int i=0;i<sizeAperture;i++){
			TreatmentNeighborhood neigborhoodnew1=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,true,0);
			neigborhoodnew1.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew1.solver=false;
			TreatmentNeighborhood neigborhoodnew2=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,true,1);
			neigborhoodnew2.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew2.solver=false;
			TreatmentNeighborhood neigborhoodnew3=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,false,0);
			neigborhoodnew3.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew3.solver=false;
			TreatmentNeighborhood neigborhoodnew4=new TreatmentNeighborhood(BestSolution,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,i,false,1);
			neigborhoodnew4.generateNeighborhood(angleChoose,apertureChoose);
			neigborhoodnew4.solver=false;
			pool.add(neigborhoodnew1);
			pool.add(neigborhoodnew2);
			pool.add(neigborhoodnew3);
			pool.add(neigborhoodnew4);
		}
		return pool;
	}
	public void movTest(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M) {

		Random r = new Random(System.currentTimeMillis());

		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();
		
		int angleSize=BestSolution.apertures.size()-1;
		int angleChoose= r.nextInt(angleSize+1);
		double delta=1;
		int apertureChoose= getAperture(BestSolution, angleChoose, delta);
		List<Aperture> listApertures = BestSolution.apertures.elementAt(angleChoose);
		int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();

		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);

		//System.out.println("Angulo:"+angleChoose+" Apertura:"+apertureChoose);
		ArrayList<TreatmentNeighborhood> pool =	generateTest(BestSolution, angleChoose, apertureChoose, sizeAperture);
		//ArrayList<TreatmentNeighborhood> pool =	generateTest(BestSolution);

		for(int i=0;i<pool.size();i++) {
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		for(int i=0;i<pool.size();i++) {
			
			if(pool.get(i).isBetter==true) {
				finalBest.add(pool.get(i));
				
			}
		}
		
		for(int i=0;i<finalBest.size();i++) {

			initial.add(finalBest.get(i).actualSolution);
			
		}
		
	}
	/**
	 * Esta funcion evlaua todas las soluciones de la vecindad de forma paralela
	 * @param initial
	 * @param BestSolution
	 * @param M
	 */
	public void neighbourhoodEvaluation(ArrayList<TreatmentPlan> neighbour) {


		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();

		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);


		for(int i=0;i<neighbour.size();i++) {
			
			calls.add(Executors.callable(neighbour.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		

		
	}
	
	/**
	 * Esta funcion evlaua todas las soluciones de la vecindad de forma paralela con el solver
	 * @param initial
	 * @param BestSolution
	 * @param M
	 */
	public void neighbourhoodEvaluationSolver(ArrayList<TreatmentPlan> neighbour) {


		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();

		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);

		int aux=neighbour.get(0).evaluationMethod;
		for(int i=0;i<neighbour.size();i++) {
			neighbour.get(i).evaluationMethod=3;
			calls.add(Executors.callable(neighbour.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		for(int i=0;i<neighbour.size();i++) {
			neighbour.get(i).evaluationMethod=aux;
			
		}
		

		
	}
	

	
	/**
	 * Genera una vecindad moviendo las hojas de una apertura en especifico
	 * @param initial
	 * @param BestSolution
	 * @param M
	 */
	public  void movIntensification(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M) {
		// creamos una copia de 
		//long Vecindad=System.currentTimeMillis();
		

		Random r = new Random(System.currentTimeMillis());
		
		
		
		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();
		
		int angleSize=BestSolution.apertures.size()-1;
		int angleChoose= r.nextInt(angleSize+1);
		double delta=1;
		int apertureChoose= getAperture(BestSolution, angleChoose, delta);
		List<Aperture> listApertures = BestSolution.apertures.elementAt(angleChoose);
		int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();
		int size=0;
		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);
		TreatmentNeighborhood[] bestForLeaf=new TreatmentNeighborhood[sizeAperture];
		//System.out.println("Angulo:"+angleChoose+" Apertura:"+apertureChoose);
		ArrayList<TreatmentNeighborhood> pool =	generateNeighborhoodAperture(BestSolution, angleChoose, apertureChoose, sizeAperture,true);
		

		//System.out.println("CreacionVecindad: "+(System.currentTimeMillis()-Vecindad));
		//long iteracion=System.currentTimeMillis();
		//System.out.println(pool.size());
		
		//Comienza primera iteracion
		for(int i=0;i<pool.size();i++) {
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		

		//calculamos la intensidad 
		for(int i=0;i<pool.size();i++) {
			
			if(pool.get(i).isBetter==true) {
				pool.get(i).solver=false;
				finalBest.add(pool.get(i));
				//System.out.println(pool.get(i).singleObjectiveValue);
			}
		}
		
		
		
		size=finalBest.size();
		finalBest.clear();
		//pool.clear();
		boolean notChange=false;
		//System.out.println("IteracionInicial: "+(System.currentTimeMillis()-iteracion));
		int iteration=1;
		while(size>0 && notChange != true) {
			iteration++;
			//iteracion=System.currentTimeMillis();
			//initial.clear();

			for(int i=0;i<pool.size();i++) {
				if(pool.get(i).isBetter==true) {
					finalBest.add(pool.get(i));
				}
			}
			for(int i=0;i<finalBest.size();i++) {
				if(bestForLeaf[finalBest.get(i).leaf]==null) {
					bestForLeaf[finalBest.get(i).leaf]=finalBest.get(i);
					//finalBest.get(i).id();
					//printIntensitiesAperture(finalBest.get(i).neigborhoodTreatment.get(0));
					
				}else {
					if(!finalBest.get(i).neigborhoodTreatment.isEmpty()) {
						if(bestForLeaf[finalBest.get(i).leaf].singleObjectiveValue>finalBest.get(i).neigborhoodTreatment.get(0).singleObjectiveValue) {
							//System.out.print("entre matriz");
							bestForLeaf[finalBest.get(i).leaf]=finalBest.get(i);
							//finalBest.get(i).id();
							//printIntensitiesAperture(finalBest.get(i).neigborhoodTreatment.get(0));
						}
					}
				}
				
			}
			pool.clear();
			

			 //System.out.println("Entre bucle "+finalBest.size());
			for(int i=0;i<finalBest.size();i++) {
				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
				finalBest.get(i).generateNeighborhood(angleChoose,apertureChoose);
				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
			}
			//System.out.println("antes de llamado"+finalBest.size());
			List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
			ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());

			for(int i=0;i<finalBest.size();i++) {
				ncalls.add(Executors.callable(finalBest.get(i)));
				
			}
			
			try {
				poolThread.invokeAll(ncalls);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			poolThread.shutdown();
			//System.out.println("sali pool bucle");
			notChange=true;
			for(int i=0;i<finalBest.size();i++) {
				if(finalBest.get(i).isBetter==true) {
					
					pool.add(finalBest.get(i));
					notChange=false;
				}
			}
			
			//System.out.println("Sali bucle "+finalBest.size());
			//System.out.println("Sali bucle "+pool.size());
			size=pool.size();
			finalBest.clear();
			//System.out.println("Iteracion: "+(System.currentTimeMillis()-iteracion));
		}
		
		//Se obtiene la mejor intensidad para los mejores resultados

		
		
		for(int i=0;i<bestForLeaf.length;i++) {
			if(bestForLeaf[i]!=null) {
				//bestForLeaf[i].solver=true;
				//System.out.println();
				bestForLeaf[i].refresh=true;
				//System.out.println(bestForLeaf[i].singleObjectiveValue);
				finalBest.add(bestForLeaf[i]);
			}
		}
	
		for(int i=0;i<finalBest.size();i++) {
			//System.out.println("Tratamiento["+i +"]:"+finalBest.get(i).actualSolution.singleObjectiveValue);
			//finalBest.get(i).actualSolution.aperturesToIntensityMatrixPerStationwithoutIntensty();
			//finalBest.get(i).actualSolution.ApertureMatrixStationToIntensityMatrixPlan();
			//evaluateSolution(finalBest.get(i).actualSolution, M, o, w);
			//System.out.println(finalBest.get(i).actualSolution.singleObjectiveValue);
			initial.add(finalBest.get(i).actualSolution);
			
		}
		//System.out.println(iteration);
		//devuelve los mejores

	}
	
	
	
	/**
	 * Genera una vecindad moviendo las hojas de una apertura en especifico
	 * @param initial
	 * @param BestSolution
	 * @param M
	 */
	public  void movAperture(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M) {
		// creamos una copia de 
		//long Vecindad=System.currentTimeMillis();
		

		Random r = new Random(System.currentTimeMillis());
		
		
		
		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();
		
		int angleSize=BestSolution.apertures.size()-1;
		int angleChoose= r.nextInt(angleSize+1);
		double delta=1;
		int apertureChoose= getAperture(BestSolution, angleChoose, delta);
		List<Aperture> listApertures = BestSolution.apertures.elementAt(angleChoose);
		int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();
		int size=0;
		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);
		TreatmentNeighborhood[] bestForLeaf=new TreatmentNeighborhood[sizeAperture];
		//System.out.println("Angulo:"+angleChoose+" Apertura:"+apertureChoose);
		ArrayList<TreatmentNeighborhood> pool =	generateNeighborhoodAperture(BestSolution, angleChoose, apertureChoose, sizeAperture);
		

		//System.out.println("CreacionVecindad: "+(System.currentTimeMillis()-Vecindad));
		//long iteracion=System.currentTimeMillis();
		//System.out.println(pool.size());
		for(int i=0;i<pool.size();i++) {
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		for(int i=0;i<pool.size();i++) {
			
			if(pool.get(i).isBetter==true) {
				finalBest.add(pool.get(i));
				
			}
		}
		
		size=finalBest.size();
		pool.clear();
		boolean notChange=false;
		//System.out.println("IteracionInicial: "+(System.currentTimeMillis()-iteracion));
		int iteration=1;
		while(size>0 && notChange != true) {
			iteration++;
			//iteracion=System.currentTimeMillis();
			//initial.clear();

			for(int i=0;i<pool.size();i++) {
				if(pool.get(i).isBetter==true) {
					finalBest.add(pool.get(i));
				}
			}
			for(int i=0;i<finalBest.size();i++) {
				if(bestForLeaf[finalBest.get(i).leaf]==null) {
					bestForLeaf[finalBest.get(i).leaf]=finalBest.get(i);
					//finalBest.get(i).id();
					//printIntensitiesAperture(finalBest.get(i).neigborhoodTreatment.get(0));
					
				}else {
					if(!finalBest.get(i).neigborhoodTreatment.isEmpty()) {
						if(bestForLeaf[finalBest.get(i).leaf].singleObjectiveValue>finalBest.get(i).neigborhoodTreatment.get(0).singleObjectiveValue) {
							//System.out.print("entre matriz");
							bestForLeaf[finalBest.get(i).leaf]=finalBest.get(i);
							//finalBest.get(i).id();
							//printIntensitiesAperture(finalBest.get(i).neigborhoodTreatment.get(0));
						}
					}
				}
				
			}
			pool.clear();
			

			 //System.out.println("Entre bucle "+finalBest.size());
			for(int i=0;i<finalBest.size();i++) {
				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
				finalBest.get(i).generateNeighborhood(angleChoose,apertureChoose);
				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
			}
			//System.out.println("antes de llamado"+finalBest.size());
			List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
			ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());

			for(int i=0;i<finalBest.size();i++) {
				ncalls.add(Executors.callable(finalBest.get(i)));
				
			}
			
			try {
				poolThread.invokeAll(ncalls);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			poolThread.shutdown();
			//System.out.println("sali pool bucle");
			notChange=true;
			for(int i=0;i<finalBest.size();i++) {
				if(finalBest.get(i).isBetter==true) {
					
					pool.add(finalBest.get(i));
					notChange=false;
				}
			}
			
			//System.out.println("Sali bucle "+finalBest.size());
			//System.out.println("Sali bucle "+pool.size());
			size=pool.size();
			finalBest.clear();
			//System.out.println("Iteracion: "+(System.currentTimeMillis()-iteracion));
		}
		
		//Se obtiene la mejor intensidad para los mejores resultados

		
		finalBest.clear();
		for(int i=0;i<bestForLeaf.length;i++) {
			if(bestForLeaf[i]!=null) {
				//bestForLeaf[i].solver=true;
				//System.out.println();
				bestForLeaf[i].refresh=true;
				//System.out.println(bestForLeaf[i].singleObjectiveValue);
				finalBest.add(bestForLeaf[i]);
			}
		}
		if(finalBest.size()!=0) {
			List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
			ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());

			for(int i=0;i<finalBest.size();i++) {
				ncalls.add(Executors.callable(finalBest.get(i)));
				
			}
			try {
				poolThread.invokeAll(ncalls);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			poolThread.shutdown();
		}
		//System.out.println();

	
		for(int i=0;i<finalBest.size();i++) {
			//System.out.println("Tratamiento["+i +"]:"+finalBest.get(i).actualSolution.singleObjectiveValue);
			//finalBest.get(i).actualSolution.aperturesToIntensityMatrixPerStationwithoutIntensty();
			//finalBest.get(i).actualSolution.ApertureMatrixStationToIntensityMatrixPlan();
			//evaluateSolution(finalBest.get(i).actualSolution, M, o, w);
			//System.out.println(finalBest.get(i).actualSolution.singleObjectiveValue);
			initial.add(finalBest.get(i).actualSolution);
			
		}
		//System.out.println(iteration);
		//devuelve los mejores

	}
	
	/**
	 * Genera una vecindad moviendo las hojas de una apertura en especifico (este es un codigo que reemplaza a movsleafsApertura ya que por alguna razon ese mov no sirve)
	 * @param initial
	 * @param BestSolution
	 * @param M
	 */
	public  void movAperture2(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M) {
		// creamos una copia de 
		//long Vecindad=System.currentTimeMillis();
		

		Random r = new Random(System.currentTimeMillis());
		
		
		
		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();
		ArrayList<TreatmentNeighborhood> changedTreatment=new ArrayList<TreatmentNeighborhood>();
		int angleSize=BestSolution.apertures.size()-1;
		int angleChoose= r.nextInt(angleSize+1);
		double delta=1;
		int apertureChoose= getAperture(BestSolution, angleChoose, delta);
		List<Aperture> listApertures = BestSolution.apertures.elementAt(angleChoose);
		int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();
		int size=0;
		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);

		
		ArrayList<TreatmentNeighborhood> pool =	generateNeighborhoodAperture(BestSolution, angleChoose, apertureChoose);
		

		for(int i=0;i<pool.size();i++) {
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		for(int i=0;i<pool.size();i++) {
			
			if(pool.get(i).isBetter==true) {
				finalBest.add(pool.get(i));
				
			}
		}
		changedTreatment.addAll(finalBest);
		size=finalBest.size();
		pool.clear();
		boolean notChange=false;
	

		while(size>0 && notChange != true) {


			for(int i=0;i<pool.size();i++) {
				if(pool.get(i).isBetter==true) {
					finalBest.add(pool.get(i));
				}
			}

			pool.clear();
			
			for(int i=0;i<finalBest.size();i++) {
				
				finalBest.get(i).generateNeighborhood(angleChoose,apertureChoose);
				
			}
			//System.out.println("antes de llamado"+finalBest.size());
			List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
			ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());

			for(int i=0;i<finalBest.size();i++) {
				ncalls.add(Executors.callable(finalBest.get(i)));
				
			}
			
			try {
				poolThread.invokeAll(ncalls);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			poolThread.shutdown();
			//System.out.println("sali pool bucle");
			notChange=true;
			for(int i=0;i<finalBest.size();i++) {
				if(finalBest.get(i).isBetter==true) {
					
					pool.add(finalBest.get(i));
					notChange=false;
				}
			}
			
			//System.out.println("Sali bucle "+finalBest.size());
			//System.out.println("Sali bucle "+pool.size());
			size=pool.size();
			if(notChange!=true)
				finalBest.clear();
			//System.out.println("Iteracion: "+(System.currentTimeMillis()-iteracion));
		}
		
		initial.add(joinTreatmentSolution(BestSolution,changedTreatment));
		

		for(int i=0;i<finalBest.size();i++) {

			initial.add(finalBest.get(i).actualSolution);
			
		}


	}	
	/**
	 * Mueve las hojas de una apertura solo una vez (Movimiento Base)
	 * @param initial
	 * @param BestSolution
	 * @param M
	 */
	public  void movOnlyLeaf(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M) {
		// creamos una copia de 
		long Vecindad=System.currentTimeMillis();
		

		Random r = new Random(System.currentTimeMillis());
		
		
		
		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();
		
		int angleSize=BestSolution.apertures.size()-1;
		int angleChoose= r.nextInt(angleSize+1);
		double delta=1;
		int apertureChoose= getAperture(BestSolution, angleChoose, delta);
		List<Aperture> listApertures = BestSolution.apertures.elementAt(angleChoose);
		int sizeAperture =listApertures.get(apertureChoose).aperture.length;
		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();
		int size=0;
		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);
		TreatmentNeighborhood[] bestForLeaf=new TreatmentNeighborhood[sizeAperture];
		//System.out.println("Angulo:"+angleChoose+" Apertura:"+apertureChoose);
		ArrayList<TreatmentNeighborhood> pool =	generateNeighborhoodAperture(BestSolution, angleChoose, apertureChoose, sizeAperture);
		

		//System.out.println("CreacionVecindad: "+(System.currentTimeMillis()-Vecindad));
		//long iteracion=System.currentTimeMillis();
		//System.out.println(pool.size());
		for(int i=0;i<pool.size();i++) {
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		for(int i=0;i<pool.size();i++) {
			
			if(pool.get(i).isBetter==true) {
				finalBest.add(pool.get(i));
				
			}
		}
		
		size=finalBest.size();
		pool.clear();
		boolean notChange=false;
		//System.out.println("IteracionInicial: "+(System.currentTimeMillis()-iteracion));
	
		
		// Aqui se actualizan las intensidaddes
	
		for(int i=0;i<finalBest.size();i++) {
			
			finalBest.get(i).refresh=true;
		}
		if(finalBest.size()!=0) {
			List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
			ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());

			for(int i=0;i<finalBest.size();i++) {
				ncalls.add(Executors.callable(finalBest.get(i)));
				
			}
			try {
				poolThread.invokeAll(ncalls);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			poolThread.shutdown();
		}
		for(int i=0;i<finalBest.size();i++) {
			//System.out.println("Tratamiento["+i +"]:"+finalBest.get(i).actualSolution.singleObjectiveValue);
			//finalBest.get(i).actualSolution.aperturesToIntensityMatrixPerStationwithoutIntensty();
			//finalBest.get(i).actualSolution.ApertureMatrixStationToIntensityMatrixPlan();
			//evaluateSolution(finalBest.get(i).actualSolution, M, o, w);
			//System.out.println(finalBest.get(i).actualSolution.singleObjectiveValue);
			initial.add(finalBest.get(i).actualSolution);
			
		}
		//System.out.println(iteration);
		//devuelve los mejores

	}
	
	/**
	 * 
	 * @param initial
	 * @param BestSolution
	 * @param M
	 * @param nMovs cantidad de
	 */
	public  void movLeafRandom(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M,int nNeighbor) {
		// creamos una copia de 

		

		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();
		
		
		double delta=1;

		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();

		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);

			//aqui no estaba el error
		ArrayList<TreatmentNeighborhood> pool =	generateNeighborhoodLeafs(BestSolution,nNeighbor);
		

		//System.out.println("CreacionVecindad: "+(System.currentTimeMillis()-Vecindad));
		//long iteracion=System.currentTimeMillis();
		//System.out.println(pool.size());
		//poolCheck(pool);
		for(int i=0;i<pool.size();i++) {
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		for(int i=0;i<pool.size();i++) {
			
			if(pool.get(i).isBetter==true) {
				pool.get(i).solver=false;
				finalBest.add(pool.get(i));
				
			}
		}
		//Aqui termina la primera  generacion
		
		
		int size=finalBest.size();
		pool.clear();
		boolean notChange=false;
		int index=0;

		

		//System.out.println("IteracionInicial: "+(System.currentTimeMillis()-iteracion));
		while(size>0 && notChange != true) {
			index++;
			//iteracion=System.currentTimeMillis();
			//initial.clear();
			
			//poolCheck(pool);
			
			int poolsize=pool.size();
			for(int i=0;i<pool.size();i++) {
				if(pool.get(i).isBetter==true) {
					
					finalBest.add(pool.get(i));
				}
			}
			//poolCheck(finalBest);
			pool.clear();
			

			 //System.out.println("Entre bucle "+finalBest.size());
			
			
			/// que deberia hacer esta parte, segun yo deberia generar los nuevos vecinos
			for(int i=0;i<finalBest.size();i++) {

				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
				finalBest.get(i).generateNeighborhood();
				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
			}
			//System.out.println("antes de llamado"+finalBest.size());
			List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
			ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());

			for(int i=0;i<finalBest.size();i++) {
				ncalls.add(Executors.callable(finalBest.get(i)));
				
			}
			
			try {
				poolThread.invokeAll(ncalls);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			poolThread.shutdown();
			//System.out.println("sali pool bucle");
			notChange=true;
			for(int i=0;i<finalBest.size();i++) {
				if(finalBest.get(i).isBetter==true) {
					
					pool.add(finalBest.get(i));
					//poolCheck(pool);
					notChange=false;
				}
			}
			
			
			//System.out.println("Sali bucle "+finalBest.size());
			//System.out.println("Sali bucle "+pool.size());
			size=pool.size();
			if(notChange!=true)
				finalBest.clear();
			
			
			
			//System.out.println("Iteracion: "+(System.currentTimeMillis()-iteracion));
		}
		

	
		for(int i=0;i<finalBest.size();i++) {
			//System.out.println("Tratamiento["+i +"]:"+finalBest.get(i).actualSolution.singleObjectiveValue);
			//finalBest.get(i).actualSolution.aperturesToIntensityMatrixPerStationwithoutIntensty();
			//finalBest.get(i).actualSolution.ApertureMatrixStationToIntensityMatrixPlan();
			//evaluateSolution(finalBest.get(i).actualSolution, M, o, w);
			//System.out.println(finalBest.get(i).actualSolution.singleObjectiveValue);
			initial.add(finalBest.get(i).actualSolution);
			
		}
		//System.out.println("sali bucle "+initial.size());
		//devuelve los mejores
		
		
		


	}
	
	
	/** Esta es la nueva version la cual considera el join y no trabaja con las hojas en pares
	 * 
	 * @param initial
	 * @param BestSolution
	 * @param M
	 * @param nMovs cantidad de
	 */
	public int[] movLeafsAperture(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M) {
		// creamos una copia de 

		
		ArrayList<TreatmentNeighborhood> changedTreatment=new ArrayList<TreatmentNeighborhood>();
		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();

		int[] apertureChoose= BestSolution.getApertureRandom();
	
		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();
		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);

		//System.out.println("Angulo:"+angleChoose+" Apertura:"+apertureChoose);
		ArrayList<TreatmentNeighborhood> pool =	generateNeighborhoodAperture(BestSolution, apertureChoose[0], apertureChoose[1]);
		

		for(int i=0;i<pool.size();i++) {
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		//Si el movimiento genero una mejora este se guarda
		for(int i=0;i<pool.size();i++) {
			
			if(pool.get(i).isBetter==true) {
				pool.get(i).solver=false;
				finalBest.add(pool.get(i));
				
			}
		}
		changedTreatment.addAll(finalBest);
		//Aqui termina la primera  generacion
		
		
		int size=finalBest.size();
		pool.clear();
		boolean notChange=false;
		int index=0;


		while(size>0 && notChange != true) {
			index++;
			//iteracion=System.currentTimeMillis();
			//initial.clear();
			
			//poolCheck(pool);
			
			int poolsize=pool.size();
			for(int i=0;i<pool.size();i++) {
				if(pool.get(i).isBetter==true) {
					
					finalBest.add(pool.get(i));
				}
			}
			//poolCheck(finalBest);
			pool.clear();
			

			 //System.out.println("Entre bucle "+finalBest.size());
			
			
			/// que deberia hacer esta parte, segun yo deberia generar los nuevos vecinos
			for(int i=0;i<finalBest.size();i++) {

				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
				finalBest.get(i).generateNeighborhood();
				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
			}
			//System.out.println("antes de llamado"+finalBest.size());
			List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
			ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());

			for(int i=0;i<finalBest.size();i++) {
				ncalls.add(Executors.callable(finalBest.get(i)));
				
			}
			
			try {
				poolThread.invokeAll(ncalls);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			poolThread.shutdown();
			//System.out.println("sali pool bucle");
			notChange=true;
			for(int i=0;i<finalBest.size();i++) {
				if(finalBest.get(i).isBetter==true) {
					
					pool.add(finalBest.get(i));
					//poolCheck(pool);
					notChange=false;
				}
			}
			
			
			//System.out.println("Sali bucle "+finalBest.size());
			//System.out.println("Sali bucle "+pool.size());
			size=pool.size();
			if(notChange!=true)
				finalBest.clear();
			
			
			
			//System.out.println("Iteracion: "+(System.currentTimeMillis()-iteracion));
		}
		

		initial.add(joinTreatmentSolution(BestSolution,changedTreatment));
		

		for(int i=0;i<finalBest.size();i++) {
			//System.out.println("Tratamiento["+i +"]:"+finalBest.get(i).actualSolution.singleObjectiveValue);
			//finalBest.get(i).actualSolution.aperturesToIntensityMatrixPerStationwithoutIntensty();
			//finalBest.get(i).actualSolution.ApertureMatrixStationToIntensityMatrixPlan();
			//evaluateSolution(finalBest.get(i).actualSolution, M, o, w);
			//System.out.println(finalBest.get(i).actualSolution.singleObjectiveValue);
			initial.add(finalBest.get(i).actualSolution);
			
		}
		return apertureChoose;
	}
	
	/** Esta es la nueva version la cual considera el join y no trabaja con las hojas en pares
	 * 
	 * @param initial
	 * @param BestSolution
	 * @param M
	 * @param nMovs cantidad de
	 */
	public  void movLeafsRandom(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M,int nNeighbor) {
		// creamos una copia de 


		ArrayList<TreatmentNeighborhood> changedTreatment=new ArrayList<TreatmentNeighborhood>();
		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();
		
		
		double delta=1;

		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();

		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);

			//aqui no estaba el error
		ArrayList<TreatmentNeighborhood> pool =	generateNeighborhoodLeafRandom(BestSolution,nNeighbor);
		

		//poolCheck(pool);
		for(int i=0;i<pool.size();i++) {
		//for(int i=0;i<1;i++) {
			pool.get(i).setPriority(Thread.MAX_PRIORITY);
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		//Si el movimiento genero una mejora este se guarda
		for(int i=0;i<pool.size();i++) {
			
			initial.add(pool.get(i).bestSolution);
				
			
		}

		
		


	}
	
	/** Esta es la nueva version la cual considera el join y no trabaja con las hojas en pares
	 * Esta revision deberia permitir quedar a las hojas que no siguieron mejorando, pero que aun asi obtienen mejora comparado con su original
	 * @param initial
	 * @param BestSolution
	 * @param M
	 * @param nMovs cantidad de
	 */
	public  void movLeafsRandomRevision(ArrayList<TreatmentPlan> initial, TreatmentPlan BestSolution, DDM M,int nNeighbor) {
		// creamos una copia de 

		
		ArrayList<TreatmentNeighborhood> changedTreatment=new ArrayList<TreatmentNeighborhood>();
		ArrayList<TreatmentNeighborhood> finalBest=new ArrayList<TreatmentNeighborhood>();
		

		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();

		ExecutorService pool2 = Executors.newFixedThreadPool(nThreads);

		ArrayList<TreatmentNeighborhood> pool =	generateNeighborhoodLeafRandom(BestSolution,nNeighbor);
		

		//System.out.println("CreacionVecindad: "+(System.currentTimeMillis()-Vecindad));
		//long iteracion=System.currentTimeMillis();
		//System.out.println(pool.size());
		//poolCheck(pool);
		for(int i=0;i<pool.size();i++) {
			calls.add(Executors.callable(pool.get(i)));
			
		}

		try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		
		
		//Si el movimiento genero una mejora este se guarda
		for(int i=0;i<pool.size();i++) {
			
			if(pool.get(i).isBetter==true) {
				pool.get(i).solver=false;
				finalBest.add(pool.get(i));
				
			}
		}
		changedTreatment.addAll(finalBest);
		//Aqui termina la primera  generacion
		
		
		int size=finalBest.size();
		pool.clear();
		boolean notChange=false;
		int index=0;


		while(size>0 && notChange != true) {
			index++;
			//iteracion=System.currentTimeMillis();
			//initial.clear();
			
			//poolCheck(pool);
			
			int poolsize=pool.size();
			for(int i=0;i<pool.size();i++) {
				if(pool.get(i).isBetter==true) {
					
					finalBest.add(pool.get(i));
				}
			}
			//poolCheck(finalBest);
			pool.clear();
			

			 //System.out.println("Entre bucle "+finalBest.size());
			
			
			/// que deberia hacer esta parte, segun yo deberia generar los nuevos vecinos
			for(int i=0;i<finalBest.size();i++) {

				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
				finalBest.get(i).generateNeighborhood();
				//System.out.println(finalBest.get(i).neigborhoodTreatment.size());
			}
			//System.out.println("antes de llamado"+finalBest.size());
			List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
			ExecutorService poolThread = Executors.newFixedThreadPool(finalBest.size());

			for(int i=0;i<finalBest.size();i++) {
				ncalls.add(Executors.callable(finalBest.get(i)));
				
			}
			
			try {
				poolThread.invokeAll(ncalls);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			poolThread.shutdown();
			//System.out.println("sali pool bucle");
			notChange=true;
			for(int i=0;i<finalBest.size();i++) {
				if(finalBest.get(i).isBetter==true) {
					
					pool.add(finalBest.get(i));
					//poolCheck(pool);
					notChange=false;
				}
			}
			
			
			//System.out.println("Sali bucle "+finalBest.size());
			//System.out.println("Sali bucle "+pool.size());
			size=pool.size();
			if(notChange!=true)
				finalBest.clear();
			
			
			
			//System.out.println("Iteracion: "+(System.currentTimeMillis()-iteracion));
		}
		

		initial.add(joinTreatmentSolution(BestSolution,changedTreatment));
		
		//initial.get(0).solverIntensity();
		//initial.get(0).evaluateSolution();
		for(int i=0;i<finalBest.size();i++) {
			//System.out.println("Tratamiento["+i +"]:"+finalBest.get(i).actualSolution.singleObjectiveValue);
			//finalBest.get(i).actualSolution.aperturesToIntensityMatrixPerStationwithoutIntensty();
			//finalBest.get(i).actualSolution.ApertureMatrixStationToIntensityMatrixPlan();
			//evaluateSolution(finalBest.get(i).actualSolution, M, o, w);
			//System.out.println(finalBest.get(i).actualSolution.singleObjectiveValue);
			initial.add(finalBest.get(i).actualSolution);
			
		}
		//System.out.println("sali bucle "+initial.size());
		//devuelve los mejores
	
		
		


	}
	
	
	
	
	
	
	
	
	
	
	
	
	public void poolCheck(ArrayList<TreatmentNeighborhood> poolArray) {
		for(int i=0;i<poolArray.size();i++) {
			for(int j=(i+1);j<poolArray.size();j++) {
				if(poolArray.get(i).idstring().equals(poolArray.get(j).idstring())) {
					String a=poolArray.get(i).idstring();
					String b=poolArray.get(j).idstring();
					
					System.out.println("error 11");
				}
			}
		}
	}
	

	//--------------------Funciones de evaluacion y factibilidad
	public  double evaluateSolution(TreatmentPlan solution, DDM M, Organs[] o, double w[]){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		
		double[] solutionVector = solution.intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		//double[] solutionVector = new double[]{7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 17.0, 17.0, 17.0, 13.0, 13.0, 17.0, 17.0, 17.0, 17.0, 17.0, 13.0, 13.0, 13.0, 17.0, 17.0, 17.0, 17.0, 17.0, 13.0, 13.0, 13.0, 13.0, 17.0, 17.0, 17.0, 17.0, 17.0, 13.0, 13.0, 13.0, 13.0, 15.0, 15.0, 15.0, 15.0, 15.0, 11.0, 11.0, 11.0, 11.0, 15.0, 15.0, 15.0, 15.0, 15.0, 11.0, 11.0, 11.0, 11.0, 15.0, 15.0, 15.0, 15.0, 15.0, 11.0, 11.0, 11.0, 11.0, 15.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0};
		//double[] solutionVector = new double[]{10,20,20,14,20,19,11,10,20,5,6,7,16,19,18,14,9,20,10,4,5,14,17,16,11,10,20,12,4,5,13,17,15,9,10,20,13,5,4,13,18,12,8,11,19,14,7,4,14,17,9,10,12,20,9,4,17,18,8,14,10,20,20,20,5,14,13,20,20,19,11,18,20,5,10,7,10,15,20,15,3,15,6,3,12,17,20,9,2,16,4,3,9,20,20,6,3,19,2,4,1,5,17,17,0,0,17,0,4,0,0,17,20,1,1,19,1,0,1,0,13,0,1,19,0,4,0,7,11,4,0,20,15,3,0,0,0,9,10,6,20,0,1,3,8,10,19,7,20,6,3,6,8,9,20,5,9,20,8,1,4,5,6,17,6,7,13,8,1,2,3,2,10,5,8,16,10,2,3,3,1,5,4,4,20,9,7,9,3,4,3,5,5,3,2,4,7,4,6,2,7,4,0,0,20,9,20,20,20,20,12,1,20,20,20,20,19,18,9,10,17,20,18,15,14,18,17,7,13,16,19,14,15,14,20,20,6,17,15,20,14,16,12,16,19,5,18,19,20,20,19,13,14,19,3,19,20,20,20,20,12,7,20,4,17,19,20,6,0,17,1,7,18,2,0,0,3,4,0,0,0,0,0,0,0,2,0,0,11,17,19,10,11,15,0,10,10,8,17,18,10,13,14,11,2,8,14,20,20,13,13,10,12,0,4,14,19,14,11,8,8,12,0,17,19,9,8,7,7,11,20,19,4,5,6,6};
		//double[] solutionVector=new double[] {10.31,20.00,20.00,14.49,20.00,19.37,11.44,10.08,20.00,5.23,6.49,7.40,16.28,18.54,17.56,13.60,9.32,20.00,10.41,3.67,4.82,13.81,16.93,15.60,10.51,9.81,20.00,11.63,4.02,4.74,12.96,16.99,14.88,8.83,9.79,20.00,13.07,4.90,3.71,13.46,17.62,12.28,8.30,11.12,19.34,13.87,7.48,4.13,13.53,16.78,9.34,9.68,12.06,20.00,8.51,4.18,17.12,18.09,7.79,14.35,10.17,20.00,20.00,20.00,4.55,14.37,13.27,20.00,20.00,19.11,10.89,18.46,20.00,4.95,10.04,6.66,9.82,15.16,19.87,14.52,3.23,14.87,6.21,3.34,11.99,16.70,20.00,8.77,2.48,16.18,3.84,2.63,9.30,20.00,20.00,5.58,2.86,19.32,2.23,3.68,0.68,5.38,17.18,17.28,0.00,0.00,16.50,0.00,4.28,0.00,0.00,17.06,20.00,1.04,0.58,19.19,0.92,0.00,0.54,0.45,13.30,0.00,1.39,18.92,0.00,4.37,0.00,7.01,11.15,4.37,0.00,20.00,15.45,2.62,0.00,0.00,0.00,8.54,9.78,5.60,20.00,0.00,1.43,3.35,8.33,10.28,19.01,6.61,20.00,5.81,2.98,6.48,8.28,8.88,20.00,5.25,8.81,20.00,8.45,1.12,3.79,5.25,5.60,16.52,6.01,7.16,12.88,7.73,0.68,1.61,3.07,1.61,10.29,5.49,7.75,16.13,9.73,1.64,3.13,2.94,1.12,5.49,3.85,3.76,20.00,8.51,6.98,9.28,2.81,3.76,3.33,5.10,5.30,2.95,2.44,3.95,7.35,3.54,5.58,1.85,6.66,3.53,0.00,0.00,20.00,9.14,19.85,20.00,20.00,20.00,12.22,1.47,20.00,20.00,20.00,20.00,18.61,17.75,8.77,9.74,16.53,19.93,18.08,15.31,13.80,18.36,16.93,6.52,13.03,16.03,19.12,13.84,14.91,13.61,20.00,19.84,6.39,16.74,15.23,19.66,13.69,15.69,12.40,15.65,18.58,4.84,18.19,18.92,20.00,20.00,18.85,12.91,13.85,19.30,3.44,19.49,20.00,20.00,20.00,20.00,12.08,7.44,19.75,3.63,17.48,18.65,20.00,6.44,0.00,17.11,1.27,7.16,18.34,1.94,0.00,0.00,2.97,3.64,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.83,0.00,0.00,10.85,17.47,18.61,9.50,11.06,14.53,0.00,9.56,9.64,8.43,17.10,17.73,10.13,12.52,13.76,10.97,1.70,7.53,13.55,20.00,20.00,13.21,12.79,9.77,12.10,0.00,3.59,14.10,19.19,13.58,10.81,8.01,8.19,11.79,0.00,17.35,19.34,9.14,8.15,6.98,7.11,11.09,20.00,18.97,3.56,5.08,5.88,5.87};
		//double[] solutionVector=new double[] {10.3097,20,20,14.4913,20,19.3684,11.4378,10.075,20,5.2253,6.49239,7.39881,16.2819,18.5441,17.563,13.6014,9.32086,20,10.4118,3.67387,4.81576,13.8092,16.9314,15.6042,10.5102,9.80699,20,11.6295,4.02062,4.74281,12.9592,16.9919,14.8824,8.83135,9.78841,20,13.0725,4.89888,3.70851,13.4564,17.6188,12.2771,8.30026,11.1225,19.3446,13.8665,7.48044,4.1291,13.5259,16.7837,9.34001,9.68003,12.0576,20,8.5066,4.17597,17.1152,18.0933,7.79289,14.3488,10.1703,20,20,20,4.55037,14.3681,13.2736,20,20,19.1062,10.8932,18.4637,20,4.95187,10.0398,6.66412,9.82459,15.1614,19.8728,14.5203,3.22511,14.8658,6.20638,3.34013,11.9883,16.6987,19.9997,8.77295,2.47631,16.1838,3.8422,2.63254,9.29702,20,20,5.58084,2.86308,19.3228,2.23349,3.67659,0.677008,5.37564,17.1813,17.2778,5.34E-10,1.04E-09,16.5044,1.38E-10,4.27637,6.97E-10,1.66E-10,17.0639,20,1.04173,0.576239,19.1878,0.924726,2.55E-10,0.541289,0.450083,13.3042,4.80E-10,1.39492,18.915,2.08E-10,4.37446,5.83E-10,7.01266,11.153,4.37109,8.41E-10,20,15.4508,2.61777,6.29E-10,7.45E-10,4.03E-10,8.54297,9.77775,5.60274,20,6.33E-10,1.43363,3.35208,8.32635,10.2826,19.0077,6.61247,20,5.80767,2.97914,6.47571,8.28413,8.88179,20,5.2468,8.80924,20,8.44605,1.12493,3.78676,5.25264,5.59573,16.5185,6.01393,7.15919,12.8792,7.73481,0.680547,1.60642,3.06524,1.60948,10.2937,5.48907,7.74604,16.1253,9.73036,1.64011,3.12747,2.93724,1.11502,5.48732,3.8463,3.75988,20,8.51367,6.97555,9.27777,2.80921,3.75663,3.32548,5.10275,5.30177,2.95004,2.43942,3.94816,7.35041,3.54059,5.58316,1.85172,6.65664,3.52539,3.64E-09,3.36E-09,20,9.14021,19.854,20,20,20,12.2175,1.4684,20,20,20,20,18.6082,17.7476,8.7686,9.73554,16.5349,19.9298,18.0823,15.3061,13.8026,18.3586,16.931,6.52291,13.0256,16.0252,19.1165,13.8403,14.9063,13.6087,20,19.8418,6.39398,16.7429,15.234,19.6629,13.6932,15.6869,12.3977,15.6517,18.5804,4.84421,18.1889,18.918,20,20,18.8455,12.9068,13.8491,19.3014,3.44182,19.4928,20,20,20,20,12.0805,7.44039,19.7482,3.63102,17.4829,18.6455,20,6.43969,3.12E-09,17.1054,1.26997,7.15592,18.3381,1.9401,5.58E-10,2.14E-10,2.96576,3.64303,1.14E-09,1.55E-09,2.63E-09,9.97E-11,1.23E-10,1.79E-10,5.91E-11,1.82695,4.63E-10,5.50E-10,10.8482,17.4672,18.6106,9.50046,11.0649,14.5312,1.64E-10,9.56407,9.63549,8.42675,17.102,17.7268,10.1277,12.5151,13.7564,10.9661,1.69777,7.53269,13.5452,20,20,13.2125,12.7909,9.77117,12.097,5.81E-10,3.5894,14.1014,19.1897,13.5777,10.805,8.00704,8.19408,11.7871,5.21E-09,17.3491,19.3405,9.136,8.14758,6.98342,7.10701,11.0933,20,18.9689,3.55873,5.08056,5.87941,5.86829	};
		//double[] solutionVector=new double[] {2,16,36,2,4,36,16,16,0,4,12,12,12,0,0,0,16,32,0,0,8,0,0,12,28,36,4,8,8,8,8,2,24,12,0,0,0,12,12,8,24,32,32,32,2,8,2,2,2,8,8,16,2,24,0,24,36,36,24,28,2,12,2,32,48,28,28,0,24,24,18,18,36,32,54,54,22,36,36,18,18,18,18,44,0,26,44,58,3,3,48,26,14,14,32,36,36,36,36,44,44,0,14,4,58,44,58,58,58,3,18,0,0,0,0,32,44,44,26,18,18,4,4,44,62,48,26,4,4,5,26,14,54,58,72,4,62,8,56,8,0,16,24,48,4,5,32,32,16,16,26,16,16,16,16,32,46,1,1,1,42,42,72,46,46,0,38,62,7,7,62,46,0,1,0,0,0,0,0,22,22,0,4,58,8,8,8,46,7,7,48,5,5,0,4,62,8,26,24,4,34,4,1,64,8,52,34,22,16,24,24,24,52,3,0,0,28,4,56,56,4,3,56,56,28,28,28,28,1,18,0,0,36,36,52,64,54,1,0,0,16,16,16,56,44,44,52,28,1,0,0,0,36,36,36,28,38,56,0,0,34,62,46,1,22,1,0,1,54,46,18,8,8,3,98,42,78,52,3,26,26,26,46,56,78,58,78,32,32,52,52,72,52,56,56,26,0,22,22,48,26,26,3,2,0,36,36,56,56,56,56,76,62,42,0,0,2,0,3,3,56,42,42,42,22,22,42,52,78,58,2,2,42,42,78,76};
		
		
		double F = 0.0, pen,score;
		
		//System.out.println("Largo solución :"+solutionVector.length);
		//Recorremos los organos
		for(int i=0; i<o.length; i++){
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			pen = 0.0;
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			//Recorremos claves de voxel por organo para su evaluación
			while(keys.hasMoreElements()){
	            key = keys.nextElement();
	            beams = aux_index.get(key);
	            intensityVoxel2Beams = 0.0;
	            for(int b=0;b<beams.size(); b++){
	            	beam = beams.get(b);
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	intensityVoxel2Beams+= solutionVector[beam] * radiation;
	         
	            }
	            
	            /*
	            if(intensityVoxel2Beams < Zmin[i]){
            		pen += w[i] * (Math.pow(Zmin[i]-intensityVoxel2Beams, 2) );
            	}
   				if(intensityVoxel2Beams > Zmax[i] ){
   					pen += w[i] * (Math.pow(intensityVoxel2Beams-Zmax[i], 2) );
   				}
   				*/
            	if(i == 0){
            		pen += w[i]* Math.pow((dd[i] - intensityVoxel2Beams),2);
            	}else{
   					pen += w[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				}
            	
	        }
			//System.out.println(i+" "+Zmax[i]+" ");
			score=pen/aux_index.size();
			solution.scores[i]=score;
			F+=pen/aux_index.size();
			//System.out.println("organo: "+solution.scores[i]);
		}
		solution.singleObjectiveValue=F;
		//System.out.println("Solucion: "+F);
		return F;
	}
	
	/**
	 * Funcion de evaluacion que minimiza la penalizacion de los organos en riesgo, no tiene en consideracion la restriccion del ptv
	 * @param solution
	 * @param M
	 * @param o
	 * @param w
	 * @return
	 */
	public  double evaluateSolutionORPen(TreatmentPlan solution, DDM M, Organs[] o, double w[]){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		
		double[] solutionVector = solution.intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		double F = 0.0, pen,score;
		
		//System.out.println("Largo solución :"+solutionVector.length);
		//Recorremos los organos
		for(int i=0; i<o.length; i++){
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			pen = 0.0;
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			//Recorremos claves de voxel por organo para su evaluación
			while(keys.hasMoreElements()){
	            key = keys.nextElement();
	            beams = aux_index.get(key);
	            intensityVoxel2Beams = 0.0;
	            for(int b=0;b<beams.size(); b++){
	            	beam = beams.get(b);
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	intensityVoxel2Beams+= solutionVector[beam] * radiation;
	         
	            }
	            
	            /*
	            if(intensityVoxel2Beams < Zmin[i]){
            		pen += w[i] * (Math.pow(Zmin[i]-intensityVoxel2Beams, 2) );
            	}
   				if(intensityVoxel2Beams > Zmax[i] ){
   					pen += w[i] * (Math.pow(intensityVoxel2Beams-Zmax[i], 2) );
   				}
   				*/
            	if(i == 0){
            		
            	}else{
   					pen += w[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				}
            	
	        }
			//System.out.println(i+" "+Zmax[i]+" ");
			
			score=pen/aux_index.size();
			solution.scores[i]=score;
			F+=pen/aux_index.size();
		}
		solution.singleObjectiveValue=F;

		return F;
	}
	/**
	 * Funcion de evaluacion que minimiza la penalizacion de los organos en riesgo,considera la restriccion del PTV
	 * @param solution
	 * @param M
	 * @param o
	 * @param w
	 * @return
	 */
	public  double evaluateSolutionORPenstPTV(TreatmentPlan solution, DDM M, Organs[] o){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		int Factible=1;
		double worstValue=999999;
		
		double[] solutionVector = solution.intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		double F = 0.0, pen,score;
		
		//System.out.println("Largo solución :"+solutionVector.length);
		//Recorremos los organos
		for(int i=0; i<o.length; i++){
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			pen = 0.0;
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			//Recorremos claves de voxel por organo para su evaluación
			while(keys.hasMoreElements()){
	            key = keys.nextElement();
	            beams = aux_index.get(key);
	            intensityVoxel2Beams = 0.0;
	            for(int b=0;b<beams.size(); b++){
	            	beam = beams.get(b);
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	intensityVoxel2Beams+= solutionVector[beam] * radiation;
	         
	            }
	            
	            /*
	            if(intensityVoxel2Beams < Zmin[i]){
            		pen += w[i] * (Math.pow(Zmin[i]-intensityVoxel2Beams, 2) );
            	}
   				if(intensityVoxel2Beams > Zmax[i] ){
   					pen += w[i] * (Math.pow(intensityVoxel2Beams-Zmax[i], 2) );
   				}
   				*/
            	if(i == 0){
	            	if(intensityVoxel2Beams<dd[i]) {
	            		Factible=0;
	            		break;
	            	}
            		
            	}else{
   					pen += solution.weights[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				}
            	
	        }
			//System.out.println(i+" "+Zmax[i]+" ");
			
			score=pen/aux_index.size();
			solution.scores[i]=score;
			F+=pen/aux_index.size();
		}
		solution.singleObjectiveValue=F;

		if(Factible==1) {
			return F;
		}
		else {
			solution.singleObjectiveValue=worstValue;
			return worstValue;
		}
	}
	
	public  double deltaAproxEvaluateSolution(TreatmentPlan solution, DDM M, Organs[] o, double w[],int deltaBeamlet){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> beamlet_to_voxel_dao_ddm = M.beamlet_to_voxel_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		ArrayList<Integer> voxels;
		String value_index_key;
		Double radiation, intensityVoxel2BeamsActual,intensityVoxel2BeamsPast;
		Integer key, voxel;
		
		double actualBeamlet = solution.intensity[deltaBeamlet]; //es el beamlet al cual se quiere calcular el delta
		double pastBeamlet=1;
		double delta=0;
		double F = 0.0, penActual,penPast,score;

		for(int i=0; i<o.length; i++){
			aux_index = beamlet_to_voxel_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			penActual = 0.0;
			penPast=0.0;
			//System.out.println("Organo: "+o[i].name);
		
			//Recorremos claves de voxel por organo para su evaluación
		
            key = deltaBeamlet;
            voxels = aux_index.get(deltaBeamlet);
            intensityVoxel2BeamsActual = 0.0;
            intensityVoxel2BeamsPast = 0.0;
            //obtenemos que voxel son afectados por este beam
            for(int b=0;b<voxels.size(); b++){
            	voxel = voxels.get(b);
            	value_index_key = voxel+"-"+key;
            	radiation = aux_values.get(value_index_key);
            //por cada voxel le vamos a restar
            	intensityVoxel2BeamsActual= actualBeamlet * radiation;
            	intensityVoxel2BeamsPast= pastBeamlet * radiation;
            	
            	if(i == 0){
            		penActual += w[i] * Math.pow((dd[i] - intensityVoxel2BeamsActual),2);
            		penPast += w[i] * Math.pow((dd[i] - intensityVoxel2BeamsPast),2);
            	}else{
            		penActual += w[i] * Math.pow(Math.max((intensityVoxel2BeamsActual - dd[i]),0),2);
            		penPast += w[i] * Math.pow(Math.max((intensityVoxel2BeamsPast - dd[i]),0),2);
    			}
         
            }
            
        	
	
        	delta=penPast-penActual;
			//System.out.println(i+" "+Zmax[i]+" ");
			score=delta/aux_values.size();
			solution.scores[i]=solution.scores[i]+score;
			F+=delta/aux_values.size();
		}
		solution.singleObjectiveValue=solution.singleObjectiveValue-F;
		//System.out.println("Solucion: "+F);
		return F;
	}
	/**
	 * Esta funcion retorna verdadero si es que la solucion es factible, en caso contrario devuelve falso
	 * @param sol
	 * @return
	 */
	public  boolean isFeasible(TreatmentPlan sol){
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
						//System.out.println("angulo"+i+"apertura"+j);
						feasible=false;
						
						return feasible; 
						
					}
					if(matrix[x][0]>matrix[x][1]) {
						//System.out.println(limit_shape[x][0]+ " " + matrix[x][0]+ " "+limit_shape[x][1]+ " " + matrix[x][1]);
						//System.out.println("angulo"+i+"apertura"+j);
						feasible=false;
						return feasible; 
					}

					
				}	
			}
				
		}
		return feasible;
	}
	
	
	/**
	 * Movimiento de union
	 * Junta las mejor apertura por par de hojas de cada uno de los tratamientos
	 * No calcula el valor de evaluacion ya que, quien lo use, tiene que calcularlo aparte 
	 * @param angleChoose
	 * @param apertureChoose
	 * @param actual
	 * @param bestForLeaf
	 * @param M
	 * @return
	 */
	public TreatmentNeighborhood joinTreatmentNeighborhood(int angleChoose,int apertureChoose,TreatmentPlan actual,TreatmentNeighborhood[] bestForLeaf, DDM M) {
		int[][] list_shape = initial_aperture_shape.get(angleChoose);
		TreatmentNeighborhood neigborhoodnew1=new TreatmentNeighborhood(actual,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity,initial_aperture_shape,numAngles,numOrgans,M,o,dd,0,true,0);
		//System.out.println("entre a la union ");
		
		//long tiempoCrearApertura=System.currentTimeMillis();
		
		for(int i=0;i<bestForLeaf.length;i++) {
			if(bestForLeaf[i]!=null) {
				for(int index=0;index<bestForLeaf[i].called;index++) {
					neigborhoodnew1.actualSolution.apertures.get(angleChoose).get(apertureChoose).changeRow(i,bestForLeaf[i].direction, bestForLeaf[i].state,list_shape[i][0],list_shape[i][1]);
				}
			}
		}
		return neigborhoodnew1;
	}
	/**
	 * Movimiento de union General
	 * Junta las mejor apertura por par de hojas de cada uno de los tratamientos
	 * No calcula el valor de evaluacion ya que, quien lo use, tiene que calcularlo aparte 
	 * @param angleChoose
	 * @param apertureChoose
	 * @param actual
	 * @param bestForLeaf
	 * @param M
	 * @return
	 */
	public TreatmentPlan joinTreatmentSolution(TreatmentPlan actual,ArrayList<TreatmentNeighborhood> neighborList) {
		TreatmentPlan newSol = new TreatmentPlan(numAngles, numOrgans);
		newSol.updateSolDAO(actual);
		
		for(int i=0;i<neighborList.size();i++) {
			TreatmentNeighborhood actualTreatment =neighborList.get(i);
			
			int [] vector1=actualTreatment.getApertureChanged();
			int [] vector2=actualTreatment.getApertureChanged2();
			newSol.apertures.get(actualTreatment.angleChoose).get(actualTreatment.apertureChoose).changeRowAperture(actualTreatment.leaf,vector1);
			

		}
		newSol.evaluateSolution();
		return newSol;
	}
	
	public int[] doubletoint(double[]doubleArray) {
		 int[] intArray = new int[doubleArray.length];
		for (int i=0; i<intArray.length; ++i)
		    intArray[i] = (int) doubleArray[i];
		return intArray;
	}
	
	public void setInitialApertureShape(TreatmentPlan solution){
		int sizeX, sizeY, x, y;
		apertureShapes=new Vector<int[][]>();
		initial_aperture_shape=new Vector<int[][]>();
		int beamlets=0;
		boolean leftLimit, rightLimit;
		
		for(int i=0;i<numAngles;i++){
			sizeX  = solution.selAngles[i].maxAbsX;
			sizeY = solution.selAngles[i].maxAbsY;
			int [][] apertureShape = new int[sizeX][sizeY];
			for(x=0; x<sizeX; x++){
				for(y=0; y<sizeY; y++){
					apertureShape[x][y] = -1;
				}
			}
			
			double[][] beamletsCoord = solution.selAngles[i].beamletsCoord;
			for(int j=0;j<solution.selAngles[i].beamlets;j++){
				int xcoord=(int)beamletsCoord[j][1]-1;
				int ycoord=(int)beamletsCoord[j][2]-1;
				apertureShape[xcoord][ycoord] = 0;
				beamlets=beamlets+1;	
			}
			apertureShapes.add(apertureShape);
			
			int[][] apertureVectorShape = new int[sizeX][2];
			for(x=0; x<sizeX; x++){ 
				leftLimit = false;
				rightLimit = false;
				for(y=0; y<sizeY; y++){
					if(!leftLimit && apertureShape[x][y] == 0){
						apertureVectorShape[x][0] = y;
						leftLimit = true;
					}
					if(leftLimit && !rightLimit && apertureShape[x][y] == 0){
						apertureVectorShape[x][1] = y;
						leftLimit = true;
					}
				}
			}
			
			initial_aperture_shape.add(apertureVectorShape);
			
		}
	}
    public int[] getposition(int number,int TotalLeafs,int leafsAperture,int[]leafShape) {
    	int [] position=new int[3];//angle,aper,leaf
    	int truncated = Math.round((number/leafsAperture) - 0.5f);
    	int apert=0;
    	int angle=truncated;
    	int aux=number-(angle*leafsAperture);
    	int leaf=0;
    	int sum=0;
    	for(int i=0;i<num_apertures;++i) {
    		sum=sum+leafShape[i];
    		
    		if(aux<sum) {
    			leaf=(aux+leafShape[i])-sum;
    			apert=i;
    			break;
    		}
    	}
    	position[0]=apert;
    	position[1]=angle ;
    	position[2]=leaf;
    	
    	return position;
    }
	
	public  int getAperture(TreatmentPlan BestSolution,int angleChoose,double delta) {
		Random r = new Random();
		int rand=0;
		int apertureChoose=0;
		double limite=0 +delta;
		double intensidad;
		while(rand<5) {
			rand++;
			apertureChoose= r.nextInt(BestSolution.apertures.get(angleChoose).size());
			intensidad=BestSolution.apertures.get(angleChoose).get(apertureChoose).intensity;
			if(intensidad>limite) {
				rand=5;
			}
		}

		return apertureChoose;
	}

	
	
	public  int getLowIntensityAperture(TreatmentPlan BestSolution,double delta) {
		Random r = new Random();

		int [] listApertureLowIntensity=BestSolution.getApertureRandom();
		int apertureChoose= r.nextInt(listApertureLowIntensity.length);
		return apertureChoose;
	}
	
	/**
	 * Funcion que devuelve una lista de numeros random sin repetir
	 * @param n
	 * @param maxRange
	 * @return
	 */
    public  int[] RandomNumbers(int n, int maxRange) {  
    	 Random gen = new Random(); 
        assert n <= maxRange : "cannot get more unique numbers than the size of the range";  
          
        int[] result = new int[n];  
        Set<Integer> used = new HashSet<Integer>();  
          
        for (int i = 0; i < n; i++) {  
              
            int newRandom;  
            do {  
                newRandom = gen.nextInt(maxRange);  
            } while (used.contains(newRandom));  
            result[i] = newRandom;  
            used.add(newRandom);  
        }  
        return result;  
    } 
	/**
	 * Funcion que devuelve un mapa a de numeros random sin repetir
	 * @param n
	 * @param maxRange
	 * @return
	 */
    public  HashMap<String, int[]> RandomNumbersMap(int n, int[] leaf) {  
    	HashMap<String, int[]> map = new HashMap<String, int[]>();
    	Random r = new Random(); 

    	
    	while(map.size()<n) {
    		int[] position=new int[4];
    		position[0] = r.nextInt(numAngles);
    		position[1] =r.nextInt(num_apertures);
    		position[2] =r.nextInt(leaf[position[0]]);
    		position[3] =r.nextInt(2);
    		
    		String key=String.valueOf(position[0])+String.valueOf(position[1])+String.valueOf(position[2])+String.valueOf(position[3]);
    		map.put(key, position);
    		
    	}
    	
        return map;  
    } 
}
