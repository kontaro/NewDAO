package MultiObj;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import IMRT_DAO.DDM;
import IMRT_DAO.Organs;
import IMRT_DAO.TreatmentPlan;
import SingleObj.Algorithm;
import SingleObj.TreatmentNeighborhood;
import gurobi.GRBException;

public class MultiObjectiveAlgorithm extends Algorithm {
	 public MultiObjectiveAlgorithm(){
		 super();
		
	}
	 public MultiObjectiveAlgorithm(int[][] beamletSet,int[] initialBACs,DDM M, Organs[] o, double[] dd, int[] beamletAngles, int init_intensity, int max_intensity,
				int max_delta, int max_iter, int max_time, int seed, int step_intensity, int numAngles, int numOrgans){
		 super(beamletSet, initialBACs, M, o, dd, beamletAngles, init_intensity, max_intensity, max_delta,  max_iter,  max_time,  seed,  step_intensity,  numAngles, numOrgans);
		
	}
	 
	 
		public TreatmentPlan PLS(DDM M,TreatmentPlan solution) throws IOException, GRBException {
		       
			setInitialApertureShape(solution); //guardamos la figura
		    
		    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
		    solution.evaluateSolution();
		   
		    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
		    ArrayList<TreatmentPlan> front=standardParetoLocalSearch(solution);
		    return front.get(0);
		}
		
		public ArrayList<TreatmentPlan> PLS2(DDM M,TreatmentPlan solution) throws IOException, GRBException {
		       
			setInitialApertureShape(solution); //guardamos la figura
		    
		    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
		    solution.evaluateSolution();
		   
		    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
		    ArrayList<TreatmentPlan> front=standardParetoLocalSearch(solution);
		    return front;
		}
		public TreatmentPlan Test(DDM M,TreatmentPlan solution) throws IOException, GRBException {
		       
			setInitialApertureShape(solution); //guardamos la figura
		    
		    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
		    solution.evaluateSolution();
		   
		    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
		    ArrayList<TreatmentPlan> front= new ArrayList<TreatmentPlan>();
		    front.add(solution);
    
		    
			navigationPhase(front);
		    return solution;
		}
		public TreatmentPlan Test2(DDM M,TreatmentPlan solution) throws IOException, GRBException {
		    double[][] matriz = {{6.119506558,5.79159529,2.0079E-09,6.87127E-13,0.33484031},{6.66042171,5.357001527,3.429922091,8.55192E-13,1.587193838},{7.519234978,4.488360749,1.555064553,3.08849E-12,3.785037792},{12.53562674,12.91069231,0.958313018,1.66664E-12,1.424697012},{7.8029117,10.94813508,2.12498E-11,1.00387E-12,3.96136739}};
			setInitialApertureShape(solution); //guardamos la figura
		    
		    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
		    solution.setIntensity(matriz);
		    solution.evaluateSolution();
		   
		    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
		    ArrayList<TreatmentPlan> front= new ArrayList<TreatmentPlan>();
		    front.add(solution);
    
		    
			navigationPhase(front);
		    return solution;
		}
		
		
		public TreatmentPlan LSFront(DDM M,TreatmentPlan solution) throws IOException, GRBException {
		       
			setInitialApertureShape(solution); //guardamos la figura
		    
		    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
		    solution.evaluateSolution();
		   
		    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
		    //solution.updateSolDAO(LocalSearchWithSolver(solution,M));
		    //solution.updateSolDAO(LocalSearchTwoStep(solution,M));
		    
		    ArrayList<TreatmentPlan> list= weightCopy(solution,weightGeneratorTestMO());
		    //ArrayList<TreatmentPlan> list= new ArrayList<TreatmentPlan>();
		    list.add(solution);
		    for(TreatmentPlan node:list) {
		    	node.updateSolDAO(LocalSearch(node,M,2));
		    }
		    
		    frontPrint(list);
		    navigationPhase(list);
		    
		    System.out.println("__________________ End____________________");
		    //System.out.println("value:"+solution.singleObjectiveValue);
		   
		    return solution;
		}
		
		
		
		public TreatmentPlan WeightedSum(DDM M,TreatmentPlan solution) throws IOException, GRBException {
			System.out.println("__________________ Start ____________________");   
			setInitialApertureShape(solution); //guardamos la figura
		    
		    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
		    solution.evaluateSolution();
		   
		    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
		    //solution.updateSolDAO(LocalSearchWithSolver(solution,M));
		    //solution.updateSolDAO(LocalSearchTwoStep(solution,M));
		    
		    ArrayList<TreatmentPlan> list= weightCopy(solution,weightGeneratorTestPTV());
		    //ArrayList<TreatmentPlan> list= new ArrayList<TreatmentPlan>();
		    list.add(solution);
		    for(TreatmentPlan node:list) {
		    	node.updateSolDAO(LocalSearch(node,M,2));
		    }
		    
		    frontPrintALL(list);
		    
		    System.out.println("__________________ End ____________________");
		    //System.out.println("value:"+solution.singleObjectiveValue);
		   
		    return list.get(0);
		}
		public TreatmentPlan initialFront(DDM M,TreatmentPlan solution) throws IOException, GRBException {
		       
			setInitialApertureShape(solution); //guardamos la figura
		    
		    solution.beamletCoord();//esto es necesario para poder saber a que beemlet pertenece una cuadricula de la apertura
		    solution.evaluateSolution();
		   
		    System.out.print("Initial: "+solution.singleObjectiveValue+" ");
		    //solution.updateSolDAO(LocalSearchWithSolver(solution,M));
		    //solution.updateSolDAO(LocalSearchTwoStep(solution,M));
		    ArrayList<double[]> weights=new ArrayList<double[]>();
			double a[]={5,1,0};// +PTV// +Rectum
			weights.add(a);
		    ArrayList<TreatmentPlan> list= weightCopy(solution,weights);
			    
			
		    
		    frontPrintALL(list);
		    navigationPhase(list);
		    
		   
		    //System.out.println("value:"+solution.singleObjectiveValue);
		   
			

		    return solution;
		}
	 /**
	  * This adaptation of the pareto ANGEL algorithm only modify
	  * @param initialSolution
	  * @return
	  */
	 public ArrayList<TreatmentPlan> standardParetoLocalSearch(TreatmentPlan initialSolution) {
		 ArrayList<TreatmentPlan> front= new ArrayList<TreatmentPlan>();
		 ArrayList<TreatmentPlan> notVisited= new ArrayList<TreatmentPlan>();
		 ArrayList<TreatmentPlan> visited= new ArrayList<TreatmentPlan>();
		 ArrayList<TreatmentPlan> neighbourhood = new ArrayList<TreatmentPlan>();
		 front.add(initialSolution);
		 notVisited.add(initialSolution);
		 
		 boolean allVisited=false;
		 int iteraciones = 0;
	        //Mientras se puedea expandir la frontera
	        
		while(allVisited==false && iteraciones <=30) {
			System.out.println(iteraciones+" "+front.size());
			if(iteraciones%1==0) {
				System.out.println("Start Iteration: "+iteraciones);
				frontPrint(front);
				System.out.println("End Iteration: "+iteraciones);
            }
			iteraciones++;
            asoPhase(notVisited, front, visited);
            //frontPrint(front);
            
            awoPhase(front);
            //frontPrint(front);
            NDP(front, neighbourhood,visited,notVisited);
            //frontPrint(front);
            if(visited.containsAll(front)) {
            	allVisited=true;
            }           
            
            
           //call solver with same weigths
         
		}
		frontPrint(front);
		//navigationPhase(front);
		return front;
		 
	 }
	 
	 
	public surface ParetoSurfaceNavigation(DDM m, TreatmentPlan initialsolution) {
		// TODO Auto-generated method stub
		ArrayList<TreatmentPlan> front= new ArrayList<TreatmentPlan>();
		front.add(initialsolution);
		preSurfaceGeneration(front);
		
		surface originalSurface=new surface(front.get(0),front.get(1),front.get(2));
		//originalSurface.printInitialList();
		//originalSurface.generatePoint(0, 2);
		//originalSurface.generatePoint(1, 2);
		originalSurface.printAllList();

		return originalSurface;
	}
	
	public surface FMOParetoSurfaceNavigation(DDM m, TreatmentPlan initialsolution1, TreatmentPlan initialsolution2) {
		// TODO Auto-generated method stub
		ArrayList<TreatmentPlan> front= new ArrayList<TreatmentPlan>();
		front.add(initialsolution1);
		front.add(initialsolution2);
		preSurfaceGeneration(front);
		
		surface originalSurface=new surface(front.get(0),front.get(1),front.get(2));
		originalSurface.printInitialList();
		originalSurface.generatePoint(0, 2);
		originalSurface.generatePoint(1, 2);
		originalSurface.printAllList();

		return originalSurface;
	}
	
	public ArrayList<TreatmentPlan> AttachedParetoSurfaceNavigation(DDM m, ArrayList<TreatmentPlan> initialsolutions) {
		// TODO Auto-generated method stub
		ArrayList<surface> initialSurfaces= new ArrayList<surface>();
		for(TreatmentPlan initialTreatment:initialsolutions) {
			ArrayList<TreatmentPlan> front= new ArrayList<TreatmentPlan>();
			front.add(initialTreatment);
			preSurfaceGeneration(front);
			surface originalSurface=new surface(front.get(0),front.get(1),front.get(2));
			
		}
		
		AttachedSurface newAttachedSurface=new AttachedSurface(initialSurfaces);
		
		return newAttachedSurface.treatmentPlanFront;
	}
	public ArrayList<TreatmentPlan> AttachedParetoSurfaceNavigationTest(DDM m) throws ClassNotFoundException, IOException {
		// TODO Auto-generated method stub
		ArrayList<TreatmentPlan> initialsolutions= new ArrayList<TreatmentPlan>();
		initialsolutions.add(serializeDataIn("tpPen-05-bac1.tpg"));
		initialsolutions.add(serializeDataIn("tpPen-05-bac4.tpg"));
		
		ArrayList<surface> initialSurfaces= new ArrayList<surface>();
		for(TreatmentPlan initialTreatment:initialsolutions) {
			ArrayList<TreatmentPlan> front= new ArrayList<TreatmentPlan>();
			front.add(initialTreatment);
			preSurfaceGeneration(front);
			surface originalSurface=new surface(front.get(0),front.get(1),front.get(2));
			initialSurfaces.add(originalSurface);
		}
		
		AttachedSurface newAttachedSurface=new AttachedSurface(initialSurfaces);
		
		return newAttachedSurface.intersectionListAttached;
	}
	public void AttachedParetoSurfaceNavigation2(DDM m, ArrayList<surface> initialsurface) {
		// TODO Auto-generated method stub
		
	}
	 private void preSurfaceGeneration(ArrayList<TreatmentPlan> front) {
		 ArrayList<double[]> weigths=extremeWeightTwoObjective();
		 ArrayList<TreatmentPlan> nav=new ArrayList<TreatmentPlan>();
		 
		 for(TreatmentPlan node:front) {
			 nav=generateNavFront(node,weigths);
			 awoPhase(nav,1);
			 nav.add(DAOCrossingMiddleGeneration(nav.get(0), nav.get(1)));
			 
			 
		 }
		 front.clear();
		 front.addAll(nav);
		 
	 }
	 

	 
	 /**
	  * This function made the first threpoint for each node in front
	  * @param front
	  */
	 private void navigationPhase(ArrayList<TreatmentPlan> front) {
		 ArrayList<double[]> weigths=weightGenerator2o();
		 ArrayList<TreatmentPlan> nav=new ArrayList<TreatmentPlan>();
		 System.out.println("--Navigation--");
		 for(TreatmentPlan node:front) {
			 nav=generateNavFront(node,weigths);
			 awoPhase(nav,1);
			 frontPrint(nav);
			 fmoCrossingGenerator(nav.get(0),nav.get(1));
			 System.out.println("-----");
		 }
		 
		 
	 }
	 private void navigationPhase3d(ArrayList<TreatmentPlan> front) {
		 ArrayList<double[]> weigths=weightGenerator3d();
		 ArrayList<TreatmentPlan> nav=new ArrayList<TreatmentPlan>();
		 for(TreatmentPlan node:front) {
			 nav=generateNavFront(node,weigths);
			 awoPhase(nav,1);
			 frontPrint(nav);
			 fmoCrossingGenerator3d(nav.get(0),nav.get(1),nav.get(2));
			 System.out.println("-----");
		 }
		 
		 
	 }
	 
	 private ArrayList<TreatmentPlan> generateNavFront(TreatmentPlan node,ArrayList<double[]> weigths) {
		 ArrayList<TreatmentPlan> newTreatments=new ArrayList<TreatmentPlan>();
		 dephtCopy(node,newTreatments,weigths);
		 return newTreatments;
	 }
	 

	private void asoPhase(ArrayList<TreatmentPlan> notVisited,ArrayList<TreatmentPlan> front,ArrayList<TreatmentPlan> visited){
		 
		 ArrayList<TreatmentPlan> neighbourhood = new ArrayList<TreatmentPlan>();
		 
		 for(TreatmentPlan nodo : notVisited ) {//Se genera la vecindad de los nodos no visitados en la frontera
             if(!visited.contains(nodo)) {
             	//We generate the neigbourhood from each node not visited
                 neighbourhood(nodo, neighbourhood,2);
                 //add the node in the visited list
                 visited.add(nodo);
                 
                 //NDP(front, neighbourhood,visited);
                 
             }
             //notVisited.clear();

         }
         frontPrint(neighbourhood);
       //then, we update the front this is different from what angel does. is like paquete does
         NDP(front, neighbourhood,visited,notVisited);
         
		 
	 }
	 /**
	  * This function run the solver for each treatment plan in the front using threads
	  * @param front
	  * @param M
	  */
	 private void awoPhase(ArrayList<TreatmentPlan> frontera){
		//Se crea una copia de la frontera, esta es evaluada con el solver y luego se junta con la frontera original 
		ArrayList<TreatmentPlan> front=new ArrayList<TreatmentPlan>();
		depthCopy(frontera,front,weightGenerator());
		int cores = Runtime.getRuntime().availableProcessors();
		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();
		ExecutorService pool2 = Executors.newFixedThreadPool(cores);
		for(int i=0;i<front.size();i++) {
			front.get(i).evaluationFunction =2;
			calls.add(Executors.callable(front.get(i)));
				
		}
		 try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		for(int i=0;i<front.size();i++) {
			front.get(i).evaluationFunction =0;			
		}
		List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
		ExecutorService poolThread = Executors.newFixedThreadPool(front.size());

		for(int i=0;i<front.size();i++) {
			ncalls.add(Executors.callable(front.get(i)));
			
		}
		
		try {
			poolThread.invokeAll(ncalls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		poolThread.shutdown();
		frontera.addAll(front);
		deleteSameValue(frontera);
		
	 }
	 /**
	  * This function run the solver for each treatment plan in the front using threads
	  * @param front
	  * @param M
	  */
	 private void awoPhase(ArrayList<TreatmentPlan> front,int a){
		//Se crea una copia de la frontera, esta es evaluada con el solver y luego se junta con la frontera original 

		int cores = Runtime.getRuntime().availableProcessors();
		List<Callable<Object>> calls = new ArrayList<Callable<Object>>();
		ExecutorService pool2 = Executors.newFixedThreadPool(cores);
		for(int i=0;i<front.size();i++) {
			front.get(i).evaluationMethod =3;
			front.get(i).evaluationFunction =3;
			calls.add(Executors.callable(front.get(i)));
				
		}
		 try {
			 pool2.invokeAll(calls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		pool2.shutdown();
		for(int i=0;i<front.size();i++) {
			front.get(i).evaluationMethod =1;
			front.get(i).evaluationFunction =3;			
		}
		List<Callable<Object>> ncalls = new ArrayList<Callable<Object>>();
		ExecutorService poolThread = Executors.newFixedThreadPool(front.size());

		for(int i=0;i<front.size();i++) {
			ncalls.add(Executors.callable(front.get(i)));
			
		}
		
		try {
			poolThread.invokeAll(ncalls);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		poolThread.shutdown();
		
	 }

	 
	 
	 private void frontPrint(ArrayList<TreatmentPlan> front) {
		// TODO Auto-generated method stub
		 System.out.println("");
		 System.out.println("--------------");
		 for(TreatmentPlan nodo:front) {
			 nodo.scorePrint();
			 System.out.println("");
		 }
		
	}
	 private void frontPrintALL(ArrayList<TreatmentPlan> front) {
			// TODO Auto-generated method stub
			 System.out.println("");
			 System.out.println("--------------");
			 for(TreatmentPlan nodo:front) {
				 nodo.scorePrint();
				 System.out.print(" ");
				 System.out.println(nodo.getLowIntensityperBeam()+" "+nodo.getLowIntensityAngle()+" "+nodo.beamOnTime());
			 }
			
		}
	 private void frontPrintWeigthObjectiveFunction(ArrayList<TreatmentPlan> front,double[] weight) {
		// TODO Auto-generated method stub
		 System.out.println("--------------");
		 for(TreatmentPlan nodo:front) {
			 System.out.println(nodo.evaluationFunction);
		 }
		
	}
	/**
	  * In this version we add the solver solution
	  * @param initialSolution
	  * @return
	  */
	 public ArrayList<TreatmentPlan> ParetoLocalSearch(TreatmentPlan initialSolution) {
		 ArrayList<TreatmentPlan> front= new ArrayList<TreatmentPlan>();
		 ArrayList<TreatmentPlan> notVisited= new ArrayList<TreatmentPlan>();
		 ArrayList<TreatmentPlan> visited= new ArrayList<TreatmentPlan>();
		 ArrayList<TreatmentPlan> neighbourhood = new ArrayList<TreatmentPlan>();
		 front.add(initialSolution);
		 
		 boolean allVisited=false;
		 int iteraciones = 0;
	        //Mientras se puedea expandir la frontera
	        
		while(!allVisited || iteraciones ==30) {
        	
			iteraciones++;
            for(TreatmentPlan nodo : notVisited ) {//Se genera la vecindad de los nodos no visitados en la frontera
                if(!visited.contains(nodo)) {
                	//We generate the neigbourhood from each node not visited
                    neighbourhood(nodo, neighbourhood,1);
                    //add the node in the visited list
                    visited.add(nodo);
                    
                    //NDP(front, neighbourhood,visited);
                    
                }
                

            }
          //then, we update the front this is different from what angel does. is like paquete does
            NDP(front, neighbourhood,visited,notVisited);
            
            
            if(visited.containsAll(front)) {
            	
            }
            
            
            //If all the elements of the front has been visit We stop the algorithm, else We continue with the next iteration.
            if(visited.containsAll(front))
                allVisited=true;
		}
		 
		return front;
		 
	 }
	public ArrayList<double[]> weightGenerator(){
		ArrayList<double[]> weights=new ArrayList<double[]>();
		double a[]={5,1,0};// +PTV
		double b[]={1,1,0};// =
		double c[]={1,5,0};// +Rectum
		weights.add(a);
		weights.add(b);
		weights.add(c);
		
		
		return weights;
	}
	public ArrayList<double[]> extremeWeightTwoObjective(){
		ArrayList<double[]> weights=new ArrayList<double[]>();
		double a[]={0,5,1};// +PTV
		double c[]={0,1,5};// +Rectum
		weights.add(a);
		weights.add(c);
		
		
		return weights;
	}
	public ArrayList<double[]> weightGenerator3d(){
		ArrayList<double[]> weights=new ArrayList<double[]>();
		double a[]={5,1,1};// +PTV
		double b[]={1,5,5};// =
		double c[]={3,2,2};// +Rectum
		weights.add(a);
		weights.add(b);
		weights.add(c);
		
		
		return weights;
	}
	public ArrayList<double[]> weightGenerator2o(){
		ArrayList<double[]> weights=new ArrayList<double[]>();
		double a[]={5,1,1};// +PTV
		double b[]={1,5,5};// +rectum +bladder
		
		weights.add(a);
		weights.add(b);

		
		return weights;
	}
	public ArrayList<double[]> weightGenerator2(){
		ArrayList<double[]> weights=new ArrayList<double[]>();
		double a[]={5,1,1};// +PTV
		double b[]={4,2,2};// =
		double c[]={3,3,3};// +Rectum
		weights.add(a);
		weights.add(b);
		weights.add(c);
		
		
		double d[]={2,4,4};
		weights.add(d);
		double f[]={1,5,5};
		weights.add(f);

		
		
		return weights;
	}
	public ArrayList<double[]> weightGenerator3(){
		ArrayList<double[]> weights=new ArrayList<double[]>();
		double a[]={5,1,1};// +PTV
		weights.add(a);

		return weights;
	}
	public ArrayList<double[]> weightGeneratorTestPTV(){
		ArrayList<double[]> weights=new ArrayList<double[]>();
		double a[]={5,1,0};// +PTV
		weights.add(a);

		return weights;
	}
	public ArrayList<double[]> weightGeneratorTestMO(){
		ArrayList<double[]> weights=new ArrayList<double[]>();
		double a[]={5,1,0};// +PTV
		double b[]={4,1,0};//
		double c[]={3,1,0};//
		double d[]={2,1,0};// =
		double e[]={1,1,0};//
		double f[]={1,2,0};//
		double g[]={1,3,0};//
		double h[]={1,4,0};//
		double i[]={1,5,0};// +Rectum
		weights.add(a);
		weights.add(b);
		weights.add(c);
		weights.add(d);
		weights.add(e);
		weights.add(f);
		weights.add(g);
		weights.add(h);
		weights.add(i);

		return weights;
	}
	public ArrayList<TreatmentPlan> adaptableParetoLocalSearch(TreatmentPlan initialSolution) {
			return null;
			 
	 }
	 public ArrayList<TreatmentPlan> kmeansParetoLocalSearch(TreatmentPlan initialSolution) {
			return null;
			 
	 }
	 
	 public ArrayList<TreatmentPlan> convexHull(ArrayList<TreatmentPlan> neigborhood){
		 ArrayList<TreatmentPlan> paretoFront= new ArrayList<TreatmentPlan>();
		 return  paretoFront;
	 }
	 /**
	  * This method fresh the front list, 
	  * @param front
	  * @param neighbourhood
	  * @param visited
	  */
	 private void NDP(ArrayList<TreatmentPlan> front, ArrayList<TreatmentPlan> neighbourhood,ArrayList<TreatmentPlan> visited,ArrayList<TreatmentPlan> notVisited) {
		
		 front.addAll(neighbourhood);
		 
		 ArrayList<TreatmentPlan> dominatedList= new ArrayList<TreatmentPlan>();
		 

        for(TreatmentPlan  p: front ) {

            for(TreatmentPlan q : front ) {
            	
                if(dominated(p,q)&& p!=q) {
                	if(dominatedList.contains(p)) {
                		break;
                	}
                	else {
                		dominatedList.add(p);
                		break;
                	}
                	
                	
                }

            }

        }
        //refresh the not visited list
        notVisited.clear();
        notVisited.addAll(front);
        notVisited.removeAll(visited);
        //then we eliminated the dominated points of the front
		front.removeAll(dominatedList);
		 
		
	}
	 
	 /**
	  * TODO: add the other leaf random implementen in the class algorithm
	  * @param node
	  * @param neighbourhood
	  * @param mov
	  */
	 public void neighbourhood(TreatmentPlan node, ArrayList<TreatmentPlan> neighbourhood , int mov) {
		 neighbourhood.clear();
		 switch(mov) {
			
			case 1:
				movAperture2(neighbourhood,node,M);
			
				break;
				
			case 2:
				movLeafsRandomRevision(neighbourhood,node,M,24);
				break;
				
			case 3:
				movLeafsAperture(neighbourhood,node,M);
				movLeafsRandom(neighbourhood,node,M,24);
				break;
			
			default:
				movLeafsRandom(neighbourhood,node,M,24);
				break;
			
			}
	 }
	 
	 /**
	  * This function return true if q solution is dominated by p
	  * @param p
	  * @param q
	  * @return
	  */
	 public boolean dominated(TreatmentPlan p,TreatmentPlan q) {

		
		//return p.isDominatedPtvRectum(q);
		return p.isDominatedTwoObjective(q);

	 }
	 

	 public void orderFront(ArrayList<TreatmentPlan> neigborhood){
		 
	}
	 public void depthCopy(ArrayList<TreatmentPlan> original,ArrayList<TreatmentPlan> copy) {
		 for(TreatmentPlan node:original) {
			 TreatmentPlan newTreatment= new TreatmentPlan(node);
			 copy.add(newTreatment);
		 }
		 
	 }
	 /**
	  * From each node in the ArrayList Original a set of new treatmentplas are made it. 
	  * The number of treatmentPlan generated depend of the size from setweight
	  * @param original
	  * @param copy
	  * @param setWeight
	  */
	 public void depthCopy(ArrayList<TreatmentPlan> original,ArrayList<TreatmentPlan> copy,ArrayList<double[]> setWeight) {
		 for(TreatmentPlan node:original) {
			 for(double[] weight:setWeight) {
				 TreatmentPlan newTreatment= new TreatmentPlan(node);
				 newTreatment.weights=weight;
				 copy.add(newTreatment);
				 
			 }
			 
		 }
		 
	 }
	 private void dephtCopy(TreatmentPlan node, ArrayList<TreatmentPlan> newTreatments, ArrayList<double[]> weigths) {
		// TODO Auto-generated method stub
		 for(double[] weight:weigths) {
			 TreatmentPlan newTreatment= new TreatmentPlan(node);
			 newTreatment.weights=weight;
			 newTreatments.add(newTreatment);
			 
		 }
		
	}
	private ArrayList<TreatmentPlan> weightCopy(TreatmentPlan node, ArrayList<double[]> weights) {
			// TODO Auto-generated method stub
		ArrayList<TreatmentPlan> newTreatments=new ArrayList<TreatmentPlan>();
		for(double[] weight:weights) {
			 TreatmentPlan newTreatment= new TreatmentPlan(node);
			 newTreatment.weights=weight;
			 newTreatments.add(newTreatment);
				 
		 }
		return newTreatments;
			
	}
	 
	private ArrayList<double[]> fmoCrossingGenerator(TreatmentPlan nodeOne,TreatmentPlan nodeTwo){
		double[] fmo1,fmo2,fmo3,fmo4,fmo5;
		//ArrayList<double[]> scores
		fmo1=nodeOne.intensity;
		//nodeOne.printIntensity();
		fmo2=nodeTwo.intensity;
		//nodeTwo.printIntensity();
		fmo3=fmoCross(fmo1,fmo2);
		evaluateSolution(fmo3, M, o, nodeOne.dd);
		fmo4=fmoCross(fmo1,fmo3);
		evaluateSolution(fmo4, M, o, nodeOne.dd);
		fmo5=fmoCross(fmo3,fmo2);
		evaluateSolution(fmo5, M, o, nodeOne.dd);
		
		return null;
		
		
	}
	private ArrayList<double[]> fmoLinearCombination(TreatmentPlan nodeOne,TreatmentPlan nodeTwo, double percent){
		double[] fmo1,fmo2,fmo3;
		//ArrayList<double[]> scores
		fmo1=nodeOne.intensity;
		//nodeOne.printIntensity();
		fmo2=nodeTwo.intensity;
		//nodeTwo.printIntensity();
		fmo3=fmoCross(fmo1,fmo2);
		evaluateSolution(fmo3, M, o, nodeOne.dd);

		
		return null;
		
		
	}
	
	private TreatmentPlan DAOCrossingMiddleGeneration(TreatmentPlan one,TreatmentPlan two) {;
		return DAOCrossingGeneration(one,two, (double) 0.5);
		
	}
	private TreatmentPlan DAOCrossingGeneration(TreatmentPlan one,TreatmentPlan two, double percent) {
		TreatmentPlan crossDAO=new TreatmentPlan(one);
		crossDAO.setIntensity(intensitiesLinearCombination(one.getIntensityApertures(), two.getIntensityApertures(), percent));
		crossDAO.evaluateSolution();
		return crossDAO;
		
	}
	
	
	/**
	 * This function made a linear combination of the two intensities
	 * Consider that the function only work with itnensities that have the same number of bac and same number  of apertures for all bacs
	 * @param a
	 * @param b
	 * @param alfa
	 * @return
	 */
	private double [][] intensitiesLinearCombination(double [][] a,double [][] b,double alfa){
		double [][]linearCombination=new double [a.length][a[0].length];
		for(int i=0;i<a.length;i++) {
			for(int j=0;j<a[i].length;j++) {
				linearCombination[i][j]=(alfa*a[i][j])+((1-alfa)*b[i][j]);
			}
		}
		return linearCombination;
	}
	
	
	private ArrayList<double[]> FMOCrossingGenerator(TreatmentPlan nodeOne,TreatmentPlan nodeTwo, int numberOfCrossing){
		double[] fmo1,fmo2,fmo3,fmo4,fmo5;
		//ArrayList<double[]> scores
		fmo1=nodeOne.intensity;
		//nodeOne.printIntensity();
		fmo2=nodeTwo.intensity;
		//nodeTwo.printIntensity();
		for(int i=0;i<numberOfCrossing;i++) {
			
			
			
		}
		
		fmo3=fmoCross(fmo1,fmo2);
		evaluateSolution(fmo3, M, o, nodeOne.dd);
		fmo4=fmoCross(fmo1,fmo3);
		evaluateSolution(fmo4, M, o, nodeOne.dd);
		fmo5=fmoCross(fmo3,fmo2);
		evaluateSolution(fmo5, M, o, nodeOne.dd);
		
		return null;
		
		
	}
	
	
	private ArrayList<double[]> fmoCrossingGenerator3d(TreatmentPlan nodeOne,TreatmentPlan nodeTwo, TreatmentPlan nodeThree){
		double[] fmo1,fmo2,fmo3,fmo4,fmo5,fmo6,fmo7,fmo8,fmo9;
		//ArrayList<double[]> scores
		fmo1=nodeOne.intensity;
		fmo2=nodeTwo.intensity;
		fmo3=nodeThree.intensity;
		
		
		fmo4=fmoCross(fmo1,fmo2);
		evaluateSolution(fmo4, M, o, nodeOne.dd);
		fmo5=fmoCross(fmo1,fmo3);
		evaluateSolution(fmo5, M, o, nodeOne.dd);
		fmo6=fmoCross(fmo3,fmo2);
		evaluateSolution(fmo6, M, o, nodeOne.dd);
		
		fmo7=fmoCross(fmo4,fmo5);
		evaluateSolution(fmo7, M, o, nodeOne.dd);
		fmo8=fmoCross(fmo4,fmo6);
		evaluateSolution(fmo8, M, o, nodeOne.dd);
		fmo9=fmoCross(fmo6,fmo5);
		evaluateSolution(fmo9, M, o, nodeOne.dd);
		
		
		return null;
		
		
	}
	
	private double[] fmoCross(double[] fmoOne,double[] fmoTwo){
		double[] newFMO = new double[fmoOne.length];
		for (int i = 0; i < newFMO.length; i++) {
			newFMO[i] = (fmoOne[i] + fmoTwo[i])/2;
		}
		return newFMO;
		
	}
	/**
	 * Cross with the percent of of the fmoOne
	 * @param fmoOne
	 * @param fmoTwo
	 * @param percent
	 * @return
	 */
	private double[] fmoCross(double[] fmoOne,double[] fmoTwo, double percent){
		double scalar1,scalar2;
		scalar1=percent;
		scalar2=1-scalar1;
		double[] newFMO = new double[fmoOne.length];
		
		for (int i = 0; i < newFMO.length; i++) {
			newFMO[i] = (scalar1*fmoOne[i] + scalar2*fmoTwo[i]);
		}
		return newFMO;
		
	}

	
	public double[]  evaluateSolution(double[] solutionVector, DDM M, Organs[] o, double[] dd){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		
		double F = 0.0, pen;
		double[] score= {0,0,0};
		
		

		
		//System.out.println("Largo solución :"+solutionVector.length);
		//Recorremos los organos
		for(int i=0; i<3; i++){
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
            		pen +=  Math.pow((dd[i] - intensityVoxel2Beams),2);
            	}else{
   					pen +=  Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				}
            	
	        }
			//System.out.println(i+" "+Zmax[i]+" ");
			score[i]=(pen/aux_index.size());
			//System.out.print(score[i]+" ");
			
		}
		System.out.println("");
		//System.out.println("Solucion: "+F);
		return score;
	}
	void deleteSameValue(ArrayList<TreatmentPlan> Front) {
        // Create a new ArrayList
        ArrayList<TreatmentPlan> newList = new ArrayList<TreatmentPlan>();
        newList.add(Front.get(0));
        // Traverse through the first list
        for (int j=1;j<Front.size();j++) {
  
            // If this element is not present in newList
            // then add it
        	boolean isInTheList=false;
        	for(int i=0;i<newList.size();i++) {
	            if (newList.get(i).isSameValue(Front.get(j))) {
	            	isInTheList=true;
	            	break;
	            }
        	}
        	if(isInTheList==false)newList.add(Front.get(j));
        }
  
        // return the new list
        Front.clear();
        Front.addAll(newList);
	}


	
	/** Funciones que solo estan por pruebas, en versiones posteriores se deben eliniar **/
	public static TreatmentPlan serializeDataIn(String fileName) throws IOException, ClassNotFoundException{
		   
		   FileInputStream fin = new FileInputStream(fileName);
		   ObjectInputStream ois = new ObjectInputStream(fin);
		   TreatmentPlan saveTreatmentPlan= (TreatmentPlan) ois.readObject();
		   ois.close();
		   return saveTreatmentPlan;
		}
	
	
	
}
