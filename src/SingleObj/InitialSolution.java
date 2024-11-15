package SingleObj;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ApertureBase.sequentialGenerator;
import IMRT_DAO.*;
import gurobi.GRBException;
public class InitialSolution {
	DDM M;
	Organs[] o;
	double[] dd;
	int[] beamletAngles;
	int init_intensity;
    int max_intensity;
    
    
    
    
    
	int max_delta;
	int max_iter;int max_time; 
	int seed; 
	int step_intensity; 
	int numAngles; 
	int numOrgans;
	int numAperture, rounded;
	public int totalbmlt;
	public int[] selAngles;
	public double[] w;
    public int[][] beamletSet;
    public int[] initialBACs;
    String pathFile;
    TreatmentPlan solution;
    public Vector<int[][]> initial_aperture_shape;  // # DAO : Conjunto de formas de las aperturas según restricciones (rango sup y rango inf)
    public Vector<int[][]> apertureShapes; 
    public Gurobi_Solver newModel;
    public sequentialGenerator leafSequencing;
	public InitialSolution(){
		
	}
	
	public InitialSolution(double[] w,int[][] beamletSet,int[] initialBACs,String pathFile,DDM M, Organs[] o, double[] dd, int[] beamletAngles, int init_intensity, int max_intensity,
			int max_delta, int max_iter, int max_time, int seed, int step_intensity, int numAngles, int numOrgans) {
		// TODO Auto-generated constructor stub
		this.w=w;
		this.pathFile=pathFile;
		this.M=M;
		this.o=o;
		this.dd=dd;
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
		solution = new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans,o,dd,M,beamletAngles);

	    solution.setWeights(w);
	    try {
			setBeams(solution);
			setInitialApertureShape();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
	}

	public TreatmentPlan GetInitialSolution() {
		heuristicAperture();
		solverIntensity();
		
		
		return solution;
		
	}
	
	public TreatmentPlan GetInitialSolution2() {
		openAperture();
		//solution.evaluateSolution();
		solverIntensity();
		
		return solution;
		
	}
	public TreatmentPlan GetInitialSolution3() {
		closeAperture();
		solverIntensity();
		
		return solution;
		
	}
	public TreatmentPlan GetInitialSolution4() {
		randomAperture();
		solverIntensity();
		
		return solution;
		
	}
	public TreatmentPlan GetInitialSolution5() {
		centerAperture();
		solverIntensity();
		
		return solution;
		
	}
	public TreatmentPlan GetInitialSolution6() {
		heuristicAperture2();
		solverIntensity();
		
		return solution;
		
	}
	
	public TreatmentPlan GetInitialSolution7() {
		sequientialAperture();
		
		solution.evaluateSolution();
		//System.out.print(solution.singleObjectiveValue+" ");
		solverIntensity2();
		//System.out.print(solution.singleObjectiveValue+" ");
		//System.out.print(solution.getLowIntensityperBeam()+" ");
		
		//bmlt();
		return solution;
		
	}
	
	//Test 
	
	public TreatmentPlan test1() {
		openAperture();
		solverIntensity(0);
		return solution;
		
	}
	public TreatmentPlan test2() {
		openAperture();
		solverIntensity(1);
		return solution;
		
	}
	public TreatmentPlan test3() {
		heuristicAperture();
		solverIntensity(0);
		return solution;
		
	}
	public TreatmentPlan test4() {
		heuristicAperture();
		solverIntensity(1);
		return solution;
		
	}
	public TreatmentPlan test4gEUD() {
		solution.evaluationFunction=3;
		heuristicAperture();
		solverIntensitygEUD();
		return solution;
		
	}
	public TreatmentPlan test5() {
		heuristicAperture2();
		solverIntensity(0);
		return solution;
		
	}
	public TreatmentPlan test6() {
		heuristicAperture2();
		solverIntensity(1);
		return solution;
		
	}
	public TreatmentPlan generateSolution() {
		try {
			readSolution(solution);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return solution;
		
	}


	
	
	
	public void bmlt() {
		for(int i=0;i<solution.aperturesBmlts.size();i++) {
			for(int j=0;j<solution.aperturesBmlts.get(i).length;j++) {
				System.out.println("param x"+(i+1)+(j+1)+" :=");
				for(int k=0;k<solution.aperturesBmlts.get(i)[j].length;k++) {
					System.out.println((k+1)+"	"+solution.aperturesBmlts.get(i)[j][k] );
				}
				System.out.println(";");
				System.out.println("");
			}
		}
	}
	
	/**
	 * Divide la intensidad en todas a la aperturas de forma proporcional por cada anglulo
	 * TODO: modificar para que dependa de la cantidad de aperturass del treatment plan
	 * @param min
	 * @return
	 */
	public void openFixSolution(int min) {
		for(int i=0; i<solution.apertures.size(); i++){
			double aux=solution.getTotalIntensityperAngle(i);
			aux=aux/5;
			if(aux<min) {
				aux=min;
			}
			for(int j=0;j<solution.apertures.get(i).size();j++) {
				solution.setIntensity(i, j, aux);
			}
		
		}
		solution.evaluateSolution();
	
	}
	public TreatmentPlan randomSolution() {
		
		ArrayList<TreatmentPlan> population=generatePopulation(10);
		
		
		return solution;
	}
	
	private ArrayList<TreatmentPlan> generatePopulation(int size) {
		ArrayList<TreatmentPlan> newPopulation=new ArrayList<TreatmentPlan>();
		for(int i=0;i<size;i++) {
			TreatmentPlan newsolution = new TreatmentPlan(init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans,o,dd,M,beamletAngles);
			newsolution.setWeights(w);
			generateSolutionRandomUniform(newsolution);
			randomIntensity(newsolution);
			newPopulation.add(newsolution);
		}
		

		return newPopulation;
	}

	public TreatmentPlan randomFixSolution() {
		randomAperture();
		solverIntensity();
		int []numberOfApertureRemove=new int[numAngles];
		int actual=20;
		int better=solution.getLowIntensityperBeam();
		solution.evaluateSolution();
		while(better>5) {
			//System.out.println(better);
			List<Callable<Object>> calls = new ArrayList<Callable<Object>>();
			ExecutorService pool2 = Executors.newFixedThreadPool(16);
			for(int i=0;i<numAngles;i++){
				
				//Eliminar aperturas con malas intensidades
				numberOfApertureRemove[i]=solution.removeBadAperturePerBeam(i, 1);
				//agregar nuebas aperturas//
				
				
			}
			ArrayList<TreatmentPlan> pool =generateListRandomApertures(solution,numberOfApertureRemove);
			
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
				actual=pool.get(i).getLowIntensityperBeam();
				if(actual<better) {
					solution=pool.get(i);
					better=actual;
				}
			} 
			
		}
		

		
		return solution;
		
	}
	
	

	/**
	 * genera una aperutra random en todas las aperturas que tengan una intensidad 0
	 */
	public void fixSolutionRandom() {
		int numberOfApertureRemove=0;
		
		while(solution.getLowIntensityperBeam()>0) {
		
			for(int i=0;i<numAngles;i++){
				
				//Eliminar aperturas con malas intensidades
				numberOfApertureRemove=solution.removeBadAperturePerBeam(i, 1);
				//agregar nuebas aperturas//
				generateRandomAperture(solution,i,numberOfApertureRemove);
				
			}
			
			
			
			
		}
		solverIntensity();
		
		//System.out.print("");

		
	}
	public void fixSolutionSQR() {
		int numberOfApertureRemove=0;
		double[] fatherOne;
		double[] fatherTwo;
		for(int i=0;i<numAngles;i++){
			//obtenemos los sqf de cada apertura
			 Vector<double[]> v =solution.indexSQFperBeam(1, 1, i, beamletAngles);
			orderVector(v);
			//generamos las aperturas hijos
			
			
			fatherOne=v.lastElement();
			fatherTwo=v.get(v.size()-2);
			ArrayList<int[][]> childsAperture=generateGeneticAperture(i,(int)fatherOne[0],(int)fatherTwo[0]);
			//Eliminar 2 aperturas con las intensidades mas bajas / despues probar con los sqf mas bajos
			solution.removeAperture(i,v);
			//agregar las nuevas aperturas
			
			Aperture aper1 = new Aperture(childsAperture.get(0), 1);
			solution.apertures.get(i).add(aper1);
			Aperture aper2 = new Aperture(childsAperture.get(1), 1);
			solution.apertures.get(i).add(aper2);

			
		}
		solverIntensity();
		


		
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
	
	
	public void randomIntensity(TreatmentPlan sol) {
		Random rd = new Random(); 
		double[][] newIntensities=new double[numAngles][];
		for(int i=0;i<numAngles;i++) {
			newIntensities[i]=rd.doubles(numAperture,0,max_intensity).toArray();
		}
		sol.setIntensity(newIntensities);
	}
	public void solverIntensity() {
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
	public void solverIntensitygEUD() {
	 	solution.aperturesToIntensityMatrixPerStationwithoutIntensty();
	    solution.ApertureMatrixStationToIntensityMatrixPlan();
	

	    //time=System.currentTimeMillis();
	    
	    Gurobi_Solver newModel;
		try {
			newModel = new Gurobi_Solver(solution,M,o,beamletAngles,dd,solution.weights,0,1);
			solution.setIntensity(newModel.newIntensity);
		    solution.setOF(newModel.objVal);
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
}
	public void solverIntensity2() {
	 	solution.aperturesToIntensityMatrixPerStationwithoutIntensty();
	    solution.ApertureMatrixStationToIntensityMatrixPlan();
	

	    //time=System.currentTimeMillis();
	    
	    Gurobi_Solver newModel;
		try {
			newModel = new Gurobi_Solver(solution,M,o,beamletAngles,dd,solution.weights,0);
			solution.setIntensity(newModel.newIntensity);
		    solution.setOF(newModel.objVal);
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
}
	
	public void solverIntensity(int min) {
	 	solution.aperturesToIntensityMatrixPerStationwithoutIntensty();
	    solution.ApertureMatrixStationToIntensityMatrixPlan();
	

	    //time=System.currentTimeMillis();
	    
	    Gurobi_Solver newModel;
		try {
			newModel = new Gurobi_Solver(solution,M,o,beamletAngles,dd,solution.weights,min);
			solution.setIntensity(newModel.newIntensity);
		    solution.setOF(newModel.objVal);
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
}
	public void solverIntensityPen(int min) {
	 	solution.aperturesToIntensityMatrixPerStationwithoutIntensty();
	    solution.ApertureMatrixStationToIntensityMatrixPlan();
	

	    //time=System.currentTimeMillis();
	    
	    Gurobi_Solver newModel;
		try {
			newModel = new Gurobi_Solver(solution,M,o,beamletAngles,dd,solution.weights,1,min);
			solution.setIntensity(newModel.newIntensity);
		    solution.setOF(newModel.objVal);
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
}
	
	public void heuristicAperture() {
		 generateMixFirstSolution(solution);
	}
	public void heuristicAperture2() {
		 generateMixFirstSolution2(solution);
	}
	public void sequientialAperture() {
		leafSequencing = new sequentialGenerator(totalbmlt, M, o, dd, beamletAngles, numAperture, rounded, numAngles, solution);
		solution.aperturesLimit=apertureShapes;
		

	}
	public void openAperture() {
		 generateSolutionOpen(solution);
	}
	public void closeAperture() {
		 generateSolutionClose(solution);
	}
	public void randomAperture() {
		 generateSolutionRandomUniform(solution);
	}
	public void centerAperture() {
		 generateSolutioncenter(solution);
	}
	public Gurobi_Solver returnModel() {
		return newModel;
	}
	
	//--------------Funciones para laobtencion de datos iniciales para el treatment plan
	public int totalBeamlet() {
		int sum=0;
		for(int i=0 ;i<numAngles;i++) {
			sum+=beamletAngles[i];
		}
		return sum;
	}
	public void setBeamlet() {
		for(int i=0 ;i<numAngles;i++) {
			beamletAngles[i]=beamletSet[initialBACs[i]][1];
		}
	}
	
	
	public  void setBeams(TreatmentPlan solution) throws IOException{
		totalbmlt=0;
		for(int i=0; i<numAngles; i++){
			solution.selAngles[i] = new Beam(beamletSet[initialBACs[i]][1],initialBACs[i], pathFile);
			totalbmlt=totalbmlt+solution.selAngles[i].beamlets;
		}
		solution.beamlets=totalbmlt;
	}
	public void setInitialApertureShape(){
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
	
	//------------------ Funciones de Generacion de aperturas
	public void readSolution(TreatmentPlan sol) throws IOException{
		int [][] aper_shape, new_aper;
		int half=0;
		Random r = new Random();
		int random_intensity;
		List<Aperture> list_aperture;
		list_aperture = new ArrayList<Aperture>();
		aper_shape = initial_aperture_shape.get(0);
		new_aper = new int[aper_shape.length][2];
		String sp=" ";
        String path = "./apertura.txt";
        File f = new File(path);
        System.out.println(f.getAbsolutePath());
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        String line = "";
        line=fileIn.readLine();
        int angulo,apertura=0,a,b;
        int flag=0,row=0;
        float intensidad=0;
        while(line != null){
            if (line.contains("Angulo: ")){
                String[] auxReader = line.split(sp);
                angulo = Integer.parseInt(auxReader[1]);
                aper_shape = initial_aperture_shape.get(angulo);
                if(!list_aperture.isEmpty()) {
                	Aperture aper = new Aperture(new_aper, intensidad);
                	list_aperture.add(aper);
                	sol.apertures.add(list_aperture);
                	list_aperture = new ArrayList<Aperture>();
                	flag=0;
        			row=0;
                }
                line=fileIn.readLine();
                
            }
            if (line.contains("Apertura N")){
                String[] auxReader = line.split(sp);
                apertura = Integer.parseInt(auxReader[2]);
                if(flag==1) {
                	Aperture aper = new Aperture(new_aper, intensidad);
        			list_aperture.add(aper);
        			flag=0;
        			row=0;
                }
                
                new_aper = new int[aper_shape.length][2];
                
                line=fileIn.readLine();
               
            }
            if (line.contains("Intensidad: ")){
                String[] auxReader = line.split(sp);
                intensidad = Float.parseFloat(auxReader[1]);
                line=fileIn.readLine();
                
            }
            
            if (line.contains("\t")){
            	
                break;
            }
                
            if (line.contains(",")){
                String[] auxReader = line.split(", ");
                
                a = Integer.parseInt(auxReader[0]);
                b = Integer.parseInt(auxReader[1]);
                new_aper[row][0] = a;
				new_aper[row][1] = b;
				row++;
                flag=1;
                //linne=fileIn.readLine();
                
            }
            if (line.contains("+")){
            	Aperture aper = new Aperture(new_aper, intensidad);
    			list_aperture.add(aper);
                break;
            }
            line=fileIn.readLine();
        }
        Aperture aper = new Aperture(new_aper, intensidad);
		list_aperture.add(aper);
        sol.apertures.add(list_aperture);
        sol.aperturesLimit=apertureShapes;
        sol.evaluateSolution();
		
		
		
		
		
	}
	
	public void generateMixFirstSolution(TreatmentPlan sol){

		int [][] aper_shape, new_aper;
		int half=0;
		Random r = new Random();
		int random_intensity;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		for(int i=0;i<numAngles;i++){
			int sizeY = sol.selAngles[i].maxAbsY;
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			

			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			half=aper_shape.length/2;
			
			//Apertura abierta por arriba
			for(int j=0; j<aper_shape.length; j++){
				new_aper[j][0] = aper_shape[j][0];
				if(j<half) {
					new_aper[j][1] = aper_shape[j][1];
				}else {
					new_aper[j][0] = -1;
					new_aper[j][1] = -1;
				}
			}
			
			Aperture aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);
			
			
			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			half=aper_shape.length/2;
			
			//Apertura abierta por abajo
			for(int j=0; j<aper_shape.length; j++){
				new_aper[j][0] = aper_shape[j][0];
				if(j<half) {
					new_aper[j][0] = -1;
					new_aper[j][1] = -1;
				}else {
					new_aper[j][1] = aper_shape[j][1];
				}
			}
			
			aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);
			
			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			half=sizeY/2;
			
			//centro vertical
			for(int j=0; j<aper_shape.length; j++){
				new_aper[j][0] = aper_shape[j][0];
				if(aper_shape[j][0]>half) {
					new_aper[j][0] = -1;
					new_aper[j][1] = -1;
				}else {
					new_aper[j][1] = half;
				}
			}
			
			 aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);
			
			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			half=sizeY/2;
			
			//centro horizontal
			for(int j=0; j<aper_shape.length; j++){
				new_aper[j][0] = aper_shape[j][0];
				if(aper_shape[j][1]<half) {
					new_aper[j][1] = -1;
					new_aper[j][0] = -1;
				}else {
					if(aper_shape[j][0]<half)new_aper[j][0]=half;
					new_aper[j][1] = aper_shape[j][1];
				}
			}
			
			aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);

			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			
			//centro completo
			for(int j=0; j<aper_shape.length; j++){
				
				new_aper[j][0] = aper_shape[j][0];
	
				new_aper[j][1] = aper_shape[j][1];
			}
			
			aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);
			
		
		sol.apertures.add(list_aperture);
		sol.aperturesLimit=apertureShapes;
		}
	}

	/**
	 * En esta version se cambian las aperturas abierto por el lado y abierto completamente por, centro por vertical y horizontal.
	 * Y una por centro completamente
	 * @param sol
	 */
	public void generateMixFirstSolution2(TreatmentPlan sol){

		int [][] aper_shape, new_aper;
		int half=0;
		Random r = new Random();
		int random_intensity;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		for(int i=0;i<numAngles;i++){
			
			int sizeY = sol.selAngles[i].maxAbsY;
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			int auxRow = ((sol.selAngles[i].maxAbsY)/2)/2;//distancia que debe separar desde el inicio o el final
			int auxColumn=(aper_shape.length/2)/2;//distancia que debe separar desde el inicio o el final
			int up=auxColumn,down=aper_shape.length-auxColumn,right=auxRow,left=(sol.selAngles[i].maxAbsY)-auxRow;
			list_aperture = new ArrayList<Aperture>();
			

			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			half=aper_shape.length/2;
			
			//Apertura abierta por arriba
			for(int j=0; j<aper_shape.length; j++){
				new_aper[j][0] = aper_shape[j][0];
				if(j<half) {
					new_aper[j][1] = aper_shape[j][1];
				}else {
					new_aper[j][0] = -1;
					new_aper[j][1] = -1;
				}
			}
			
			Aperture aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);
			
			
			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			half=aper_shape.length/2;
			
			//Apertura abierta por abajo
			for(int j=0; j<aper_shape.length; j++){
				new_aper[j][0] = aper_shape[j][0];
				if(j<half) {
					new_aper[j][0] = -1;
					new_aper[j][1] = -1;
				}else {
					new_aper[j][1] = aper_shape[j][1];
				}
			}
			
			aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);
			
			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			half=sizeY/2;
			
			//centro vertical
			for(int j=up; j<down; j++){
				new_aper[j][0] = aper_shape[j][0];
				new_aper[j][1] = aper_shape[j][1];
			}
			
			 aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);
			
			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			
			
			//centro horizontal
			for(int j=0; j<aper_shape.length; j++){
				if(aper_shape[j][0]>right) {
					new_aper[j][0] = aper_shape[j][0];
				
				}else {
					new_aper[j][0] = right;
				}
				//centra por la derecha
				if(aper_shape[j][1]<left) {
					new_aper[j][1] = aper_shape[j][1];

				}else {
					new_aper[j][1] = left;
				}
			}
			
			aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);

			new_aper = new int[aper_shape.length][2];
			random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
			
			//centro completo
			for(int j=up; j<down; j++){

				
				//centra por la izquierda
				if(aper_shape[j][0]>right) {
					new_aper[j][0] = aper_shape[j][0];
				
				}else {
					new_aper[j][0] = right;
				}
				//centra por la derecha
				if(aper_shape[j][1]<left) {
					new_aper[j][1] = aper_shape[j][1];

				}else {
					new_aper[j][1] = left;
				}
			}
			
			aper = new Aperture(new_aper, random_intensity);
			list_aperture.add(aper);
			
		
		sol.apertures.add(list_aperture);
		sol.aperturesLimit=apertureShapes;
		}
	}

	
	public  void generateSolutionOpen(TreatmentPlan sol){
		int [][] aper_shape, new_aper;

		Random r = new Random();
		int random_intensity;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		for(int i=0;i<numAngles;i++){
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			for(int z=0;z<5;z++){
				new_aper = new int[aper_shape.length][2];
				random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
				
				
				//Solución inicial completamente abierta
				for(int j=0; j<aper_shape.length; j++){
					new_aper[j][0] = aper_shape[j][0];
					
					new_aper[j][1] = aper_shape[j][1];
				}
				
				Aperture aper = new Aperture(new_aper, 1);
				list_aperture.add(aper);
			}
			sol.apertures.add(list_aperture);
			sol.aperturesLimit=apertureShapes;
		}
		
	}
	

	
	
	
	public void generateSolutionClose(TreatmentPlan sol){
		int [][] aper_shape, new_aper;
		Random r = new Random();
		int random_intensity;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		for(int i=0;i<numAngles;i++){
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			for(int z=0;z<5;z++){
				new_aper = new int[aper_shape.length][2];
				random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 

				//Solución inicial completamente cerrada
				for(int j=0; j<aper_shape.length; j++){
					new_aper[j][0] = -1;
					new_aper[j][1] = -1;
					
				}
				
				Aperture aper = new Aperture(new_aper, random_intensity);
				list_aperture.add(aper);
			}
			sol.apertures.add(list_aperture);
			sol.aperturesLimit=apertureShapes;
		}
		
	}
	
	public  void generateSolutionRandomUniform(TreatmentPlan sol){
		int [][] aper_shape, new_aper;

		Random r = new Random();
		int random_intensity,range;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		for(int i=0;i<numAngles;i++){
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			for(int z=0;z<5;z++){
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
				list_aperture.add(aper);
			}
			sol.apertures.add(list_aperture);
			sol.aperturesLimit=apertureShapes;
		}
		
	}
	
	
	
	private ArrayList<TreatmentPlan> generateListRandomApertures(TreatmentPlan solution, int[] numberOfApertureRemove) {
		
		ArrayList<TreatmentPlan> pool =new ArrayList<TreatmentPlan>(); 
		for(int i=0;i<16;i++) {
			TreatmentPlan newNeigbor=new TreatmentPlan(solution);
			for(int j=0;j<numberOfApertureRemove.length;j++) {
				generateRandomAperture(newNeigbor,j,numberOfApertureRemove[j]);
			}
			newNeigbor.evaluationFunction=2;
			pool.add(newNeigbor);
		}
		return pool;
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
	
	
	public  void generateSolutioncenter(TreatmentPlan sol){
		int [][] aper_shape, new_aper;

		Random r = new Random();
		int random_intensity;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		for(int i=0;i<numAngles;i++){
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			int auxRow = ((sol.selAngles[i].maxAbsY)/2)/2;//distancia que debe separar desde el inicio o el final
			int auxColumn=(aper_shape.length/2)/2;//distancia que debe separar desde el inicio o el final
			int up=auxColumn,down=aper_shape.length-auxColumn,right=auxRow,left=(sol.selAngles[i].maxAbsY)-auxRow;
			for(int z=0;z<5;z++){
				new_aper = new int[aper_shape.length][2];
				random_intensity = r.nextInt((max_intensity / step_intensity) +1) * step_intensity; 
				
				
				
				//Solución inicial completamente centrada
				for(int j=up; j<down; j++){

					
					//centra por la izquierda
					if(aper_shape[j][0]>right) {
						new_aper[j][0] = aper_shape[j][0];
					
					}else {
						new_aper[j][0] = right;
					}
					//centra por la derecha
					if(aper_shape[j][1]<left) {
						new_aper[j][1] = aper_shape[j][1];

					}else {
						new_aper[j][1] = left;
					}
				}
				
				Aperture aper = new Aperture(new_aper, random_intensity);
				list_aperture.add(aper);
			}
			sol.apertures.add(list_aperture);
			sol.aperturesLimit=apertureShapes;
		}
		
	}
	
	/**
	 * This fuction generate new aperture swapping two father apertures.
	 * @param BestSolution
	 * @param angleChoose
	 * @param i
	 * @param j
	 * @return
	 */
	public  ArrayList<int[][]> generateGeneticAperture(int angleChoose, int i, int j){
		//return solution.singlePointCrossover(angleChoose, i, j);
		return solution.twoPointCrossover(angleChoose, i, j);
	}
	
	
	

	
	
	/**
	 * No esta implementada , la idea es que la matriz tenga una mayor apertura en el centro que en los extremos
	 * @param sol
	 */
	public  void generateSolutionRandomNormal(TreatmentPlan sol){
		int [][] aper_shape, new_aper;

		Random r = new Random();
		int random_intensity,range;
		List<Aperture> list_aperture;
		sol.apertures.clear();
		for(int i=0;i<numAngles;i++){
			aper_shape = initial_aperture_shape.get(i); //initial shape for angle
			list_aperture = new ArrayList<Aperture>();
			for(int z=0;z<5;z++){
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
				list_aperture.add(aper);
			}
			sol.apertures.add(list_aperture);
			sol.aperturesLimit=apertureShapes;
		}
		
	}
	
	//----- hilos
	public void run() {
		
	}
	public String bmltMLC() {
		String exit="";
		for(int i=0;i<solution.aperturesBmlts.size();i++) {
			for(int j=0;j<solution.aperturesBmlts.get(i).length;j++) {
				exit=exit+"param x"+(i+1)+(j+1)+" := ";
				for(int k=0;k<solution.aperturesBmlts.get(i)[j].length;k++) {
					exit=exit+(k+1)+"	"+solution.aperturesBmlts.get(i)[j][k]+"\n";
				}
				exit=exit+";\n";
			}
		}
		return exit;
	}
	public TreatmentPlan MLCHeuristic(Beam[] b,int rounded) throws IOException{
		ExportMLC(b,rounded);
		solution.removeBadAperture(5);
		solution.evaluateSolution();
		System.out.println(" "+solution.singleObjectiveValue+" "+solution.getBOT());
		solution.solverIntensity();
		System.out.println(" "+solution.singleObjectiveValue+" "+solution.getBOT());
		return null;
		
	}
	public TreatmentPlan ExportMLC(Beam[] b,int rounded) throws IOException {
		leafSequencing = new sequentialGenerator(totalbmlt, M, o, dd, beamletAngles, numAperture, rounded, numAngles, solution);
		solution.aperturesLimit=apertureShapes;
		solution.evaluateSolution();
	 	solution.aperturesToIntensityMatrixPerStationwithoutIntensty();
	    solution.ApertureMatrixStationToIntensityMatrixPlan();
		String parameterFile =  "extra.dat";
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(1)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
        }
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param  w := ", bwParametersFile);
        for (int i=0;i<solution.weights.length;i++){
            int j = i+1;
            writeLine(j + " " + solution.weights[i]+ "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param  EUD0 := ", bwParametersFile);
        for (int i=0;i<dd.length;i++){
            int j = i+1;
            writeLine(j + " " + dd[i]+ "\t", bwParametersFile);
        }

        writeLine(";\n", bwParametersFile);
        
        writeLine("param R1:= " + M.OrganVoxels[0]+ ";\n", bwParametersFile);
        writeLine("param R2:= " + M.OrganVoxels[1]+ ";\n", bwParametersFile);
        writeLine("param R3:= " + M.OrganVoxels[2]+ ";\n", bwParametersFile);
        
        writeLine("param  bmlt := ", bwParametersFile);
        for (int i=0;i<b.length;i++){
            int j = i+1;
            writeLine(j + " " + b[i].beamlets+ "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        writeLine("param  bmltAcum := ", bwParametersFile);
        int sum=0;
        for (int i=0;i<b.length;i++){
        	sum+=b[i].beamlets;
            int j = i+1;
            writeLine(j + " " + sum+ "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        bwParametersFile.close();
        
         parameterFile =  "aper.dat";
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(1)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
        }
        bwParametersFile =createBufferedWriter(parameterFile);
        writeLine(bmltMLC(), bwParametersFile);
        bwParametersFile.close();
        
        parameterFile =  "MLC.mod";
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(1)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
        }
        bwParametersFile =createBufferedWriter(parameterFile);
        writeLine("option solver gurobi ;\n", bwParametersFile);
        writeLine("param R1 ;\n", bwParametersFile);
        writeLine("param R2 ;\n", bwParametersFile);
        writeLine("param R3 ;\n", bwParametersFile);
        writeLine("param bmlt{1 .. 5} ;\n", bwParametersFile);
        writeLine("param bmltAcum{1 .. 5} ;\n", bwParametersFile);
        writeLine("var dPTVHD{1 .. R1} ;\n", bwParametersFile);
        writeLine("var dRectum{1 .. R2} ;\n", bwParametersFile);
        writeLine("var dBladder{1 .. R3} ;\n", bwParametersFile);
        writeLine("param ddmPTVHD{1 .. R1, 1 .. bmltAcum[5]} ;\n", bwParametersFile);
        writeLine("param ddmRECTUM{1 .. R2, 1 .. bmltAcum[5]}; ;\n", bwParametersFile);
        writeLine("param ddmBLADDER{1 .. R3, 1 .. bmltAcum[5]}; ;\n", bwParametersFile);
        
        for (int i=0;i<solution.apertures.get(0).size();i++){
            int j = 1;
            writeLine("param x"+j+(i+1)+" {1 .. bmlt["+j+"]} ;\n", bwParametersFile);
        }
        for (int i=0;i<solution.apertures.get(1).size();i++){
            int j = 2;
            writeLine("param x"+j+(i+1)+" {1 .. bmlt["+j+"]} ;\n", bwParametersFile);
        }
        for (int i=0;i<solution.apertures.get(2).size();i++){
            int j = 3;
            writeLine("param x"+j+(i+1)+" {1 .. bmlt["+j+"]} ;\n", bwParametersFile);
        }
        for (int i=0;i<solution.apertures.get(3).size();i++){
            int j = 4;
            writeLine("param x"+j+(i+1)+" {1 .. bmlt["+j+"]} ;\n", bwParametersFile);
        }
        for (int i=0;i<solution.apertures.get(4).size();i++){
            int j = 5;
            writeLine("param x"+j+(i+1)+" {1 .. bmlt["+j+"]} ;\n", bwParametersFile);
        }
        
        writeLine("var intensity1{1 .. "+solution.apertures.get(0).size() +"} >= 0, <=20;\n", bwParametersFile);
        writeLine("var intensity2{1 .. "+solution.apertures.get(1).size() +"} >= 0, <=20;\n", bwParametersFile);
        writeLine("var intensity3{1 .. "+solution.apertures.get(2).size() +"} >= 0, <=20;\n", bwParametersFile);
        writeLine("var intensity4{1 .. "+solution.apertures.get(3).size() +"} >= 0, <=20;\n", bwParametersFile);
        writeLine("var intensity5{1 .. "+solution.apertures.get(4).size() +"} >= 0, <=20;\n", bwParametersFile);
        
        writeLine("param EUD0{1 .. 3} ;\n", bwParametersFile);
        
        writeLine("var v1{1 .. R1} ;\n", bwParametersFile);
        writeLine("var v2{1 .. R2} >= 0 ;\n", bwParametersFile);
        writeLine("var v3{1 .. R3} >= 0 ;\n", bwParametersFile);
        writeLine("param w{1 .. 3} ;\n", bwParametersFile);
        
        writeLine("minimize Total_Cost: (w[1]/R1)*(sum{i in 1 .. R1} v1[i]^2) +(w[2]/R2)*(sum{i in 1 .. R2} v2[i]^2) + (w[3]/R3)*(sum{i in 1 .. R3} v3[i]^2)", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        
        writeLine("s.t.\n", bwParametersFile);
        writeLine("voxelsDeviationPTVHD {i in 1..R1}: v1[i] = EUD0[1] - dPTVHD[i];\n", bwParametersFile);
        writeLine("voxelsDeviationRect {i in 1..R2}:  v2[i] >= dRectum[i] - EUD0[2];\n", bwParametersFile);
        writeLine("voxelsDeviationBladder {i in 1..R3}: v3[i]>= dBladder[i] - EUD0[3] ;\n", bwParametersFile);
        
        
        writeLine("voxelsIrradiationPTVHD {i in 1..R1}: dPTVHD[i] = ", bwParametersFile);
        
        writeLine("sum{j in 1..bmltAcum[1]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(0).size();i++){
            int j = i+1;
            writeLine("intensity1["+j+"]* x1"+j+"[j] * ddmPTVHD[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(0).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        
        writeLine("sum{j in (bmltAcum[1]+1)..bmltAcum[2]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(1).size();i++){
            int j = i+1;
            writeLine("intensity2["+j+"]* x2"+j+"[j-bmltAcum[1]] * ddmPTVHD[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(1).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        writeLine("sum{j in (bmltAcum[2]+1)..bmltAcum[3]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(2).size();i++){
            int j = i+1;
            writeLine("intensity3["+j+"]* x3"+j+"[j-bmltAcum[2]] * ddmPTVHD[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(2).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        writeLine("sum{j in (bmltAcum[3]+1)..bmltAcum[4]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(3).size();i++){
            int j = i+1;
            writeLine("intensity4["+j+"]* x4"+j+"[j-bmltAcum[3]] * ddmPTVHD[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(3).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        writeLine("sum{j in (bmltAcum[4]+1)..bmltAcum[5]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(4).size();i++){
            int j = i+1;
            writeLine("intensity5["+j+"]* x5"+j+"[j-bmltAcum[4]] * ddmPTVHD[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(4).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  ;\n", bwParametersFile);
        writeLine("\n", bwParametersFile);
        
        
        writeLine("voxelsIrradiationRectum {i in 1..R2}: dRectum[i]= ", bwParametersFile);
        
        writeLine("sum{j in 1..bmltAcum[1]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(0).size();i++){
            int j = i+1;
            writeLine("intensity1["+j+"]* x1"+j+"[j] * ddmRECTUM[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(0).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        
        writeLine("sum{j in (bmltAcum[1]+1)..bmltAcum[2]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(1).size();i++){
            int j = i+1;
            writeLine("intensity2["+j+"]* x2"+j+"[j-bmltAcum[1]] * ddmRECTUM[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(1).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        writeLine("sum{j in (bmltAcum[2]+1)..bmltAcum[3]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(2).size();i++){
            int j = i+1;
            writeLine("intensity3["+j+"]* x3"+j+"[j-bmltAcum[2]] * ddmRECTUM[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(2).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        writeLine("sum{j in (bmltAcum[3]+1)..bmltAcum[4]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(3).size();i++){
            int j = i+1;
            writeLine("intensity4["+j+"]* x4"+j+"[j-bmltAcum[3]] * ddmRECTUM[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(3).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        writeLine("sum{j in (bmltAcum[4]+1)..bmltAcum[5]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(4).size();i++){
            int j = i+1;
            writeLine("intensity5["+j+"]* x5"+j+"[j-bmltAcum[4]] * ddmRECTUM[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(4).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  ;\n", bwParametersFile);
        writeLine("\n", bwParametersFile);
        
        
        writeLine("voxelsIrradiationBladder {i in 1..R3}: dBladder[i]= ", bwParametersFile);
        
        writeLine("sum{j in 1..bmltAcum[1]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(0).size();i++){
            int j = i+1;
            writeLine("intensity1["+j+"]* x1"+j+"[j] * ddmBLADDER[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(0).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        
        writeLine("sum{j in (bmltAcum[1]+1)..bmltAcum[2]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(1).size();i++){
            int j = i+1;
            writeLine("intensity2["+j+"]* x2"+j+"[j-bmltAcum[1]] * ddmBLADDER[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(1).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        writeLine("sum{j in (bmltAcum[2]+1)..bmltAcum[3]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(2).size();i++){
            int j = i+1;
            writeLine("intensity3["+j+"]* x3"+j+"[j-bmltAcum[2]] * ddmBLADDER[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(2).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        writeLine("sum{j in (bmltAcum[3]+1)..bmltAcum[4]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(3).size();i++){
            int j = i+1;
            writeLine("intensity4["+j+"]* x4"+j+"[j-bmltAcum[3]] * ddmBLADDER[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(3).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  +\n", bwParametersFile);
        writeLine("sum{j in (bmltAcum[4]+1)..bmltAcum[5]}( ", bwParametersFile);
        for (int i=0;i<solution.apertures.get(4).size();i++){
            int j = i+1;
            writeLine("intensity5["+j+"]* x5"+j+"[j-bmltAcum[4]] * ddmBLADDER[i,j]  ", bwParametersFile);
            if(j!=solution.apertures.get(4).size()) {
            	writeLine("+\n", bwParametersFile);
            }
        }
        writeLine(")  ;\n", bwParametersFile);
        writeLine("\n", bwParametersFile);
        
        
        
        bwParametersFile.close();
        
        parameterFile =  "MLC.run";
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(1)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
        }
        bwParametersFile =createBufferedWriter(parameterFile);
        writeLine("model MLC.mod;\n", bwParametersFile);
        
        writeLine("data "+initialBACs[0]+"_"+initialBACs[1] +"_"+ initialBACs[2]+"_"+initialBACs[3] +"_"+initialBACs[4]+"_DDM_PTVHD.dat;\n", bwParametersFile);
        writeLine("data "+initialBACs[0]+"_"+initialBACs[1] +"_"+ initialBACs[2]+"_"+initialBACs[3] +"_"+initialBACs[4]+"_DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+initialBACs[0]+"_"+initialBACs[1] +"_"+ initialBACs[2]+"_"+initialBACs[3] +"_"+initialBACs[4]+"_DDM_RECTUM.dat;\n", bwParametersFile);
        
        writeLine("data aper.dat;\n", bwParametersFile);
        writeLine("data extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display Total_Cost;\n", bwParametersFile);
        
        writeLine("display intensity1;\n", bwParametersFile);
        
        
        writeLine("display intensity2;\n", bwParametersFile);
        
        
        writeLine("display intensity3;\n", bwParametersFile);
        
        writeLine("display intensity4;\n", bwParametersFile);
        
        writeLine("display intensity5;\n", bwParametersFile);
        
        bwParametersFile.close();
        
    
        //Running Process
        String scriptFile = "MLC.run";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);

        
        BufferedReader br = new BufferedReader(isr);
        String line;
        double [][]intensity=new double [5][];
        int i=0,j=0,flag=0;
        for(int x=0;x<solution.apertures.size();x++) {intensity[x]=new double [solution.apertures.get(x).size()];}
        while ((line = br.readLine()) != null) {
        	if(line.contains("objective")) {
        		System.out.println(line);
        	}
        	if(line.contains("intensity1")) {
        		i=0;
        		j=0;
        		flag=1;
        	}
        	if(line.contains("intensity2")) {
        		i=1;
        		j=0;
        	}
        	if(line.contains("intensity3")) {
        		i=2;
        		j=0;
        	}
        	if(line.contains("intensity4")) {
        		i=3;
        		j=0;
        	}
        	if(line.contains("intensity5")) {
        		i=4;
        		j=0;
        	}
        	
        		
        	if(flag==1 && (!line.contains(":=") && !line.contains(";"))&& !line.isEmpty()) {
        		String[] splited = line.split("\\s+");
        		intensity[i][j]=Double.parseDouble(splited[(splited.length-1)]);  
        		j++;
        	}

        }
        solution.setIntensity(intensity);

        br.close();
        
		return null;
	}
	
	public void prueba() throws IOException {
		String scriptFile = "MLC.run";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);

        
        BufferedReader br = new BufferedReader(isr);
        String line;

        double [][]intensity=new double [5][11];
        int i=0,j=0,flag=0;
        //for(int x=0;x<solution.apertureSize();x++) {intensity[x]=new double [solution.apertures.get(x).size()]}
        while ((line = br.readLine()) != null) {
        	if(line.contains("objective")) {
        		System.out.println(line);
        	}
        	if(line.contains("intensity1")) {
        		i=0;
        		j=0;
        		flag=1;
        	}
        	if(line.contains("intensity2")) {
        		i=1;
        		j=0;
        	}
        	if(line.contains("intensity3")) {
        		i=2;
        		j=0;
        	}
        	if(line.contains("intensity4")) {
        		i=3;
        		j=0;
        	}
        	if(line.contains("intensity5")) {
        		i=4;
        		j=0;
        	}
        	
        		
        	if(flag==1 && (!line.contains(":=") && !line.contains(";"))&& !line.isEmpty()) {
        		String[] splited = line.split("\\s+");
        		intensity[i][j]=Double.parseDouble(splited[(splited.length-1)]);  
        		j++;
        	}

        }
        br.close();
	}
	
    public void writeLine(String l, BufferedWriter bw) throws IOException{
    	
	String row;
	row = l; 
	bw.write(row);
	

    }
    private BufferedWriter createBufferedWriter(String route) throws IOException {
        BufferedWriter bw;
        File auxFile = new File(route);
                if (!auxFile.exists()) {
                    bw = new BufferedWriter(new FileWriter(route));
                }else{
                    bw = new BufferedWriter(new FileWriter(route,true));
                }
        return bw;     
    }
}
