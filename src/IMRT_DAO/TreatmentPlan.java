package IMRT_DAO;

/* LogFunction
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.Vector;

import SingleObj.Gurobi_Solver;
import gurobi.GRBException;

import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;

import java.util.concurrent.ConcurrentHashMap;

import java.util.Map;
/**
 *
 * @author gcab623
 * 
 * @Modificaciones Denisse Katherine Maicholl Mauricio
 */
public class TreatmentPlan extends Thread implements Serializable {
	public int evaluationFunction=0;
	public int evaluationMethod=0;
    public int beams;
    public DDM M;
    Organs[] o;
    int[] beamletAngles;;
    public int beamlets;
    public double slope;
    public double[] intensity;
    public Beam[] selAngles;
    public double[] weights;
    public double[] gEUD;
    public double[] a=new double[] {-10,8,2};
    public double singleObjectiveValue;
    public boolean visited;
	int[][] stateAperture;//Guarda cuantas veces se ha escogido una aperture
	ArrayList<int[]> lowIntensityApertureNotChanged;//Matriz que tiene todas las aperturas con baja intensidad que no se han modificado 
    public double scores[];//katty


    public double epsilonDominance[]; // % of the objective value that is used to relax/stress the dominance condition. It might be positive or negative
    public int[][] beamletCoordID;
    public Vector<int[][]> aperturesLimit;   // # DAO : vector con beam por angulo
    public Vector<double[][]> intesities;       // # DAO : vector de intesidades
    public Vector<List<Aperture>> apertures; // # DAO : vector para cada estaci�n que dentro guarda ua lista de las N aperturas por cada una de ellas
    public int init_intensity;               // # DAO : Intesidad inicial para las aperturas
    public int max_intesity;                 // # DAO : Intesidad m�ima para las aperturas
    public int max_delta;                    // # DAO : Delta m�ximo para la variaci�n de intesidad por apertura
    public int max_iter;                     // # DAO : Cantidad m�xima de iteraciones del algoritmo
    public int max_time;                     // # DAO : tiempo m�ximo de ejecuci�n del algoritmo
    public int seed;                         // # DAO : semilla para partir la busqueda local
    public int step_intensity;               // # DAO : Step size for aperture intensity (2)
    public Vector<double[][]>aperturesBmlts;
    public double[] dd;
    
    public TreatmentPlan(int angles,int numOrgans) {
    	this.beams = angles;
	 	this.selAngles = new Beam[angles];
	 	this.intesities = new Vector<double[][]>();
	 	this.apertures = new Vector<List<Aperture>>();
	 	this.aperturesBmlts = new Vector<double[][]>();
	 	this.aperturesLimit=new Vector<int[][]>();
	 	this.scores=new double[numOrgans];
    }
    

    /************************************************************
	 ************** Constructor sobrecargado para DAO ************
	 *************************************************************/
	 public TreatmentPlan(int init, int max_i, int delta, int iter, int time, int seed, int step, int angles,int numOrgans) {
	 	this.beams = angles;
	 	this.selAngles = new Beam[angles];
	 	this.intesities = new Vector<double[][]>();
	 	this.apertures = new Vector<List<Aperture>>();
	 	this.aperturesBmlts = new Vector<double[][]>();
	 	this.aperturesLimit=new Vector<int[][]>();
	 	this.init_intensity = init;
	 	this.max_intesity = max_i;
	 	this.max_delta = delta;
	 	this.max_iter = iter;
	 	this.max_time = time;
	 	this.seed = seed;
	 	this.step_intensity = step;
	 	this.singleObjectiveValue = 0;
	 	this.scores=new double[numOrgans];
	 	
	 	
	}
	 public TreatmentPlan(int init, int max_i, int delta, int iter, int time, int seed, int step, int angles,int numOrgans,Organs[] o, double[] dd,DDM M,int[] beamletAngles) {
		 	this.beams = angles;
		 	this.selAngles = new Beam[angles];
		 	this.intesities = new Vector<double[][]>();
		 	this.apertures = new Vector<List<Aperture>>();
		 	this.aperturesBmlts = new Vector<double[][]>();
		 	this.aperturesLimit=new Vector<int[][]>();
		 	this.init_intensity = init;
		 	this.max_intesity = max_i;
		 	this.max_delta = delta;
		 	this.max_iter = iter;
		 	this.max_time = time;
		 	this.seed = seed;
		 	this.step_intensity = step;
		 	this.singleObjectiveValue = 0;
		 	this.scores=new double[numOrgans];
			this.o=o;
			this.dd=dd;
			this.M=M;
		 	this.beamletAngles=beamletAngles;
		 	
		 	
		}
	 public TreatmentPlan(TreatmentPlan original) {
		 	this.beams = original.beams;
		 	this.selAngles = new Beam[beams];
		 	this.intesities = new Vector<double[][]>();
		 	this.apertures = new Vector<List<Aperture>>();
		 	this.aperturesBmlts = new Vector<double[][]>();
		 	this.aperturesLimit=new Vector<int[][]>();
		 	this.init_intensity = original.init_intensity;
		 	this.max_intesity = original.max_intesity;
		 	this.max_delta = original.max_delta;
		 	this.max_iter = original.max_iter;
		 	this.max_time = original.max_time;
		 	this.seed = original.seed;
		 	this.step_intensity =original.step_intensity;
		 	this.singleObjectiveValue =0;
		 	this.scores=new double[original.scores.length];
		 	updateSolDAO(original);
	 }
	 
	public void copyMatrix2(Vector<int[][]> Original,Vector<int[][]> copy) {
		copy.clear();
		if(Original.size()>5) {
			System.out.print("error 11");
		}
		
		for(int[][] auxAperture:Original) {
			int [][]aux= new int[auxAperture.length][auxAperture[0].length];
			for(int i=0;i<auxAperture.length;i++) {
				for(int j=0;j<auxAperture[0].length;j++) {
					aux[i][j]= auxAperture[i][j];
				}
			}
			copy.add(aux);
			
			
			
		}
		if(copy.size()>5) {
			System.out.print("error 11");
		}
	
	}
	public void copyMatrix(Vector<int[][]> Original,Vector<int[][]> copy) {
		copy.clear();
		if(Original.size()>5) {
			System.out.print("error 11");
		}
		int t=Original.size();
		int index=0;
		for(int x=0;x<Original.size();x++) {
			
			int [][]aux= new int[Original.get(x).length][Original.get(x)[0].length];
			for(int i=0;i<Original.get(x).length;i++) {
				for(int j=0;j<Original.get(x)[0].length;j++) {
					aux[i][j]= Original.get(x)[i][j];
				}
			}
			index++;
			copy.add(aux);
			if(copy.size()>5) {
				System.out.print("error 11");
			}
			
			
		}
		if(copy.size()>5) {
			System.out.print("error 11");
		}
	
	}
	public double[][] copyMatrix(int[][] Original) {

		  double [][]copy= new double[Original.length][Original[0].length];
			for(int i=0;i<Original.length;i++) {
				for(int j=0;j<Original[0].length;j++) {
					copy[i][j]= Original[i][j];
				}
			}
			return copy;
		
	
	}
	public int beamletId(int angle, int leaf, int positionChanged){

    	int index=selAngles[angle].index;

		double[][] beamletsCoord = selAngles[angle].beamletsCoord;
		for(int j=0;j<selAngles[angle].beamlets;j++){
			
			int xcoord=(int)beamletsCoord[j][1]-1;
			int ycoord=(int)beamletsCoord[j][2]-1;
			if(xcoord==leaf&&ycoord==positionChanged) {
				break;
			}
			else
				index=index+1;
		}
    	
		
		return index-1;
		
	}
	public void aperturesToIntensityMatrixPerStation(){
    	Vector<List<Aperture>> vectorStations = apertures;
    	//Vector<int[][]> vectorIntensities = tp.intesities;
    	intesities.clear();
    	//recorrido de estaciones
    	for(int i=0; i<vectorStations.size(); i++){
    		List<Aperture> listApertures = vectorStations.elementAt(i);  
    		double[][] intensitiesMatrix = copyMatrix(aperturesLimit.get(i));
    		
    		//recorrido de lista de aperturas por estaci�n
            for(int j=0; j<listApertures.size(); j++ ){
            	int[][] aperture = listApertures.get(j).aperture;
            	double intensity = listApertures.get(j).intensity;
            	
            	//recorrido de fila de matriz de intensidades
            	for(int row=0; row<intensitiesMatrix.length; row++){
            		if(row==11 && aperturesLimit.size()>5){
            			int flag=11;
            			System.out.print("error"+ flag);
            		}
            		int leftLimit = aperture[row][0];
            		int rightLimit = aperture[row][1];
            		
            		//recorrido de columnas de matrix asinando intensidades seg�n limites
            		for(int column=leftLimit; (column<=rightLimit)&&(leftLimit!=-1); column++){
                		
            			intensitiesMatrix[row][column] =intensitiesMatrix[row][column] + intensity;
            			
                	}
            		
            	}
           
            }
            intesities.add(intensitiesMatrix);
        }
    	
    }
	public void aperturesToIntensityMatrixPerStationwithoutIntensty(){
    	Vector<List<Aperture>> vectorStations = apertures;
    	//Vector<int[][]> vectorIntensities = tp.intesities;
    	intesities.clear();
    	//recorrido de estaciones
    	for(int i=0; i<vectorStations.size(); i++){
    		List<Aperture> listApertures = vectorStations.elementAt(i);  
    		
    		
    		//recorrido de lista de aperturas por estaci�n
            for(int j=0; j<listApertures.size(); j++ ){
            	int[][] aperture = listApertures.get(j).aperture;
            	double[][] intensitiesMatrix = copyMatrix(aperturesLimit.get(i));
            	
            	//recorrido de fila de matriz de intensidades
            	for(int row=0; row<intensitiesMatrix.length; row++){
            		int leftLimit = aperture[row][0];
            		int rightLimit = aperture[row][1];
            		
            		//recorrido de columnas de matrix asinando intensidades seg�n limites
            		for(int column=leftLimit; (column<=rightLimit)&&(leftLimit!=-1); column++){
                		
            			intensitiesMatrix[row][column] =intensitiesMatrix[row][column] + 1;
            			
                	}
            		
            	}
            	intesities.add(intensitiesMatrix);
            }
            
        }
    	
    }
	public void changeAperture(int[][] newAperture,int angle, int positionAperture) {
		apertures.get(angle).get(positionAperture).aperture=newAperture;
		//apertures.get(angle).add(newAperture);
		
	}
	
	
	//---- Point Crossover function
	
	
	/**
	 * Cruza la aperturas en un punto
	 * @param angleChoose
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @return
	 */
	public ArrayList<int[][]>  singlePointCrossover(int angleChoose,int idApertureOne,int idApertureTwo) {
		Random r = new Random(System.currentTimeMillis());
		Aperture father0ne;
		Aperture fatherTwo;
		ArrayList<int[][]> childs=new ArrayList<int[][]>();
		int size=0;
		father0ne=apertures.get(angleChoose).get(idApertureOne);
		fatherTwo=apertures.get(angleChoose).get(idApertureTwo);
		size=father0ne.aperture.length;
		int pointCrossover= r.nextInt(size);
		int[][] childOne=onecrossApertures(father0ne,fatherTwo,pointCrossover,size);
		int[][] childTwo=onecrossApertures(fatherTwo,father0ne,pointCrossover,size);
		childs.add(childOne);
		childs.add(childTwo);
		return childs;
	}
	/**
	 * Cruza la aperturas en un punto para las tres representaciones de cotrutz
	 * @param angleChoose
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @return
	 */
	public ArrayList<int[][]>  singlePointCrossover2(int angleChoose,int idApertureOne,int idApertureTwo) {
		Random r = new Random(System.currentTimeMillis());
		Aperture father0ne;
		Aperture fatherTwo;
		ArrayList<int[][]> childs=new ArrayList<int[][]>();
		int size=0;
		father0ne=apertures.get(angleChoose).get(idApertureOne);
		fatherTwo=apertures.get(angleChoose).get(idApertureTwo);
		size=father0ne.aperture.length;
		int pointCrossover= r.nextInt(size);
		int[][] childOne=onecrossApertures(father0ne,fatherTwo,pointCrossover,size);
		int[][] childTwo=onecrossApertures(fatherTwo,father0ne,pointCrossover,size);
		childs.add(childOne);
		childs.add(childTwo);
		return childs;
	}
	
	/**
	 * FALTA AGREGAR EL SEGUNDO PUNTO
	 * @param angleChoose
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @return
	 */
	public ArrayList<int[][]>  twoPointCrossover(int angleChoose,int idApertureOne,int idApertureTwo) {
		Random r = new Random(System.currentTimeMillis());
		Aperture father0ne;
		Aperture fatherTwo;
		ArrayList<int[][]> childs=new ArrayList<int[][]>();
		int size=0;
		father0ne=apertures.get(angleChoose).get(idApertureOne);
		fatherTwo=apertures.get(angleChoose).get(idApertureTwo);
		size=father0ne.aperture.length;
		int pointCrossover= r.nextInt(size);
		int[][] childOne=twoCrossApertures(father0ne,fatherTwo,pointCrossover,pointCrossover,size);
		int[][] childTwo=twoCrossApertures(fatherTwo,father0ne,pointCrossover,pointCrossover,size);
		childs.add(childOne);
		childs.add(childTwo);
		return childs;
	}
	
	/**
	 * por hacer
	 * @param angleChoose
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @return
	 */
	public ArrayList<int[][]>  uniformPointCrossover(int angleChoose,int idApertureOne,int idApertureTwo) {
		Random r = new Random(System.currentTimeMillis());
		Aperture father0ne;
		Aperture fatherTwo;
		ArrayList<int[][]> childs=new ArrayList<int[][]>();
		int size=0;
		father0ne=apertures.get(angleChoose).get(idApertureOne);
		fatherTwo=apertures.get(angleChoose).get(idApertureTwo);
		size=father0ne.aperture.length;
		int pointCrossover= r.nextInt(size);
		int[][] childOne=uniformCrossApertures(father0ne,fatherTwo,size);
		int[][] childTwo=uniformCrossApertures(fatherTwo,father0ne,size);
		childs.add(childOne);
		childs.add(childTwo);
		return childs;
	}
	
	/**
	 * This function crossover to apertures  from the solution and generate a new aperture. Whe use single point crossover
	 * the point crossover is chose ramdonly
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @param pointCrossover
	 */
	public int[][] onecrossApertures(Aperture father0ne,Aperture fatherTwo,int pointCrossover, int size) {
		int[][] child =  new int[father0ne.aperture.length][father0ne.aperture[0].length]; ;
		
		
		for(int i=0;i<pointCrossover;i++) {
			child[i]=father0ne.aperture[i].clone();
		}
		
		for(int i=pointCrossover;i<size;i++) {
			child[i]=fatherTwo.aperture[i].clone();
		}
		return child;

		
		
	}
	/**
	 * This function crossover to apertures on the left side  from the solution and generate a new aperture. Whe use single point crossover
	 * the point crossover is chose ramdonly
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @param pointCrossover
	 */
	public int[][] crossoverApertureLeft(Aperture father0ne,Aperture fatherTwo,int pointCrossover, int size) {
		int[][] child =  new int[father0ne.aperture.length][father0ne.aperture[0].length]; ;
		
		
		for(int i=0;i<pointCrossover;i++) {
			child[i]=father0ne.aperture[i].clone();
		}
		
		for(int i=pointCrossover;i<size;i++) {
			child[i]=fatherTwo.aperture[i].clone();
		}
		return child;

		
		
	}
	/**
	 * This function crossover to apertures on the rigth side  from the solution and generate a new aperture. Whe use single point crossover
	 * the point crossover is chose ramdonly
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @param pointCrossover
	 */
	public int[][] crossoverApertureRigth(Aperture father0ne,Aperture fatherTwo,int pointCrossover, int size) {
		int[][] child =  new int[father0ne.aperture.length][father0ne.aperture[0].length]; ;
		
		
		for(int i=0;i<pointCrossover;i++) {
			child[i]=father0ne.aperture[i].clone();
		}
		
		for(int i=pointCrossover;i<size;i++) {
			child[i]=fatherTwo.aperture[i].clone();
		}
		return child;

		
		
	}
	
	/**
	 * This function crossover to apertures on the rigth side  from the solution and generate a new aperture. Whe use single point crossover
	 * the point crossover is chose ramdonly
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @param pointCrossover
	 */
	public int[][] crossoverIntensity(Aperture father0ne,Aperture fatherTwo,int pointCrossover, int size) {
		int[][] child =  new int[father0ne.aperture.length][father0ne.aperture[0].length]; ;
		
		
		Random rand = new Random();

        /*Generar un n�mero aleatorio entre 0 y 5 (ambos incluidos)*/
        int aperture= rand.nextInt((apertures.get(0).size()+1));
		return child;

		
		
	}
	
	/**
	 * This function crossover to apertures  from the solution and generate a new aperture. Whe use two point crossover
	 * the point crossover is chose ramdonly
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @param pointCrossover
	 */
	public int[][] twoCrossApertures(Aperture father0ne,Aperture fatherTwo,int firstPointCrossover,int secondPointCrossover, int size) {
		int[][] child =  new int[father0ne.aperture.length][father0ne.aperture[0].length]; ;
		
		
		for(int i=0;i<firstPointCrossover;i++) {
			child[i]=father0ne.aperture[i].clone();
		}
		
		for(int i=firstPointCrossover;i<secondPointCrossover;i++) {
			child[i]=fatherTwo.aperture[i].clone();
		}
		for(int i=secondPointCrossover;i<size;i++) {
			child[i]=father0ne.aperture[i].clone();
		}
		return child;

		
		
	}
	
	/**
	 * This function crossover to apertures  from the solution and generate a new aperture. Whe use two point crossover
	 * the point crossover is chose ramdonly
	 * @param idApertureOne
	 * @param idApertureTwo
	 * @param pointCrossover
	 */
	public int[][] uniformCrossApertures(Aperture father0ne,Aperture fatherTwo, int size) {
		int[][] child =  new int[father0ne.aperture.length][father0ne.aperture[0].length]; ;
		Random r = new Random(System.currentTimeMillis());
		int cross= r.nextInt(2);//
		for(int i=0;i<size;i++) {
			cross= r.nextInt(2);
			if(cross==1)
				child[i]=father0ne.aperture[i].clone();
			else
				child[i]=fatherTwo.aperture[i].clone();
		}
		
		
		return child;

		
		
	}
	
	
	
    public int[] intensityMatrixStationToIntensityMatrixPlan(Vector<int[][]> stationIntensities){
    	int numColumns = (stationIntensities.size())*(stationIntensities.elementAt(0).length)*(stationIntensities.elementAt(0)[0].length);
    	int [] solutionVector = new int[numColumns];
    	int columns = 0;
    	
    	//recorrido de estaciones
    	for(int i=0; i<stationIntensities.size(); i++){
    		int[][] intensityAperture = stationIntensities.elementAt(i);
    		
    		for(int row=0; row < intensityAperture.length; row++){
    			for(int column=0; column < intensityAperture[row].length; column++){
    				solutionVector[columns] = intensityAperture[row][column];
        		}
    		}
    	}
    	return solutionVector;
    }
    public int[] bealetsForAngle() {
    	int[] aux=new int[selAngles.length];
    	for(int j=0;j<aux.length;j++){
			aux[j]=selAngles[j].beamlets;
		}
    	
    	return  aux;
    }
    
    
    public int apertureSize() {
    	int size=apertures.size();
    	size=size*apertures.get(0).size();
    	return size;
    }
    
    /**
     * Funcion que crea una estructura que realciona el beamlet con - angulo,hoja,posicion
     * @return
     */
	public int[][] beamletCoord(){
	    	
	    	int[][] BeamletVector = new int[beamlets][3] ;//contiene angulo, coordx y coord y del beamlet, la posicion esta realcionada al beamlet
	    	//int[][][] Beamletcube=new int[selAngles.length][][];
	    	int index=0;
	    	for(int i=0;i<selAngles.length;i++){

				double[][] beamletsCoord = selAngles[i].beamletsCoord;
				for(int j=0;j<selAngles[i].beamlets;j++){
					int xcoord=(int)beamletsCoord[j][1]-1;
					int ycoord=(int)beamletsCoord[j][2]-1;

					BeamletVector[index][0]=i;
					BeamletVector[index][1]=xcoord;
					BeamletVector[index][2]=ycoord;
					
					index=index+1;
				}
	    	}
	    	//recorrido de estaciones
	    	beamletCoordID=BeamletVector;
			return BeamletVector;
	
	    	
	    }
	public int searchBeamletID(int angle, int leaf, int positionChanged) {
		int beamlet=positionChanged; // original era 0
		for(int i=0;i<beamletCoordID.length;i++) {
			if(beamletCoordID[i][0]==angle &&beamletCoordID[i][1]==leaf &&beamletCoordID[i][2]==positionChanged) {
				beamlet=i;
				return beamlet;
			}
		}
		
		
		return beamlet;
		
	}
    
    
    public void intensityMatrixStationToIntensityMatrixPlan(){
    	
    	intensity=new double[beamlets] ;
    	int index=0;
    	for(int i=0;i<intesities.size();i++){
    		double[][] intensities;
    		intensities=intesities.get(i);
			double[][] beamletsCoord = selAngles[i].beamletsCoord;
			for(int j=0;j<selAngles[i].beamlets;j++){
				int xcoord=(int)beamletsCoord[j][1]-1;
				int ycoord=(int)beamletsCoord[j][2]-1;
				//System.out.print(xcoord+" ");
				//System.out.println(ycoord);
				//System.out.println(beamlets);
				//System.out.println(intensities.length);
				intensity[index]=intensities[xcoord][ycoord];
				index=index+1;
			}
    	}
    	//recorrido de estaciones

    	
    }
    public void IntensityMatrixPlan(){
    	
    	
    	
    	
    	intensity=new double[beamlets] ;
    	for(int t=0;t<beamlets;t++) {intensity[t]=0;}
    	int index=0;
    	for(int a=0;a<aperturesBmlts.size();a++) {
	    	double[][] intensities;
	    	
	    	for(int i=0;i<aperturesBmlts.get(a).length;i++){
	    		intensities=aperturesBmlts.get(a);
				
				for(int j=0;j<selAngles[a].beamlets;j++){
					
					//System.out.print(xcoord+" ");
					//System.out.println(ycoord);
					//System.out.println(beamlets);
					//System.out.print((j+index)+"-"+i+"-"+j+" ");
					intensity[(j+index)]=intensity[(j+index)]+intensities[i][j]*apertures.get(a).get(i).intensity;
					
				}
				//System.out.println("");
	    	}
	    	//System.out.println("");
	    	index=index+selAngles[a].beamlets;
	    	//System.out.println(index);
    	}
    	

    	
    }
    public void ApertureMatrixStationToIntensityMatrixPlan(){
    	
    	aperturesBmlts.clear();
    	int index=0;
    	for(int x=0;x<beams;x++) {
    		double[][] newAperture = new double[apertures.get(x).size()][selAngles[x].beamlets];
    		
	    	for(int i=0;i<apertures.get(x).size();i++){
	    		
	    		double[][] intensities;
	    		intensities=intesities.get(index);
	    		
				double[][] beamletsCoord = selAngles[x].beamletsCoord;
				for(int j=0;j<selAngles[x].beamlets;j++){
					int xcoord=(int)beamletsCoord[j][1]-1;
					int ycoord=(int)beamletsCoord[j][2]-1;
					//System.out.print(j+"-"+xcoord+"-"+ycoord+" ");
					newAperture[i][j]=intensities[xcoord][ycoord];
					
				}
				//System.out.println("");
				index=index+1;
	    	}
	    	aperturesBmlts.add(newAperture);
	    	
    	}
    	//recorrido de estaciones

    	
    }
    public void setWeights(double weight[]) {
    	this.weights=weight.clone();
    }
    
 
    public void updateSolDAO(TreatmentPlan s) {
  	
    	
    	
        this.beamlets = s.beamlets;
        this.beams = s.beams;
        this.selAngles = new Beam[this.beams];
        this.singleObjectiveValue = s.singleObjectiveValue;
        this.init_intensity = s.init_intensity;
	 	this.max_intesity = s.max_intesity;
	 	this.max_delta = s.max_delta;
	 	this.max_iter = s.max_iter;
	 	this.max_time = s.max_time;
	 	this.seed = s.seed;
	 	this.step_intensity = s.step_intensity;
	 	this.weights=s.weights.clone();
	 	this.dd=s.dd;
        this.o=s.o;
        this.M=s.M;
        this.intensity = s.intensity.clone();
        this.beamletCoordID=s.beamletCoordID;
	 	this.intesities = new Vector<double[][]>();
	 	this.beamletAngles=s.beamletAngles;
	 	this.apertures = new Vector<List<Aperture>>();
	 

	 	copyMatrix( s.aperturesLimit, this.aperturesLimit);

	 	System.arraycopy(s.selAngles, 0, this.selAngles, 0, this.beams);
	 	//System.arraycopy(s.intensity,0,this.intensity,0,this.intensity.length);
	 	
	 	for (int i = 0; i < s.apertures.size(); i++) {
	 		List<Aperture> nuevo =  new ArrayList<Aperture>();
	 		for(int j=0;j<s.apertures.get(i).size();j++) {
	 			Aperture nuevaApertura=new Aperture(s.apertures.get(i).get(j));
	 			nuevo.add(nuevaApertura);
	 		}
	 		apertures.add(nuevo);
	 	}
	 	
	 	System.arraycopy(s.scores, 0, this.scores, 0, s.scores.length);
        
	 	

    }
    
    public void updateApertureDAO(TreatmentPlan s) {
  	

	 	System.arraycopy(s.selAngles, 0, this.selAngles, 0, this.beams);
	 	//System.arraycopy(s.intensity,0,this.intensity,0,this.intensity.length);
	 	
	 	for (int i = 0; i < s.apertures.size(); i++) {
	 		List<Aperture> nuevo =  new ArrayList<Aperture>();
	 		for(int j=0;j<s.apertures.get(i).size();j++) {
	 			Aperture nuevaApertura=new Aperture(s.apertures.get(i).get(j));
	 			nuevo.add(nuevaApertura);
	 		}
	 		apertures.add(nuevo);
	 	}

    }
    
    /**
     * Cambia el valos de intensidades de todas las aperturas dada una matriz
     * @param newIntensity
     */
    public void setIntensity(double[][] newIntensity) {
		for(int i=0; i<apertures.size(); i++){
			
			List<Aperture> list_apertures = apertures.get(i);
			//recorrido de aperturas por angulo
			for(int j=0; j<list_apertures.size(); j++){
				//list_apertures.get(j).intensity=Math.round(newIntensity[i][j]);
				list_apertures.get(j).intensity=newIntensity[i][j];
				
			}
			
		}
    }
    
    /**
     * Cambia el valor de una intensidad especifica en el angulo i , en la apertura j con el valor value
     * @param i
     * @param j
     * @param value
     */
    public void setIntensity(int i, int j, double value) {
		apertures.get(i).get(j).intensity=value;

    }
    public void getIntensity() {
    	for(int i=0; i<apertures.size(); i++){
			
				List<Aperture> list_apertures = apertures.get(i);
				//recorrido de aperturas por angulo
				for(int j=0; j<list_apertures.size(); j++){
					//list_apertures.get(j).intensity=Math.round(newIntensity[i][j]);
					System.out.print(list_apertures.get(j).intensity+" ");
					
				}
				System.out.println();
			}
	    }
    public double[][] getIntensityApertures() {
    	double[][]intensiyMatrix=new double[apertures.size()][];
    	for(int i=0; i<apertures.size(); i++){
    		
    		
			List<Aperture> list_apertures = apertures.get(i);
			
			intensiyMatrix[i]=new double[list_apertures.size()];
			//recorrido de aperturas por angulo
			for(int j=0; j<list_apertures.size(); j++){
				//list_apertures.get(j).intensity=Math.round(newIntensity[i][j]);
				intensiyMatrix[i][j]=list_apertures.get(j).intensity;
				
			}
		}
    	return intensiyMatrix;
    }
    public double getBOT() {
    	double total=0;
    	for(int i=0; i<apertures.size(); i++){
			
				List<Aperture> list_apertures = apertures.get(i);
				//recorrido de aperturas por angulo
				for(int j=0; j<list_apertures.size(); j++){
					//list_apertures.get(j).intensity=Math.round(newIntensity[i][j]);
					total=total+list_apertures.get(j).intensity;
					
				}
				
			}
		return total;
	    }
    /**
     * Devuelve la cantidad de aperturas bajo el limite que tiene las aperturas
     * @return
     */
    public int getLowIntensityperBeam() {
    	int totalLowIntensity=0;
    	for(int i=0; i<apertures.size(); i++){
			
				List<Aperture> list_apertures = apertures.get(i);
				//recorrido de aperturas por angulo
				int numberLow=0;
				for(int j=0; j<list_apertures.size(); j++){
					if(list_apertures.get(j).intensity<1){
						numberLow=numberLow+1;
					}
					
					
				}
				//System.out.println(numberLow);
				totalLowIntensity=totalLowIntensity+numberLow;
			}
    	return totalLowIntensity;
	}
    
    public double getTotalIntensityperAngle(int angle) {
    	double totalIntensity=0;
    	
			
		List<Aperture> list_apertures = apertures.get(angle);
		//recorrido de aperturas por angulo
		int numberLow=0;
		for(int j=0; j<list_apertures.size(); j++){
			
			totalIntensity=totalIntensity+list_apertures.get(j).intensity;
			
			
			
		}
				
			
    	return totalIntensity;
	}
    
    /**
     * Devuelve cuantos angulos tienen aperturas con intensidad bbajo el limite
     * @return
     */
    public int getLowIntensityAngle() {
    	
    	int angleLow=0;
    	for(int i=0; i<apertures.size(); i++){
			
				List<Aperture> list_apertures = apertures.get(i);
				//recorrido de aperturas por angulo
				int numberLow=0;
				for(int j=0; j<list_apertures.size(); j++){
					if(list_apertures.get(j).intensity<1){
						numberLow=numberLow+1;
					}
					
					
				}
				if(numberLow==5)
					angleLow=angleLow+1;
				//System.out.println(numberLow);
				
			}
    	return angleLow;
	}
    
    /**
     * eliminas las aperturas que tengaun un valor menor a una intensidad dada de un beam especifico
     * @return la cantidad de aperturas eliminadas
     */
    public int removeBadAperturePerBeam(int beamPosition, int minimunIntensity) {
    	
		int numberLow=0;
		List<Aperture> list_apertures = apertures.get(beamPosition);
		List<Aperture> list_aperturesToRemove=new ArrayList<Aperture>();
		//recorrido de aperturas por angulo

		for(int j=0; j<list_apertures.size(); j++){
			if(list_apertures.get(j).intensity<minimunIntensity){
				list_aperturesToRemove.add(list_apertures.get(j));
				numberLow=numberLow+1;
			}
		}
		list_apertures.removeAll(list_aperturesToRemove);
		
		
    	return numberLow ; 
    }
    
    /**
     * eliminas las aperturas que tengaun un valor menor a una intensidad dada de un beam especifico
     * @return la cantidad de aperturas eliminadas
     */
    public int removeBadAperture( int numAperPerBeam) {
    	
		int numberLow=0;

		//recorrido de aperturas por angulo
		for(int i=0;i<apertures.size();i++) {
			
			List<Aperture> list_apertures = apertures.get(i);
			 Collections.sort(list_apertures,Aperture.intensityComparator);
			List<Aperture> list_aperturesToRemove=new ArrayList<Aperture>();
			for(int j=0; j<list_apertures.size(); j++){
				if(j>=(numAperPerBeam)){
					list_aperturesToRemove.add(list_apertures.get(j));
					numberLow=numberLow+1;
				}
			}
			list_apertures.removeAll(list_aperturesToRemove);
		}

		
		
    	return numberLow ; 
    }


    
    /**
     * Reemplaza las intensidades con valores menor a un limite especificado
     */
    public void setBadIntensityPerAperture(int limitIntensity, int newIntensity) {
    	for(int i=0; i<apertures.size(); i++){
		
			List<Aperture> list_apertures = apertures.get(i);
			//recorrido de aperturas por angulo
	
			for(int j=0; j<list_apertures.size(); j++){
				if(list_apertures.get(j).intensity<limitIntensity){
					list_apertures.get(j).intensity=newIntensity;
					
				}
			}
    	}
		

    }
    /**
     * Elimina una apertura especifica  de un angulo especifico
     * 
     */
    public void removeAperture(int beamPosition, int indexAperture) {
    	
		
		apertures.get(beamPosition).remove(indexAperture);

    }
    /**
     * agregar una apertura especifica  de un angulo especifico
     * 
     */
    public void addAperture(List<Aperture> list_aperture) {
    	
		
    	apertures.add(list_aperture);

    }
    /**
     * Elimina una apertura especifica  de un angulo especifico
     * 
     */
    public void removeAperture(int beamPosition, Vector<double[]> v) {
    	List<Aperture> list_apertures = apertures.get(beamPosition);
    	List<Aperture> list_aperturesToRemove=new ArrayList<Aperture>();
    	list_aperturesToRemove.add(apertures.get(beamPosition).get((int) v.get(0)[0]));
    	list_aperturesToRemove.add(apertures.get(beamPosition).get((int) v.get(1)[0]));
		list_apertures.removeAll(list_aperturesToRemove);

    }
    
    
    /**
     * Obtienen la intensidad maxima del beam
     * @param beamPosition
     * @return
     */
    public double getMaxIntensityBeam(int beamPosition) {
    
				List<Aperture> list_apertures = apertures.get(beamPosition);
				//recorrido de aperturas por angulo
				double maxIntensity=0;
				for(int j=0; j<list_apertures.size(); j++){
					if(list_apertures.get(j).intensity>maxIntensity){
						maxIntensity=list_apertures.get(j).intensity;
					}
					
					
				}
				return maxIntensity;
	    }
    /**
     * Obtiene el segmentQualityFactor de una apertura en esfecifico
     * @param segmentWeight
     * @param intensityWeight
     * @param beamPosition
     * @param aperturePosition
     * @param beamletAngle
     * @return
     */
    public double segmentQualityFactor(float segmentWeight,float intensityWeight, int beamPosition, int aperturePosition, int[] beamletAngle) {
    	 double SQF=0; 
    	 int segmentith=apertures.get(beamPosition).get(aperturePosition).getTotalApertureOpen();
    	 int segmentTotal=beamletAngle[beamPosition];
    	 double weightith=apertures.get(beamPosition).get(aperturePosition).intensity;
    	 double weightMax=getMaxIntensityBeam(beamPosition);
    	 
    	 SQF=(segmentWeight*(segmentith/segmentTotal))+(intensityWeight*(weightith/weightMax));
    	 return SQF;
    }
    public double segmentQualityFactor(float segmentWeight,float intensityWeight, int beamPosition, int aperturePosition, int segmentTotal,double weightMax) {
    	double SQF=0; 
    	int segmentith=apertures.get(beamPosition).get(aperturePosition).getTotalApertureOpen();
    	double weightith=apertures.get(beamPosition).get(aperturePosition).intensity;

    	SQF=(segmentWeight*(segmentith/segmentTotal))+(intensityWeight*(weightith/weightMax));
    	return SQF;
   }
    
    /**
     * 
     * @param segmentWeight
     * @param intensityWeight
     * @param beamPosition
     * @param aperturePosition
     * @param segmentTotal
     * @param sizeBeam
     * @return
     */
    public double[] SQFperBeam(float segmentWeight,float intensityWeight, int beamPosition,int segmentTotal ,int sizeBeam) {
     double[] SQFBeam=new double[sizeBeam];
     double SQF=0;   
   	 double weightMax=getMaxIntensityBeam(beamPosition);
   	System.out.println();
   	for(int j=0; j<sizeBeam; j++){
   		SQF=segmentQualityFactor(segmentWeight, intensityWeight, beamPosition, j, segmentTotal, weightMax);
		//list_apertures.get(j).intensity=newIntensity[i][j];
   		SQFBeam[j]=SQF;
   		System.out.print(SQF+" ");
		
	}
  
   	 return SQFBeam;
   }
    
    
    /**
     * 
     * @param segmentWeight
     * @param intensityWeight
     * @param beamPosition
     * @param aperturePosition
     * @param segmentTotal
     * @param sizeBeam
     * @return
     */
    public Vector<double[]> indexSQFperBeam(float segmentWeight,float intensityWeight, int beamPosition,int beamletAngle[]) {
	    Vector<double[]> v = new Vector<double[]>();
	    int segmentTotal=0;
	    segmentTotal=beamletAngle[beamPosition];
	    double SQF=0;   
	   	double weightMax=getMaxIntensityBeam(beamPosition);
	   
	   	List<Aperture> list_apertures = apertures.get(beamPosition);
	   	for(int j=0; j<list_apertures.size(); j++){
	   		double[] SQFBeam=new double[2];
	   		SQF=segmentQualityFactor(segmentWeight, intensityWeight, beamPosition, j, segmentTotal, weightMax);
			//list_apertures.get(j).intensity=newIntensity[i][j];
	   		SQFBeam[0]=j;
	   		SQFBeam[1]=SQF;
	   		v.add(SQFBeam);
	   		
			
		}
  
   	 return v;
   }
    public void solutionSQF(float segmentWeight,float intensityWeight, int []beamletAngle) {
	 int segmentTotal=0;

	for(int i=0; i<apertures.size(); i++){
		
		List<Aperture> list_apertures = apertures.get(i);
		list_apertures.size();
		//recorrido de aperturas por angulo

		segmentTotal=beamletAngle[i];
		SQFperBeam(segmentWeight, intensityWeight, i, segmentTotal, list_apertures.size());
		//list_apertures.get(j).intensity=newIntensity[i][j];
		//System.out.print(b);
		
		
		
	}
   	 
   	 //System.out.print(b);
   }
    
    
    
    
    
    public void setIntensityRounded(double[][] newIntensity) {
		for(int i=0; i<apertures.size(); i++){
			
			List<Aperture> list_apertures = apertures.get(i);
			//recorrido de aperturas por angulo
			for(int j=0; j<list_apertures.size(); j++){
				list_apertures.get(j).intensity=Math.round(newIntensity[i][j]);
				//list_apertures.get(j).intensity=newIntensity[i][j];
				
			}
			
		}
    }
    public void setOF(double newObjectiveValue) {
		this.singleObjectiveValue=newObjectiveValue;
    }
    

    public double truncateNumber(double n, int decimalplace)   
    {   
    //moves the decimal to the right   
    	n = n* Math.pow(10, decimalplace);   
    //determines the floor value  
    	n = Math.floor(n);   
    //dividing the floor value by 10 to the power decimalplace  
    	n = n / Math.pow(10, decimalplace);   
    //prints the number after truncation  
    	return n;
    }   
    public void print() {
    	for(double i:scores) {
    		System.out.print(i+",");
    	}
    	System.out.print(" ");
    	System.out.println(singleObjectiveValue);
    	
    }
	public void printValues() {
		for (int i = 0; i < selAngles.length; i++) {
			System.out.print(selAngles[i].index+" ");
			
		}
		System.out.println(scores[1]+" "+scores[2]);
		
		
	}
	public int[] getBAC() {
		int[] auxBAC = new int[selAngles.length];
		for (int i = 0; i < selAngles.length; i++) {
			auxBAC[i] = selAngles[i].index;
			
		}
		return auxBAC;
		// TODO Auto-generated method stub
		
	}
     public boolean  isSameValue(TreatmentPlan q) {
    	for(int i=0;i<scores.length;i++) {
    		if(scores[i]!=q.scores[i]) {
    			return false;
    		}
    	}
    	 
		return true;
    	 
     }
	 public boolean isDominatedTwoObjective(TreatmentPlan q) {
		double PTVp,OARp,PTVq,OARq;
		PTVp=scores[0];
		PTVq=q.scores[0];
		OARp=scores[1]+scores[2];
		OARq=q.scores[1]+q.scores[2];
		PTVp=truncateNumber(PTVp,2);
		PTVq=truncateNumber(PTVq,2);
		OARp=truncateNumber(OARp,2);
		OARq=truncateNumber(OARq,2);
		
		
		if((PTVp>PTVq && OARp>=OARq) || (PTVp>=PTVq && OARp>OARq)) {
		    return true;
		}
		
		return false;

	 }
	 public boolean isDominatedPtvRectum(TreatmentPlan q) {
		double PTVp,OARp,PTVq,OARq;
		PTVp=scores[0];
		PTVq=q.scores[0];
		OARp=scores[1];
		OARq=q.scores[1];
		
		if((PTVp>PTVq && OARp>=OARq) || (PTVp>=PTVq && OARp>OARq)) {
		    return true;
		}
		
		return false;

	 }
    public boolean isDominated(TreatmentPlan s) {
        boolean isDominated = true; //means, s dominates 'this'
        //We assume that first index in gEDU  corresponds to the target
        //which is equal for all the solutions
        for (int i = 1; i < this.gEUD.length - 1; i++) {
            if (this.gEUD[i] < s.gEUD[i] * (1 - s.epsilonDominance[i]) ) {
                isDominated = false;
                break;
            }
        }
        return (isDominated);
    }
    public boolean isDominatedEvaluetion(TreatmentPlan s) {
        boolean isDominated = true; //means, s dominates 'this'
        //We assume that first index in gEDU  corresponds to the target
        //which is equal for all the solutions
        for (int i = 0; i < this.scores.length; i++) {
            if (this.scores[i]  < s.scores[i] ) {
                isDominated = false;
                break;
            }
        }
        return (isDominated);
    }
    public boolean isDominatedEvaluetionORPen(TreatmentPlan s) {
        boolean isDominated = true; //means, s dominates 'this'
        //We assume that first index in gEDU  corresponds to the target
        //which is equal for all the solutions
        for (int i = 1; i < this.scores.length; i++) {
            if (this.scores[i]  < s.scores[i] ) {
                isDominated = false;
                break;
            }
        }
        return (isDominated);
    }
    
    
	public void roundIntensity() {
		for(int i=0;i<apertures.size();++i) {
			for(int j=0;j<apertures.get(i).size();++j) {
				apertures.get(i).get(j).intensity=Math.round(apertures.get(i).get(j).intensity);
			}
		}
		
		// TODO Auto-generated method stub
		
	}
	public void printIntensity() {
		for(int i=0;i<apertures.size();++i) {
			System.out.println("angulo: "+i);
			for(int j=0;j<apertures.get(i).size();++j) {
				System.out.print(apertures.get(i).get(j).intensity+" ");
			}
			System.out.println();
		}
		
		// TODO Auto-generated method stub
		
	}
	
	public void numIntensity() {
		int total=0;
		for(int i=0;i<apertures.size();++i) {
			
			for(int j=0;j<apertures.get(i).size();++j) {
				if(apertures.get(i).get(j).intensity>0)total=total+1;
			}
			
		}
		System.out.println("numApertures: "+total );
		// TODO Auto-generated method stub
		
	}
	public void beamOnTimePrint() {
		int BOT=0;
		for(int i=0;i<apertures.size();++i) {
			
			for(int j=0;j<apertures.get(i).size();++j) {
				
				BOT=(int) (BOT+apertures.get(i).get(j).intensity);
			}
			
		}
		System.out.println("Beam on Time: "+BOT );
		// TODO Auto-generated method stub
		
	}
	public int beamOnTime() {
		int BOT=0;
		for(int i=0;i<apertures.size();++i) {
			
			for(int j=0;j<apertures.get(i).size();++j) {
				
				BOT=(int) (BOT+apertures.get(i).get(j).intensity);
			}
			
		}
		return BOT ;
		// TODO Auto-generated method stub
		
	}
	public void scorePrint(){
		for(int i=0;i<scores.length;i++) {
			System.out.print(scores[i]+" ");
		}
		//System.out.println();
	}
	//------------evualate function
	
	public double evaluateSolution( ){
		//long startTime = System.currentTimeMillis();
		switch(evaluationFunction) {
		case 0:
			quadraticPenalization();
			break;
		case 1:
			gEUD();
			break;
		case 3:
			//quadraticPenalizationOARs();
			quadraticPenalizationOARstPTV();
			break;
		default:
			
		}
		//long totalTime = System.currentTimeMillis() - startTime;
        
        //System.out.println("Tiempo total evaluate: " + totalTime);

		return singleObjectiveValue;
	}
	public double deltaEvaluateSolution( int  idchangedBeamlet, double[]oldIntensity){
		switch(evaluationFunction) {
		case 0:
			quadraticPenalization(idchangedBeamlet,oldIntensity );
			break;
		case 1:
			//gEUD();
			break;
		case 3:
			//quadraticPenalizationOARs();
			quadraticPenalizationOARstPTV(idchangedBeamlet,oldIntensity );
			break;
		default:
			
		}
		
		
		return singleObjectiveValue;
	}
	public double quadraticPenalization( ){
		
		// Se asignan las intensidades para cada apertura
	    aperturesToIntensityMatrixPerStation();
	    intensityMatrixStationToIntensityMatrixPlan();
	    
	    // Se obtienen los indices de los beamlets que afectan a los voxels de cada organo (DDM de cada organo)
	    // Cada valor del arraylist es un Hashtable que representa un organo
	    // El Hashtable tiene como clave el id del voxel y como valor un arraylist con los índices de los beamlets que afectan a ese voxel
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm; //ORGANO, VOXEL, []Beam que lo afectan
		// Se obtienen los valores de las tasas de los beamlets que afectan a los voxels de cada organo (valores de la DDM de cada organo)
		// Cada valor del arraylist es un Hashtable que representa un organo
		// El Hashtable tiene como clave el id del voxel y beamlet (id_voxel-id_beamlet) y como valor un numero double que representa la tasa
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm; //Voxel-beam = intensidad
		// Variable auxiliar que se va a utilizar para obtener el valor del arrayList "index_dao_ddm"
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		// Variable auxiliar que se va a utilizar para obtener el valor del arrayList "value_dao_ddm" 
		Hashtable<String, Double> aux_values;
		// Variable auxiliar utilizada para iterar los valores de "aux_index"
		Enumeration<Integer> keys;
		// Variable auxiliar que contiene los id's de los beamlets que afectan a un voxel
		ArrayList<Integer> beams;
		// Variable auxiliar utilizada para concatenar el id del voxel y el beamlet que lo afecta ("id_voxel-id_beamlet")
		String value_index_key;
		// Radiation: Variable auxiliar que contiene la tasa de radiación en la que un voxel es afectada por un beamlet
		// intensityVoxel2Beams: Variable auxiliar que almacena la dosis total que va a recibir un voxel
		Double radiation, intensityVoxel2Beams;
		// key: variable auxiliar que guarda el id del voxel
		// beam: variable auxiliar que almacena el id del beamlet
		Integer key, beam;
		// Variable deprecada
		Double[][] intensityVoxel2BeamsP;

		
		double[] solutionVector = intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		int solutionVectorSize=solutionVector.length;
		double F = 0.0, pen,score;
		intensityVoxel2BeamsP=new Double[o.length][];
		/// Se recorre cada organo
		for(int i=0; i<o.length; i++){
			// Se obtienen los beamlets (indices) que afectan a cada voxel
			aux_index = index_dao_ddm.get(i);
			// Se obtienen las intensidades que aporta cada beamlet a cada voxel
			aux_values = value_dao_ddm.get(i);
			pen = 0.0;
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			intensityVoxel2BeamsP[i]=new Double[solutionVector.length];
			//Recorremos claves de voxel por organo para su evaluaci�n
			for(int b=0;b<solutionVectorSize; b++){
            	intensityVoxel2BeamsP[i][b]=0.0;
            }
			while(keys.hasMoreElements()){
				// Se obtiene el id del voxel
	            key = keys.nextElement();
	            // Se obtiene un array que contiene los id's de los beamlet que afectan al voxel con clave "key"
	            beams = aux_index.get(key);
	            // Se inicializa el valor de "intensityVoxel2Beams" 
	            intensityVoxel2Beams = 0.0;
	            
	            
	            for(int b=0;b<beams.size(); b++){
	            	beam = beams.get(b);
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	intensityVoxel2Beams+= solutionVector[beam] * radiation; //segundos* ratio de radiacion
	            	//Representa [organo][beam]
	            	intensityVoxel2BeamsP[i][beam]+= solutionVector[beam] * radiation;
	            }
	            
            	if(i == 0){
            		pen += weights[i] * Math.pow((dd[i] - intensityVoxel2Beams),2);
            	}else{
   					pen += weights[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				}
            	
	        }
			//System.out.println(intensityVoxel2BeamsP);
			//System.out.println(i+" "+Zmax[i]+" ");
			
			score=pen/aux_index.size();
			scores[i]=score/weights[i];
			
			F+=pen/aux_index.size();
		}
		singleObjectiveValue=F;
		

		return F;

	}
	public double quadraticPenalization(int idchangedBeamlet, double[]oldIntensity){
		
	    aperturesToIntensityMatrixPerStation();
	    intensityMatrixStationToIntensityMatrixPlan();
	    // Se obtienen los indices de los beamlets que afectan a los voxels de cada organo (DDM de cada organo)
	    // Cada valor del arraylist es un Hashtable que representa un organo
	    // El Hashtable tiene como clave el id del voxel y como valor un arraylist con los índices de los beamlets que afectan a ese voxel
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm; //ORGANO, VOXEL, []Beam que lo afectan
		// Se obtienen los valores de las tasas de los beamlets que afectan a los voxels de cada organo (valores de la DDM de cada organo)
		// Cada valor del arraylist es un Hashtable que representa un organo
		// El Hashtable tiene como clave el id del voxel y beamlet (id_voxel-id_beamlet) y como valor un numero double que representa la tasa
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		// Variable auxiliar que se va a utilizar para obtener el valor del arrayList "value_dao_ddm" 
		Hashtable<String, Double> aux_values;
		// Variable auxiliar que contiene los id's de los beamlets que afectan a un voxel
		ArrayList<Integer> beams;
		// Variable auxiliar utilizada para concatenar el id del voxel y el beamlet que lo afecta ("id_voxel-id_beamlet")
		String value_index_key;
		// Variable auxiliar que contiene la tasa de radiación en la que un voxel es afectada por un beamlet
		Double radiation;
		// variable obsoleta
		//Double[][] intensityVoxel2BeamsP;
		//intensityVoxel2BeamsP=new Double[o.length][];
		/// Se recorre cada organo
		// No se utiliza en este metodo
		double[] scoresOld = this.scores.clone(); 
	
		// Se obtienen los indices de los voxels afectados por cada beamlet
	    // Cada valor del arraylist es un Hashtable que representa un organo
	    // El Hashtable tiene como clave el id del beamlet y como valor un arraylist con los índices de los voxels
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> beamlet_to_voxel_dao_ddm = M.beamlet_to_voxel_dao_ddm; 
		// Variable auxiliar que se va a utilizar para obtener el valor del arrayList "beamlet_to_voxel_dao_ddm"
		Hashtable<Integer, ArrayList<Integer>> beamletToVoxels;
		
		Hashtable<Integer, ArrayList<Integer>> voxelToBeamlets;
		ArrayList<Integer> voxels;
		Double  intensityVoxel2BeamsNew , intensityVoxel2BeamsOld , scoreNew;
		Integer aux_beam, voxel;
		
		double[] aux_solutionVector = intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)

		double FNew = 0.0, penNew,penOld, deltaPen;
			

		/// Se recorre cada organo
		for(int i=0; i<o.length; i++){
			
			// Se obtiene el Hashtable que toma como clave el id del beamlet y que almacena ArrayList's con los id's de los voxels
			beamletToVoxels = beamlet_to_voxel_dao_ddm.get(i);
			// Se obtiene Se obtiene el Hashtable que toma como clave el id voxel-beamlet y que almacena los valores de las tasas de radiacion
			aux_values = value_dao_ddm.get(i);
			// Se obtiene el Hashtable que toma como clave el id del voxel y que almacena ArrayList's con los id's de los beamlets
			voxelToBeamlets = index_dao_ddm.get(i);
			penNew = 0.0;
			penOld = 0.0;
			deltaPen= 0.0;
			
			// Se obtienen los voxels afectados por el beamlet con el id "idchangedBeamlet"
			voxels = beamletToVoxels.get(idchangedBeamlet);
				
			// Se iteran todos los voxels afectados por el beamlet con el id "idchangedBeamlet"
			for(int v=0;v<voxels.size(); v++){
				
				// Se inicializan las variables en 0 para ser sumadas posteriormente
		        intensityVoxel2BeamsNew = 0.0;
		       	intensityVoxel2BeamsOld = 0.0;
				// Se obtiene el id del voxel
				voxel = voxels.get(v);
		      	
				// Se obtienen todos los id's de los beamlets que afectan al voxel "voxel"
		       	beams = voxelToBeamlets.get(voxel);
		
		       	for(int b = 0; b <beams.size() ; b++) {
		       		// Se obtiene el id de uno de los beamlets que afectan al voxel
	        		aux_beam = beams.get(b);
	        		// Se asigna a la variable "value_index_key" la clave para obtener un valor del Hashtable "aux_values"
		           	value_index_key = voxel+"-"+aux_beam;
		           	// Se obtiene la tasa de radiacion que aporta el beamlet con id almacenado en la variable "aux_beam"
		           	radiation = aux_values.get(value_index_key);
		           	
		           	// Se suma el aporte del beamlet al voxel con id almacenado en la variable "voxel" para la nueva y vieja intensidad
		           	intensityVoxel2BeamsNew += aux_solutionVector[aux_beam] * radiation; // segundos * ratio de radiacion
		           	intensityVoxel2BeamsOld += oldIntensity[aux_beam] * radiation;	// segundos * ratio de radiacion
           			           			        			  
		        }
		       	
		        if(i == 0){
		        	// Se calcula el error cuadratico medio para el PTV. 
		           	penNew += weights[i] * Math.pow(((dd[i] - intensityVoxel2BeamsNew) ),2);
		           	penOld += weights[i] * Math.pow((dd[i] - intensityVoxel2BeamsOld),2);
		           	
		        }else{
		        	// Se calcula el error cuadratico medio para el OARS
		   	        penNew += weights[i] * Math.pow(Math.max(( (intensityVoxel2BeamsNew - dd[i])),0),2);
		   			penOld += weights[i] * Math.pow(Math.max((intensityVoxel2BeamsOld - dd[i]),0),2);
		   		}
		    }
			
			
			// Se calcula la diferencia entre la penalización anterior y la nueva
			deltaPen= penOld-penNew;
			
			scoreNew= deltaPen/voxelToBeamlets.size();

			FNew += deltaPen/voxelToBeamlets.size();
		}
			
			
		// Se calcula la nueva funcion objetivo
		FNew = singleObjectiveValue -FNew;
		//System.out.println("Solucion antigua: "+auxF);
		//System.out.println("Solucion nueva moyano F: "+F);
		//System.out.println("Solucion mia : "+FNew);
		//double diff = F - FNew;
		//System.out.println("Diferencia moyano - mia  : "+diff);
		//System.out.println("-----------------------------");
		//singleObjectiveValue=F;
		singleObjectiveValue=FNew;
		
		return FNew;
		//return F;

	}	
	/**
	 * This function only calculate the penalization of the OARs. This not take in consideration the PTV restriction
	 * @return
	 */
	
	public double quadraticPenalizationOARs( ){
		aperturesToIntensityMatrixPerStation();
	    intensityMatrixStationToIntensityMatrixPlan();
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		Double[][] intensityVoxel2BeamsP;
		
		
		double[] solutionVector = intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		int solutionVectorSize=solutionVector.length;
		double F = 0.0, pen,score;
		intensityVoxel2BeamsP=new Double[o.length][];
		for(int i=1; i<o.length; i++){
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			pen = 0.0;
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			intensityVoxel2BeamsP[i]=new Double[solutionVector.length];
			//Recorremos claves de voxel por organo para su evaluaci�n
			for(int b=0;b<solutionVectorSize; b++){
            	intensityVoxel2BeamsP[i][b]=0.0;
            }
			while(keys.hasMoreElements()){
	            key = keys.nextElement();
	            beams = aux_index.get(key);
	            intensityVoxel2Beams = 0.0;
	            
	            
	            for(int b=0;b<beams.size(); b++){
	            	beam = beams.get(b);
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	intensityVoxel2Beams+= solutionVector[beam] * radiation;
	            	
	            	intensityVoxel2BeamsP[i][beam]+= solutionVector[beam] * radiation;
	         
	            }
	            
	            // Here the penalization of the recommend dosis with the real dosis is calculate for each beamlet
   				pen += weights[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				
            	
	        }
			//System.out.println(intensityVoxel2BeamsP);
			//System.out.println(i+" "+Zmax[i]+" ");
			
			score=pen/aux_index.size();
			scores[i]=score/weights[i];
			F+=pen/aux_index.size();
		}
		scores[0]=0;
		singleObjectiveValue=F;
		//System.out.println("Solucion: "+F);
		return F;
	}
	/**
	 * This function  calculate the penalization of the OARs. 
	 * This take in consideration the PTV restriction, if the restriction is not complain return a worst value
	 * @return
	 */
	
	public double quadraticPenalizationOARstPTV( ){
		aperturesToIntensityMatrixPerStation();
	    intensityMatrixStationToIntensityMatrixPlan();
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		Double[][] intensityVoxel2BeamsP;
		int Factible=1;
		double worstValue=999999;
		
		double[] solutionVector = intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		int solutionVectorSize=solutionVector.length;
		double F = 0.0, pen,score;
		intensityVoxel2BeamsP=new Double[o.length][];
		for(int i=0; i<o.length; i++){
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			pen = 0.0;
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			intensityVoxel2BeamsP[i]=new Double[solutionVector.length];
			//Recorremos claves de voxel por organo para su evaluaci�n
			for(int b=0;b<solutionVectorSize; b++){
            	intensityVoxel2BeamsP[i][b]=0.0;
            }
			while(keys.hasMoreElements()){
	            key = keys.nextElement();
	            beams = aux_index.get(key);
	            intensityVoxel2Beams = 0.0;
	            
	            
	            for(int b=0;b<beams.size(); b++){
	            	beam = beams.get(b);
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	intensityVoxel2Beams+= solutionVector[beam] * radiation;
	            	
	            	intensityVoxel2BeamsP[i][beam]+= solutionVector[beam] * radiation;
	            }
	            if(i==0) {
	            	if(intensityVoxel2Beams<( dd[i]- 0.1)) {
	            		Factible=0;
	            		
	            		break;
	            	}
	            	
	            }else {
	            	// Here the penalization of the recommend dosis with the real dosis is calculate for each beamlet
	   				pen += weights[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
	   				
				}
	            
	            
            	
	        }
			//System.out.println(intensityVoxel2BeamsP);
			//System.out.println(i+" "+Zmax[i]+" ");
			
			score=pen/aux_index.size();
			scores[i]=score/weights[i];
			F+=pen/aux_index.size();
		}
		scores[0]=0;
		singleObjectiveValue=F;
		//System.out.println("Solucion: "+F);
		if(Factible==1) {
			//System.out.println("Solucion: "+F);
			return F;

		}
		else {
			singleObjectiveValue=worstValue;
			//System.out.println("Solucion: "+F);
			return worstValue;
		}
	}
	
	public double quadraticPenalizationOARstPTV(int  idchangedBeamlet,double[]oldIntensity){
		aperturesToIntensityMatrixPerStation();
	    intensityMatrixStationToIntensityMatrixPlan();
		
	    ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<String, Double> aux_values;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation;
		int Factible=1;
		double worstValue=999999;
	
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> beamlet_to_voxel_dao_ddm = M.beamlet_to_voxel_dao_ddm; 
		Hashtable<Integer, ArrayList<Integer>> beamletToVoxels;
		Hashtable<Integer, ArrayList<Integer>> voxelToBeamlets;
		ArrayList<Integer> voxels;
		Double  intensityVoxel2BeamsNew , intensityVoxel2BeamsOld , scoreNew;
		Integer aux_beam, voxel;
		
		double[] aux_solutionVector = intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)

		double FNew = 0.0, penNew,penOld, deltaPen;
			
		/// Se recorre cada organo
		for(int i=0; i<o.length; i++){
			beamletToVoxels = beamlet_to_voxel_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			voxelToBeamlets = index_dao_ddm.get(i);
			penNew = 0.0;
			penOld = 0.0;
			deltaPen= 0.0;
			
			voxels = beamletToVoxels.get(idchangedBeamlet);
				
			for(int v=0;v<voxels.size(); v++){
				voxel = voxels.get(v);
		        intensityVoxel2BeamsNew = 0.0;
		       	intensityVoxel2BeamsOld =0.0;
		      	       	
		       	beams = voxelToBeamlets.get(voxel);
	     		
		
		       	for(int b = 0; b <beams.size() ; b++) {
		 
	        		aux_beam = beams.get(b);
		           	value_index_key = voxel+"-"+aux_beam;
		           	
		           	radiation = aux_values.get(value_index_key);
		           	
		           	intensityVoxel2BeamsNew += aux_solutionVector[aux_beam] * radiation; //segundos* ratio de radiacion
		           	intensityVoxel2BeamsOld += oldIntensity[aux_beam] * radiation;	
           			           			        			  
		        }
			       
			        if(i==0) {
		            	if((intensityVoxel2BeamsNew<(dd[i] -0.1))) {
		            		Factible=0;
		            		break;
		            	}
	            	
	            }else {
	            	penNew += weights[i] * Math.pow(Math.max(( (intensityVoxel2BeamsNew - dd[i])),0),2);
			   		penOld += weights[i] * Math.pow(Math.max((intensityVoxel2BeamsOld - dd[i]),0),2);
	   				
				}
		            	
		    }
		            	            	           	
			deltaPen= penOld-penNew;
			
			//scoreNew= deltaPen/voxelToBeamlets.size();

			FNew += deltaPen/voxelToBeamlets.size();
		}

		if(Factible==1) {
			FNew = singleObjectiveValue - FNew;
			return FNew;

		}
		else {
			FNew = worstValue;
			return worstValue;
		}

	}
	public double quadraticPenalizationOAR( ){
		aperturesToIntensityMatrixPerStation();
	    intensityMatrixStationToIntensityMatrixPlan();
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		Double[][] intensityVoxel2BeamsP;
		int Factible=1;
		double worstValue=999999;
		
		double[] solutionVector = intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		int solutionVectorSize=solutionVector.length;
		double F = 0.0, pen,score;
		intensityVoxel2BeamsP=new Double[o.length][];
		for(int i=0; i<o.length; i++){
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			pen = 0.0;
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			intensityVoxel2BeamsP[i]=new Double[solutionVector.length];
			//Recorremos claves de voxel por organo para su evaluaci�n
			for(int b=0;b<solutionVectorSize; b++){
            	intensityVoxel2BeamsP[i][b]=0.0;
            }
			while(keys.hasMoreElements()){
	            key = keys.nextElement();
	            beams = aux_index.get(key);
	            intensityVoxel2Beams = 0.0;
	            
	            
	            for(int b=0;b<beams.size(); b++){
	            	beam = beams.get(b);
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	intensityVoxel2Beams+= solutionVector[beam] * radiation;
	            	
	            	intensityVoxel2BeamsP[i][beam]+= solutionVector[beam] * radiation;
	            		
	            		
	            	
	         
	            }
	            if(i==0) {

	            	pen += 1 * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
	            }else {
	            	// Here the penalization of the recommend dosis with the real dosis is calculate for each beamlet
	   				pen += weights[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
	   				
				}
	            
	            
            	
	        }
			//System.out.println(intensityVoxel2BeamsP);
			//System.out.println(i+" "+Zmax[i]+" ");
			
			score=pen/aux_index.size();
			scores[i]=score/weights[i];
			F+=pen/aux_index.size();
		}
		scores[0]=0;
		singleObjectiveValue=F;
		//System.out.println("Solucion: "+F);
		if(Factible==1) {
			//System.out.println("Solucion: "+F);
			return F;

		}
		else {
			singleObjectiveValue=worstValue;
			//System.out.println("Solucion: "+F);
			return worstValue;
		}
	}
	
	public double gEUD() {
		aperturesToIntensityMatrixPerStation();
	    intensityMatrixStationToIntensityMatrixPlan();
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		
	
		double[] solutionVector =intensity; 
		//System.out.println("Largo soluci�n :"+solutionVector.length);
		//Recorremos los organos
		for(int i=0; i<o.length; i++){
			ArrayList<Double> Doses=new ArrayList<Double>();
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			//Recorremos claves de voxel por organo para su evaluaci�n
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
	            Doses.add(intensityVoxel2Beams);
	            
        
	        }
			scores[i]=gEUDOrgan(Doses,i,aux_index.size());
			singleObjectiveValue+=scores[i];

		}
		return singleObjectiveValue;

		
	}
	public double evaluateSolutionPen( DDM M, Organs[] o, double w[],double[] solutionVector){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		double scores[]= {0,0,0};
		
		double F = 0.0, pen,score;
		
		//System.out.println("Largo soluci�n :"+solutionVector.length);
		//Recorremos los organos
		for(int i=0; i<o.length; i++){
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
			pen = 0.0;
			//System.out.println("Organo: "+o[i].name);
			keys = aux_index.keys();
			//Recorremos claves de voxel por organo para su evaluaci�n
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
            		//pen += w[i] * Math.pow((dd[i] - intensityVoxel2Beams),2);
            	}else{
   					pen += w[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				}
            	
	        }
			//System.out.println(i+" "+Zmax[i]+" ");
			score=pen/aux_index.size();
			
			scores[i]=score;
			F+=pen/aux_index.size();
			System.out.println("organo: "+scores[i]);
		}
		//solution.singleObjectiveValue=F;
		System.out.println("Solucion: "+F);
		return F;
	}
	
	
	
	
	
	public double gEUDOrgan(ArrayList<Double> doses,int organ,int m) {
		double gEUDr=0;
		
		for(int i=0;i<doses.size();i++) {
			gEUDr+=Math.pow(doses.get(i),a[organ]);
		}
		gEUDr=gEUDr/m;
		double exp=(double)(1/a[organ]);
		return Math.pow(gEUDr,exp);
		
	}
	
	public void solverIntensity() {
	 	aperturesToIntensityMatrixPerStationwithoutIntensty();
	    ApertureMatrixStationToIntensityMatrixPlan();
	

	    //time=System.currentTimeMillis();
	    
	    Gurobi_Solver newModel;
		try {
			newModel = new Gurobi_Solver(this,M,o,beamletAngles,dd,weights,0);
			setIntensity(newModel.newIntensity);
		    setOF(newModel.objVal);
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
}
	public void solverIntensity(int min) {
	 	aperturesToIntensityMatrixPerStationwithoutIntensty();
	    ApertureMatrixStationToIntensityMatrixPlan();
	

	    //time=System.currentTimeMillis();
	    
	    Gurobi_Solver newModel;
		try {
			
			newModel = new Gurobi_Solver(this,M,o,beamletAngles,dd,weights,min);

			setIntensity(newModel.newIntensity);
		    setOF(newModel.objVal);
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
}
	private void solverIntensityPen() {
		aperturesToIntensityMatrixPerStationwithoutIntensty();
	    ApertureMatrixStationToIntensityMatrixPlan();
	

	    //time=System.currentTimeMillis();
	    
	    Gurobi_Solver newModel;
		try {
			
			newModel = new Gurobi_Solver(this,M,o,beamletAngles,dd,weights,1,1);

			setIntensity(newModel.newIntensity);
		    setOF(newModel.objVal);
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	public void solverIntensityPen(int min) {
		aperturesToIntensityMatrixPerStationwithoutIntensty();
	    ApertureMatrixStationToIntensityMatrixPlan();
	

	    //time=System.currentTimeMillis();
	    
	    Gurobi_Solver newModel;
		try {
			
			newModel = new Gurobi_Solver(this,M,o,beamletAngles,dd,weights,1,min);

			setIntensity(newModel.newIntensity);
		    setOF(newModel.objVal);
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	//----threds
	
	public void run() {
		
		//System.out.println("method: "+evaluationMethod);
		
		switch(evaluationMethod) {
			case 1:
				evaluateSolution();
				break;
			case 2:
				solverIntensity();
				break;
			case 3:
				solverIntensityPen();
				break;
				
			default:
				evaluateSolution();
		}
		
	}





	/**
	 * Genera un vector con la posicion de las aperturas que tienen una intensidad bajo un valor delta
	 * @param delta
	 */
	public void getListLowIntensity(double delta) {
		lowIntensityApertureNotChanged=new ArrayList<int[]>();
    	for(int i=0; i<apertures.size(); i++){
			
				List<Aperture> list_apertures = apertures.get(i);
				//recorrido de aperturas por angulo
				
				for(int j=0; j<list_apertures.size(); j++){
					if(list_apertures.get(j).intensity<=delta){
						int[] aux=new int[2];
						aux[0]=i;
						aux[1]=j;
						lowIntensityApertureNotChanged.add(aux);

					}
					
				}

			}

	}
	
	
	/**
	 * Devuelve cuantas aperturas con baja intensidad no han sido modificadas
	 * @return
	 */
	public int getLowIntensityApertureNotChanged() {
	
		return lowIntensityApertureNotChanged.size();
	}
	
	
	/**
	 * Devuelve la posicion de una apertura escogida al azar, excepto cuando se tienen aperturas con intensidades bajas. 
	 * En este caso devuelve al azar una apertura de las que se encuentran con intensidad baja.
	 * @return
	 */
	public int[] getApertureRandom() {
		Random r = new Random(System.currentTimeMillis());
		int position[]=new int[2];



		if(getLowIntensityApertureNotChanged()==0) {
			int angleChoose= r.nextInt(apertures.size());
			int apertureChoose= r.nextInt(apertures.get(angleChoose).size());
			position[0]=angleChoose;
			position[1]=apertureChoose;
		}
		else {

			position=lowIntensityApertureNotChanged.get(r.nextInt(lowIntensityApertureNotChanged.size()));

		}
		return position;
	}
	
	public void deleteLowIntensityApertureNotChanged(int[] position) {
		int aux=0;
		if(!lowIntensityApertureNotChanged.isEmpty()) {
			for(int i=0;i<lowIntensityApertureNotChanged.size();i++)
				if(lowIntensityApertureNotChanged.get(i)[0]==position[0] && lowIntensityApertureNotChanged.get(i)[1]==position[1]) {
					aux=i;
					break;
				}
			lowIntensityApertureNotChanged.remove(aux);
		}
		
	}
	
	
	
	
	
	
}
