package MOFMO;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import IMRT_DAO.DDM;
import IMRT_DAO.Organs;
import IMRT_DAO.TreatmentPlan;
import SingleObj.Algorithm;
import SingleObj.InitialSolution;
import gurobi.GRBException;

public class FMO {
	public int numOrgans;
    public int numAngles;                   //#Beams in TreatmentPlan
    public Organs[] o;
    public String pathFile = "";
    public int option;
    public int selectionCriterion;
    public int[][] Vx; 				       //Vx index. % of the volume receiving more than x% of the prescribed dose
    						  				       //One row per organ (usually only for target regions)
    public int[][] Dx; 				       //Dx index. Minimum dose of the hottest x% of the volume.
    						  				       //One row per organ (usually only for OAR regions)
    public boolean randomInitialAngles;
    public int[] initialBACs;
    public int[][] beamletSet;
    
    public int num_apertures = 5;            // # DAO : Number of apertures per station
	public int init_intensity;               // # DAO : Initial intensity for apertures
    public int max_intensity;                // # DAO : Maximum intensity for apertures
    public int max_delta;                    // # DAO : Delta máximo para la variación de intesidad por apertura
    public int max_iter;              		// # DAO : Cantidad máxima de iteraciones del algoritmo
    public int max_time;                     // # DAO : tiempo máximo de ejecución del algoritmo
    public int seed;                         // # DAO : semilla para partir la busqueda local
    public int step_intensity;				// # DAO : Step size for aperture intensity (2)
    public Vector<int[][]> initial_aperture_shape;  // # DAO : Conjunto de formas de las aperturas según restricciones (rango sup y rango inf)
    public Vector<int[][]> apertureShapes;   // # DAO: Aperturas segun restricciones del problema (en forma de matriz)
    public int [] zmax,zmin;
    public double[] dd;
    public int totalbmlt;
    public int[] selAngles;
    public double[] w;
    public int[] beamletAngles;
    public DDM M;

	public FMO(String inputFile) throws IOException {

	    readInputFile(inputFile);
	    beamletSet = new int[numAngles][4];
		beamletAngles = new int[numAngles];
	    String beamInfoDir = pathFile+"beamsInfo.txt";
	    //beamletSet = loadNumBixels(beamInfoDir);
	    beamletSet = loadNumBixelsTRT(beamInfoDir);
	    

	    M = new DDM(o, initialBACs, pathFile,beamletAngles,totalbmlt);
	    
	    //Obtener dosis de cada organo
	    setDD();
	    
	    w = new double[]{0,0,1};
	    //-----Generar solucion Inicial------/
	    
		
	    //evaluateSolution(solution, M, o, w);
	    //System.out.println("1");
	    //ArrayList<ArrayList<Double>> Doses=DoseVolume(solution, M, o);
	    //doseprint(Doses,"FMOTRT2_");
	    //System.out.println("2");
	    
	    
	    evaluateSolution(M, o, w);
	    
	    
	    
	    
	    //ArrayList<double []> initialFluenceMaps= readInitialFluenceMaps(); 
	    
	    //surfaceFMO surface=new surfaceFMO(initialFluenceMaps);
		
	    //saveDosesVolumes(surface);
		
		
		
		;
		
	    
	    
	    
	    
	    
	    
	}
	
	/**
	 * Genera y guarda las dosis de los mapa de influencia de la surface dada
	 * @param surface
	 */
	private void saveDosesVolumes(surfaceFMO surface) {
		
		for (int i=0;i<surface.allSolution.size();i++) {
			String init="CERR-FMO-"+String.valueOf(i)+"-";
			doseprint(DoseVolume(surface.allSolution.get(i), M, o),init);
		}
		
		
	}
	private ArrayList<double[]> readInitialFluenceMaps() {
		// TODO Auto-generated method stub
		return null;
	}
	public void readInputFile(String dir) throws IOException{
        
    	String sp="\\s+";
        //String dir = "./inputFile.txt";
        File f = new File(dir);
        System.out.println(f.getAbsolutePath());
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        String line = "";
        line=fileIn.readLine();
        //First read numbero of Organs and angles
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                numOrgans = Integer.parseInt(auxReader[0]);
                numAngles = Integer.parseInt(auxReader[1]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        //Go to the next input line
        while(line != null){
            if (!line.contains("%")){
                break;
            }
            line=fileIn.readLine();
        }
 
        //Info Organs
        o = new Organs[numOrgans];
        while(line != null){
            if (!line.contains("%")){
                for (int y=0;y<numOrgans;y++){
                    String[] auxReader = line.split(sp);
                    o[y]=new Organs(
                        auxReader[0], 
                        Integer.parseInt(auxReader[1]), 
                        Integer.parseInt(auxReader[2]),
                        Integer.parseInt(auxReader[3]),
                        Integer.parseInt(auxReader[4]),
                        Integer.parseInt(auxReader[5]),
                        Integer.parseInt(auxReader[6]),
                        Integer.parseInt(auxReader[7]),
                        Integer.parseInt(auxReader[8]),
                        Integer.parseInt(auxReader[9]), 
                        Integer.parseInt(auxReader[10]),
                        Integer.parseInt(auxReader[11]), 
                        Integer.parseInt(auxReader[12]), 
                        Boolean.parseBoolean(auxReader[13]));
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        //get filepath
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                pathFile=auxReader[0];
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get option
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                option=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get selectionCriterion (1 = LS, 2 = nextDescent, 3= gradientBas)
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                selectionCriterion=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get max iterations LS
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                max_iter=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get max time LS
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                max_time=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get Vx Index. % of the volume receiving more than x% of the prescribed dose
        Vx = new int[numOrgans][];
        while(line != null){
            if (!line.contains("%")){
                while(!line.contains("%")){
                    String[] auxReader = line.split(sp);
                    Vx[Integer.parseInt(auxReader[0])]=new int[auxReader.length-1];
                    for(int i=1;i<auxReader.length;i++){
                        Vx[Integer.parseInt(auxReader[0])][i-1]=Integer.parseInt(auxReader[i]);
                    }
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        //get Dx Index. Minimum dose of the hottest x% of the volume
        Dx = new int[numOrgans][];
        while(line != null){
            if (!line.contains("%")){
                while(!line.contains("%")){
                    String[] auxReader = line.split(sp);
                    Dx[Integer.parseInt(auxReader[0])]=new int[auxReader.length-1];
                    for(int i=1;i<auxReader.length;i++){
                        Dx[Integer.parseInt(auxReader[0])][i-1]=Integer.parseInt(auxReader[i]);
                    }
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        
        //get initial BACs
        initialBACs = new int[numAngles];
        while(line != null){
            if (!line.contains("%")){
                while(!line.contains("%")){
                	String[] auxReader = line.split(sp);
                    for (int i=0; i<numAngles;i++){
                        initialBACs[i]= Integer.parseInt(auxReader[i]);
                    }
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        
        //Initial value aperture intensity
        while(line != null){
            if (!line.contains("%")){
            	String[] auxReader = line.split(sp);
            	init_intensity=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        //Maximum intensity for apertures
        while(line != null){
            if (!line.contains("%")){
            	String[] auxReader = line.split(sp);
            	max_intensity=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        //Máximo delta para la variación de intensidad en las aperturas
        while(line != null){
            if (!line.contains("%")){
            	String[] auxReader = line.split(sp);
            	max_delta=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
      //Step size for aperture intensity
        while(line != null){
            if (!line.contains("%")){
            	String[] auxReader = line.split(sp);
            	step_intensity=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        fileIn.close();
	}
	 /**
	  * Modificacion para los casos nuevos
	  * @param beamsInfoDir
	  * @return
	  * @throws IOException
	  */
	 public static int[][] loadNumBixelsTRT(String beamsInfoDir) throws IOException{
	        int[][] beamlets = new int[360][4];
	        //int[][] beamlets = new int[numBeams][4];
	        int auxBeamletIndex=0;
	        
	        File beamInfo = new File(beamsInfoDir);
	        if (! beamInfo.exists()) {
	            System.err.println("Couldn't find 'beamsInfo.txt' file");
	        }else{
	            String line ="";
	            String[] auxReader=null;   
	            File f= new File(beamsInfoDir);
	            BufferedReader fileIn= new BufferedReader(new FileReader(f));
	            for (int i=0;i<356;i++){
	                line=fileIn.readLine();
	                auxReader = line.split("\t");
	                int angle=(int) Double.parseDouble(auxReader[0]);
	                beamlets[angle][0]=(int) Double.parseDouble(auxReader[0]); //beamIndex
	                beamlets[angle][1]=(int) Double.parseDouble(auxReader[1]); //numBeamlets
	                beamlets[angle][2]= auxBeamletIndex;                       //firstBeamletIndex
	                beamlets[angle][3]=(int) Double.parseDouble(auxReader[2]) - 1; //lastBeamletIndex
	                i=angle;
	                auxBeamletIndex = auxBeamletIndex + (int) Double.parseDouble(auxReader[1]);
	            }
	            fileIn.close();
	        }
	        return (beamlets);
	 }
	 public void setDD() {
			dd = new double[o.length];
			for(int i=0;i<o.length;i++) {
				if(o[i].isTarget) {
					dd[i]=o[i].doseLB;
				}else {
					dd[i]=o[i].doseUB;
				}
			}
	}
	 
	public void listDosePrint(ArrayList<ArrayList<ArrayList<Double>>> listDoses,String init) {
		for (int i = 0; i < listDoses.size(); i++) {
			String namelist=init+ String.valueOf(i);
			doseprint(listDoses.get(i),namelist);
		}
	}
	 
	public void doseprint(ArrayList<ArrayList<Double>> Doses,String init){
		System.out.println("Working Directory = " + System.getProperty("user.dir"));
		PrintWriter writer;
		try {
			writer = new PrintWriter((init+"0.txt"), "UTF-8");
			ArrayList<Double> arr=Doses.get(0);
			for(Double aux:arr) {
				writer.println(Double.toString(aux));
			}
			writer.close();
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			writer = new PrintWriter((init+"1.txt"), "UTF-8");
			ArrayList<Double> arr=Doses.get(1);
			for(Double aux:arr) {
				writer.println(Double.toString(aux));
			}
			writer.close();
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			writer = new PrintWriter((init+"2.txt"), "UTF-8");
			ArrayList<Double> arr=Doses.get(2);
			for(Double aux:arr) {
				writer.println(Double.toString(aux));
			}
			writer.close();
		} catch (FileNotFoundException| UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
		
		
		
	public ArrayList<ArrayList<Double>> DoseVolume(double[] fluenceMap, DDM M, Organs[] o){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		double pen=0;
		ArrayList<ArrayList<Double>> organ=new ArrayList<ArrayList<Double>>();
		double[] solutionVector = fluenceMap; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)

		//System.out.println("Largo solución :"+solutionVector.length);
		//Recorremos los organos
		for(int i=0; i<o.length; i++){
			ArrayList<Double> Doses=new ArrayList<Double>();
			aux_index = index_dao_ddm.get(i);
			aux_values = value_dao_ddm.get(i);
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
	            
            	if(i == 0){
            	}else{
   					pen = w[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				}
	            
	            Doses.add(intensityVoxel2Beams);

        
	        }
			organ.add(Doses);

		}
		System.out.println("OF: "+String.valueOf(pen));
		return organ;
	}
	
	
	//---------------- Funcion de evaluacion--------------------------------------
	public double evaluateSolution( DDM M, Organs[] o, double w[]){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		ArrayList<Double> intensityRealDose=new ArrayList<Double>();
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		double[] solutionVector = new double[]{0.0032716,0.00105346,19.9989,17.0936,19.9919,5.48065,0.000761112,1.22194,19.9993,19.9341,9.03373,19.987,19.9988,19.9991,5.47391,3.17271,17.3202,19.9988,19.9581,18.2314,17.6408,19.9993,18.4743,6.29248,3.77701,19.4714,19.9962,19.9651,19.5352,19.9978,19.9989,19.24,6.59939,4.9474,19.9998,19.9925,19.968,19.9126,19.9963,19.798,19.9991,6.23821,4.53865,19.9996,19.9938,19.9503,19.9946,18.3809,19.4162,19.0411,5.74537,4.38015,18.7806,19.9966,8.65349,19.9954,19.99,18.2084,5.50703,5.02063,16.2003,19.9989,19.9753,19.9397,4.75063,4.60621,0.000434602,19.9994,0.000885641,19.9991,0.0325965,15.8909,14.7907,9.55396,9.80396,13.2667,19.9765,8.7949,18.1114,16.0388,10.409,10.2629,12.2066,4.72244,19.9712,9.49999,14.5707,15.845,14.415,10.6796,11.6538,5.40328,19.9688,13.0347,19.6745,15.2463,15.4161,10.921,11.1171,8.51937,0.00189968,19.9565,5.89103,17.374,15.209,18.3493,11.2274,10.7149,11.7285,0.0010501,19.952,19.9652,19.9976,19.6223,18.4756,11.6491,11.3298,12.3894,3.04464,19.9733,19.993,19.9869,19.8816,19.9567,19.998,10.4537,0.0182648,19.9867,19.965,0.0300567,0.0555712,19.957,18.3137,19.926,19.9809,9.07695,0.00600348,4.19471,19.9994,17.5258,19.9759,11.3797,0.0537545,3.9166,4.08603,18.0912,19.9978,15.0718,19.9478,19.8876,1.66195,5.39217,6.22552,16.6633,18.1925,9.36954,0.0049993,19.9669,11.9671,0.011836,6.1872,7.88034,15.6753,16.8341,8.5559,0.00260411,19.9481,19.9722,0.00411085,7.05969,8.05739,13.9679,17.3154,3.77335,0.0055701,19.9101,19.9481,14.7912,7.70609,14.0893,13.3281,15.8896,4.70618,0.0257914,19.9205,19.9697,0.414837,16.5207,14.4942,13.6891,16.8552,0.00181885,19.9959,3.45076,19.9976,4.94135,0.557815,0.328738,19.9621,0.0155394,0.0196692,0.0489466,0.615781,19.9365,0.126253,0.026799,0.0131954,19.9048,3.79554,6.79762,14.6965,5.53623,6.30233,19.998,18.3793,19.9934,8.32018,6.50102,12.2816,14.9597,0.000484785,19.9977,19.9429,19.9438,0.0117425,5.86801,3.52048,12.7305,14.4875,1.36508,7.8361,19.9564,18.6078,1.91659,6.97553,3.74408,14.2687,15.5876,0.74356,10.6546,19.9503,1.63619,0.0385567,5.23795,4.14135,14.2882,15.4392,2.30323,14.4496,19.9619,0.115126,4.08743,2.99868,5.00659,16.2821,16.3914,4.90609,2.15757,19.9597,19.9705,1.72575,2.38917,6.56455,16.6498,19.9978,5.32083,19.9975,0.00945323,3.90588,4.8851,19.998,19.9305,9.49458,19.9822,0.0265098,0.0125084,0.0136265,0.0183483,0.0199789,19.967,19.9894,4.43056,0.453624,4.68989,11.0411,19.9554,0.0379024,19.8926,19.9795,19.9522,5.65526,8.20711,8.03516,13.7366,6.91752,19.9929,19.9385,19.9728,19.9787,7.28855,7.32052,10.3224,11.1599,17.9813,0.000633798,19.9673,19.7243,14.029,7.93888,7.64038,11.6848,11.2073,15.0399,0.00154801,19.9697,8.04304,12.7305,5.77088,8.57665,12.459,11.5401,12.7952,19.9798,4.0109,13.3349,6.22237,8.36473,12.8103,10.6064,9.04039,19.9911,0.0246703,5.0079,10.5231,13.2067,9.91458};
		//double[] solutionVector = new double[]{10,20,20,14,20,19,11,10,20,5,6,7,16,19,18,14,9,20,10,4,5,14,17,16,11,10,20,12,4,5,13,17,15,9,10,20,13,5,4,13,18,12,8,11,19,14,7,4,14,17,9,10,12,20,9,4,17,18,8,14,10,20,20,20,5,14,13,20,20,19,11,18,20,5,10,7,10,15,20,15,3,15,6,3,12,17,20,9,2,16,4,3,9,20,20,6,3,19,2,4,1,5,17,17,0,0,17,0,4,0,0,17,20,1,1,19,1,0,1,0,13,0,1,19,0,4,0,7,11,4,0,20,15,3,0,0,0,9,10,6,20,0,1,3,8,10,19,7,20,6,3,6,8,9,20,5,9,20,8,1,4,5,6,17,6,7,13,8,1,2,3,2,10,5,8,16,10,2,3,3,1,5,4,4,20,9,7,9,3,4,3,5,5,3,2,4,7,4,6,2,7,4,0,0,20,9,20,20,20,20,12,1,20,20,20,20,19,18,9,10,17,20,18,15,14,18,17,7,13,16,19,14,15,14,20,20,6,17,15,20,14,16,12,16,19,5,18,19,20,20,19,13,14,19,3,19,20,20,20,20,12,7,20,4,17,19,20,6,0,17,1,7,18,2,0,0,3,4,0,0,0,0,0,0,0,2,0,0,11,17,19,10,11,15,0,10,10,8,17,18,10,13,14,11,2,8,14,20,20,13,13,10,12,0,4,14,19,14,11,8,8,12,0,17,19,9,8,7,7,11,20,19,4,5,6,6};
		//double[] solutionVector=new double[] {10.31,20.00,20.00,14.49,20.00,19.37,11.44,10.08,20.00,5.23,6.49,7.40,16.28,18.54,17.56,13.60,9.32,20.00,10.41,3.67,4.82,13.81,16.93,15.60,10.51,9.81,20.00,11.63,4.02,4.74,12.96,16.99,14.88,8.83,9.79,20.00,13.07,4.90,3.71,13.46,17.62,12.28,8.30,11.12,19.34,13.87,7.48,4.13,13.53,16.78,9.34,9.68,12.06,20.00,8.51,4.18,17.12,18.09,7.79,14.35,10.17,20.00,20.00,20.00,4.55,14.37,13.27,20.00,20.00,19.11,10.89,18.46,20.00,4.95,10.04,6.66,9.82,15.16,19.87,14.52,3.23,14.87,6.21,3.34,11.99,16.70,20.00,8.77,2.48,16.18,3.84,2.63,9.30,20.00,20.00,5.58,2.86,19.32,2.23,3.68,0.68,5.38,17.18,17.28,0.00,0.00,16.50,0.00,4.28,0.00,0.00,17.06,20.00,1.04,0.58,19.19,0.92,0.00,0.54,0.45,13.30,0.00,1.39,18.92,0.00,4.37,0.00,7.01,11.15,4.37,0.00,20.00,15.45,2.62,0.00,0.00,0.00,8.54,9.78,5.60,20.00,0.00,1.43,3.35,8.33,10.28,19.01,6.61,20.00,5.81,2.98,6.48,8.28,8.88,20.00,5.25,8.81,20.00,8.45,1.12,3.79,5.25,5.60,16.52,6.01,7.16,12.88,7.73,0.68,1.61,3.07,1.61,10.29,5.49,7.75,16.13,9.73,1.64,3.13,2.94,1.12,5.49,3.85,3.76,20.00,8.51,6.98,9.28,2.81,3.76,3.33,5.10,5.30,2.95,2.44,3.95,7.35,3.54,5.58,1.85,6.66,3.53,0.00,0.00,20.00,9.14,19.85,20.00,20.00,20.00,12.22,1.47,20.00,20.00,20.00,20.00,18.61,17.75,8.77,9.74,16.53,19.93,18.08,15.31,13.80,18.36,16.93,6.52,13.03,16.03,19.12,13.84,14.91,13.61,20.00,19.84,6.39,16.74,15.23,19.66,13.69,15.69,12.40,15.65,18.58,4.84,18.19,18.92,20.00,20.00,18.85,12.91,13.85,19.30,3.44,19.49,20.00,20.00,20.00,20.00,12.08,7.44,19.75,3.63,17.48,18.65,20.00,6.44,0.00,17.11,1.27,7.16,18.34,1.94,0.00,0.00,2.97,3.64,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.83,0.00,0.00,10.85,17.47,18.61,9.50,11.06,14.53,0.00,9.56,9.64,8.43,17.10,17.73,10.13,12.52,13.76,10.97,1.70,7.53,13.55,20.00,20.00,13.21,12.79,9.77,12.10,0.00,3.59,14.10,19.19,13.58,10.81,8.01,8.19,11.79,0.00,17.35,19.34,9.14,8.15,6.98,7.11,11.09,20.00,18.97,3.56,5.08,5.88,5.87};
		//double[] solutionVector=new double[] {20.0,20.0,0.0,12.0,8.0,12.0,12.0,4.0,0.0,12.0,8.0,4.0,8.0,4.0,8.0,8.0,0.0,20.0,16.0,4.0,0.0,4.0,0.0,4.0,4.0,0.0,20.0,16.0,4.0,0.0,4.0,0.0,4.0,4.0,0.0,20.0,20.0,4.0,0.0,8.0,4.0,4.0,4.0,0.0,20.0,20.0,8.0,0.0,8.0,4.0,4.0,4.0,4.0,20.0,12.0,4.0,12.0,8.0,8.0,12.0,4.0,20.0,16.0,16.0,16.0,20.0,20.0,0.0,20.0,20.0,8.0,20.0,16.0,20.0,20.0,16.0,0.0,16.0,20.0,12.0,16.0,20.0,16.0,12.0,4.0,16.0,20.0,8.0,16.0,16.0,12.0,8.0,4.0,20.0,20.0,8.0,12.0,16.0,8.0,8.0,20.0,4.0,12.0,16.0,4.0,8.0,8.0,0.0,4.0,0.0,0.0,12.0,20.0,4.0,4.0,4.0,0.0,0.0,0.0,0.0,16.0,4.0,4.0,4.0,0.0,4.0,0.0,0.0,4.0,4.0,0.0,4.0,4.0,8.0,0.0,0.0,0.0,4.0,0.0,0.0,20.0,0.0,4.0,8.0,12.0,8.0,16.0,20.0,20.0,12.0,8.0,16.0,20.0,16.0,20.0,20.0,20.0,20.0,12.0,4.0,12.0,20.0,12.0,16.0,20.0,20.0,16.0,12.0,8.0,12.0,20.0,12.0,16.0,12.0,12.0,16.0,12.0,8.0,16.0,20.0,12.0,12.0,12.0,8.0,20.0,12.0,12.0,20.0,20.0,12.0,12.0,16.0,8.0,20.0,12.0,12.0,16.0,8.0,16.0,8.0,16.0,12.0,0.0,0.0,20.0,16.0,20.0,20.0,20.0,20.0,16.0,8.0,20.0,20.0,20.0,20.0,16.0,16.0,16.0,16.0,20.0,20.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,20.0,12.0,16.0,16.0,20.0,20.0,16.0,20.0,20.0,20.0,12.0,16.0,12.0,16.0,20.0,16.0,20.0,20.0,20.0,20.0,20.0,16.0,20.0,20.0,16.0,20.0,20.0,20.0,20.0,20.0,16.0,12.0,20.0,16.0,16.0,16.0,20.0,0.0,0.0,12.0,4.0,12.0,20.0,0.0,0.0,0.0,4.0,4.0,0.0,4.0,4.0,0.0,0.0,0.0,0.0,8.0,0.0,0.0,8.0,12.0,8.0,4.0,4.0,4.0,4.0,8.0,4.0,4.0,12.0,8.0,4.0,4.0,8.0,20.0,0.0,4.0,12.0,16.0,16.0,12.0,8.0,8.0,20.0,20.0,0.0,12.0,20.0,16.0,8.0,8.0,4.0,16.0,0.0,16.0,20.0,16.0,8.0,4.0,4.0,20.0,20.0,20.0,16.0,8.0,8.0,8.0,};
		//double[] solutionVector=new double[] {8.66659E-05,2.81824E-05,2.00000E+01,1.81798E+01,1.99998E+01,5.32946E+00,2.19788E-05,1.00943E+00,2.00000E+01,1.99981E+01,1.22066E+01,1.99995E+01,2.00000E+01,2.00000E+01,5.13651E+00,3.18559E+00,1.72321E+01,2.00000E+01,1.99986E+01,1.75493E+01,1.84140E+01,2.00000E+01,1.83622E+01,5.93818E+00,3.78481E+00,1.94713E+01,1.99999E+01,1.99987E+01,1.99871E+01,1.99999E+01,2.00000E+01,1.92736E+01,6.22730E+00,4.98393E+00,2.00000E+01,1.99998E+01,1.99989E+01,1.97366E+01,1.99999E+01,1.98668E+01,2.00000E+01,5.90717E+00,4.58216E+00,2.00000E+01,1.99998E+01,1.99983E+01,1.99999E+01,1.85675E+01,1.95452E+01,1.89558E+01,5.46571E+00,4.44144E+00,1.87579E+01,1.99999E+01,1.27726E+01,1.99999E+01,1.99998E+01,1.79764E+01,5.29205E+00,5.11615E+00,1.61360E+01,2.00000E+01,1.99991E+01,1.97020E+01,4.65785E+00,4.60800E+00,1.14362E-05,2.00000E+01,2.26258E-05,2.00000E+01,1.99768E-03,1.64462E+01,1.43811E+01,8.90131E+00,9.95950E+00,1.32564E+01,1.99994E+01,6.69730E+00,1.72113E+01,1.59688E+01,1.00558E+01,1.03678E+01,1.22210E+01,4.90413E+00,1.99990E+01,8.77176E+00,1.50511E+01,1.59691E+01,1.42812E+01,1.07533E+01,1.17155E+01,5.68341E+00,1.99991E+01,7.70529E+00,1.90610E+01,1.57439E+01,1.54177E+01,1.10302E+01,1.11886E+01,8.83708E+00,5.12429E-05,1.99987E+01,2.25569E+00,1.67965E+01,1.58609E+01,1.87162E+01,1.13855E+01,1.08241E+01,1.20707E+01,2.75813E-05,1.99987E+01,1.99986E+01,1.99999E+01,1.93513E+01,1.81108E+01,1.19226E+01,1.13351E+01,1.28603E+01,2.81094E+00,1.99993E+01,1.99998E+01,1.99997E+01,1.99981E+01,1.94971E+01,1.99999E+01,1.11605E+01,5.76400E-04,1.99997E+01,1.99990E+01,6.99829E-04,1.44835E-03,1.99989E+01,9.81492E+00,1.99978E+01,1.99995E+01,1.16235E+01,1.58590E-04,5.03146E+00,2.00000E+01,1.74192E+01,1.99993E+01,1.52888E+01,1.66890E+00,3.28249E+00,4.62977E+00,1.81082E+01,1.99999E+01,1.49125E+01,1.99980E+01,1.99945E+01,1.62506E+00,4.70711E+00,6.39842E+00,1.66878E+01,1.82210E+01,9.31557E+00,1.34773E-04,1.99991E+01,1.17832E+01,1.20441E+00,5.14738E+00,8.22850E+00,1.55880E+01,1.69758E+01,8.42361E+00,6.71592E-05,1.99980E+01,1.99969E+01,8.94215E-05,5.89475E+00,8.35565E+00,1.37826E+01,1.75213E+01,3.65566E+00,1.46632E-04,1.99969E+01,1.99947E+01,1.43813E+01,6.05907E+00,1.36115E+01,1.31915E+01,1.59977E+01,4.67702E+00,4.58297E-04,1.99980E+01,1.99995E+01,2.64250E+00,1.27550E+01,1.28478E+01,1.34051E+01,1.69693E+01,4.67386E-05,1.99999E+01,1.59450E+00,1.99999E+01,4.51745E+00,1.46109E-01,5.22657E-03,1.99990E+01,4.28931E-04,5.24533E-04,1.32399E-03,5.06300E-03,1.99985E+01,2.24258E-03,7.85072E-04,3.60910E-04,1.99919E+01,4.35750E+00,7.68083E+00,1.44603E+01,5.81057E+00,5.77213E+00,1.99999E+01,8.45050E+00,1.99998E+01,7.13393E+00,6.16955E+00,1.19808E+01,1.48798E+01,1.28485E-05,1.99999E+01,1.99980E+01,1.99984E+01,3.18329E-03,5.34809E+00,3.91408E+00,1.25506E+01,1.44333E+01,1.33090E+00,8.30772E+00,1.99986E+01,1.11191E+01,2.11995E+00,6.69814E+00,4.22935E+00,1.42699E+01,1.54557E+01,8.41244E-01,1.04787E+01,1.99983E+01,2.08129E+00,1.47129E-02,5.18101E+00,4.69238E+00,1.43508E+01,1.53723E+01,2.33905E+00,1.45893E+01,1.99991E+01,2.20942E-02,3.12271E+00,3.44540E+00,5.88804E+00,1.63957E+01,1.63540E+01,5.04138E+00,1.85650E+00,1.99990E+01,1.99990E+01,8.48160E-01,3.27772E+00,7.84280E+00,1.67746E+01,1.99999E+01,5.66483E+00,1.99999E+01,3.53517E-04,4.47333E+00,6.36304E+00,1.99999E+01,1.99999E+01,1.00333E+01,1.99995E+01,7.02116E-04,3.18849E-04,3.41082E-04,4.82013E-04,5.29625E-04,1.99993E+01,1.99997E+01,7.81489E+00,1.73147E+00,5.22254E+00,1.14571E+01,1.99988E+01,1.42641E-03,1.99970E+01,1.99995E+01,1.99955E+01,7.35366E+00,7.78922E+00,8.57736E+00,1.35429E+01,6.60424E+00,1.99998E+01,1.99981E+01,1.99944E+01,1.99997E+01,7.96753E+00,6.68443E+00,1.07529E+01,1.09947E+01,1.76224E+01,1.66186E-05,1.99991E+01,1.24234E+01,1.46716E+01,8.04241E+00,6.96345E+00,1.20863E+01,1.10437E+01,1.46402E+01,4.20585E-05,1.99990E+01,9.73853E+00,1.22405E+01,6.04387E+00,7.83196E+00,1.28775E+01,1.13753E+01,1.24647E+01,1.99994E+01,4.95872E+00,1.17946E+01,6.11143E+00,7.93591E+00,1.31467E+01,1.05024E+01,8.75922E+00,1.99998E+01,7.34439E-01,5.19371E+00,1.03024E+01,1.34907E+01,9.85130E+00};
		//double[] solutionVector=new double[] {[7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 14.0, 14.0, 10.0, 10.0, 10.0, 10.0, 14.0, 14.0, 14.0, 14.0, 10.0, 10.0, 10.0, 10.0, 14.0, 14.0, 14.0, 14.0, 14.0, 10.0, 10.0, 10.0, 10.0, 16.0, 16.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 16.0, 16.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 16.0, 16.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 16.0, 16.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 12.0, 5.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
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
	            if(i==2) {
	            	//intensityRealDose.add(intensityVoxel2Beams);
	            	//System.out.println(intensityVoxel2Beams);
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
   					pen += Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				}
            	
	        }
			score=pen;
			F+=w[i]*(pen/aux_index.size());
			System.out.println("organo "+String.valueOf(i)+": "+score);
			//System.out.println("organo: "+w[i]);
		}

		System.out.println("OF: "+F);
		return F;
	}

			
}
