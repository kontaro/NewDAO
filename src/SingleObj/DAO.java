package SingleObj;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;


import IMRT_DAO.*;
import gurobi.GRBException;


public class DAO {
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
    public int max_delta;                    // # DAO : Delta mximo para la variacin de intesidad por apertura
    public int max_iter;              		// # DAO : Cantidad mxima de iteraciones del algoritmo
    public int max_time;                     // # DAO : tiempo mximo de ejecucin del algoritmo
    public int seed;                         // # DAO : semilla para partir la busqueda local
    public int step_intensity;				// # DAO : Step size for aperture intensity (2)
    public Vector<int[][]> initial_aperture_shape;  // # DAO : Conjunto de formas de las aperturas segn restricciones (rango sup y rango inf)
    public Vector<int[][]> apertureShapes;   // # DAO: Aperturas segun restricciones del problema (en forma de matriz)
    public int [] zmax,zmin;
    public double[] dd;
    public int totalbmlt;
    public int[] selAngles;
    public double[] w;
    public int[] beamletAngles;
    public TreatmentPlan sol;
	
	    public DAO(){
	    	
	    }
		public DAO(String inputFile) throws IOException {
		
		//------ Lectura de archivos y creacion de DDM  y seteo de parametros------
		

		
	    readInputFile(inputFile);
	    beamletSet = new int[numAngles][4];
		beamletAngles = new int[numAngles];
	    String beamInfoDir = pathFile+"beamsInfo.txt";
	    //beamletSet = loadNumBixels(beamInfoDir);
	    beamletSet = loadNumBixelsTRT(beamInfoDir);
	    for(int i = 0;i<1;i++) {//con 14 haces todos los angulos
	    	//printBAC() ;
	    	for(int j = 0;j<1;++j) {// 30
	    		long time=System.currentTimeMillis();
			    setBeamlet();
			    int totalbmlt = totalBeamlet();
			    //....Generacion DDM
			    
			    DDM M; 
			    M = new DDM(o, initialBACs, pathFile,beamletAngles,totalbmlt);
			    
			    //Obtener dosis de cada organo
			    setDD();
			    
			    w = new double[]{0,1,1};
			    //-----Generar solucion Inicial------/
			    
			    InitialSolution newInitial=new InitialSolution(w,beamletSet,initialBACs,pathFile,M,o,dd,beamletAngles,init_intensity, 20, max_delta, max_iter, max_time, seed, 2, numAngles,numOrgans);
			   
			    //TreatmentPlan solution = newInitial.test4();
			    
			    TreatmentPlan solution = newInitial.test4gEUD();
			    solution.evaluationMethod=1;
			    solution.evaluationFunction=3;
			    solution.evaluateSolution();
			    
				double initialvalue=solution.singleObjectiveValue;
				Algorithm heuristica= new Algorithm(beamletSet,initialBACs,M,o,dd,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
				try {
					heuristica.LSPen(M, solution);
					//heuristica.LSmath(M, solution);
					//heuristica.VND(M, solution);
				} catch (GRBException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				solution.evaluateSolution();
				solution.scorePrint();
				serializeDataOut(solution,"tpPen-05-bac4.tpg");
				System.out.println("result: "+solution.singleObjectiveValue+" "+solution.getLowIntensityperBeam()+" "+solution.getLowIntensityAngle()+" "+(System.currentTimeMillis()-time)+" "+solution.beamOnTime()+" "+heuristica.iterations);	;
				
				ArrayList<ArrayList<Double>> Doses=DoseVolume(solution, M, o);
				doseprint(Doses,"tpPen-05-bac4",i,j);	
	    	}
	    	for(int h =0;h<initialBACs.length;h++) {
				initialBACs[h]=initialBACs[h]+5;
				System.out.print(initialBACs[h]+" ");
			}
	    	
	    }
	    

	    
	    //-----Seleccionar heuristica DAO-----///
	    selectionCriterion=0;
	    

	}
		
		
		public DAO(String inputFile,TreatmentPlan initial) throws IOException {
			
		//------ Lectura de archivos y creacion de DDM  y seteo de parametros------
		

		
	    readInputFile(inputFile);
	    beamletSet = new int[numAngles][4];
		beamletAngles = new int[numAngles];
	    String beamInfoDir = pathFile+"beamsInfo.txt";
	    beamletSet = loadNumBixels(beamInfoDir);
	    
	    for(int i = 0;i<1;i++) {//con 14 haces todos los angulos
	    	//printBAC() ;
	    	for(int j = 0;j<7;++j) {// 30
	    		long time=System.currentTimeMillis();
			    setBeamlet();
			    int totalbmlt = totalBeamlet();
			    //....Generacion DDM
			    DDM M; 
			    M = new DDM(o, initialBACs, pathFile,beamletAngles,totalbmlt);
			    
			    //Obtener dosis de cada organo
			    setDD();
			    
			    w = new double[]{5,1,1};
			    //-----Generar solucion Inicial------/
			    
			    InitialSolution newInitial=new InitialSolution(w,beamletSet,initialBACs,pathFile,M,o,dd,beamletAngles,init_intensity, 20, max_delta, max_iter, max_time, seed, 2, numAngles,numOrgans);
	
			    
			    TreatmentPlan solution = newInitial.test4();

			
				double initialvalue=solution.singleObjectiveValue;
				Algorithm heuristica= new Algorithm(beamletSet,initialBACs,M,o,dd,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
				try {
					heuristica.LSmath(M, solution);
					//heuristica.rVNS(M, solution);
				} catch (GRBException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				solution.evaluateSolution();
				solution.scorePrint();
				serializeDataOut(solution);
				System.out.println(solution.singleObjectiveValue+" "+solution.getLowIntensityperBeam()+" "+solution.getLowIntensityAngle()+" "+(System.currentTimeMillis()-time)+" "+solution.beamOnTime()+" "+heuristica.iterations);
				//ArrayList<ArrayList<Double>> Doses=DoseVolume(solution, M, o);
				//doseprint(Doses,"FMO_");
				solution.numIntensity();
				;
				
				
	    	}
			for(int h =0;h<initialBACs.length;h++) {
				initialBACs[h]=initialBACs[h]+5;
				System.out.print(initialBACs[h]+" ");
			}
			System.out.println();
	    }
	    
	    

	    
	    //-----Seleccionar heuristica DAO-----///
	    selectionCriterion=0;
	    

	}
		
		/**
		 * 
		 * @param inputFile direccion de archivo con parametros
		 * @param initialBACs Configuracion de BAC para DAO
		 * @param w Pesos de la funcion mono objetico [PTV,Bladder,Rectucm]
		 * @param ob Indica que funcion objetivo se va a utilizar
		 * @return
		 * @throws IOException
		 */
		
		public DAO(String inputFile,int initialBACs[],double[] w)  throws IOException {
			
			//------ Lectura de archivos y creacion de DDM  y seteo de parametros------
			

			
		    readInputFile(inputFile);
		    beamletSet = new int[numAngles][4];
			beamletAngles = new int[numAngles];
		    String beamInfoDir = pathFile+"beamsInfo.txt";
		    
		    //....Carga de beamlets segun BAC usado
		    this.initialBACs=initialBACs;
		    beamletSet = loadNumBixels(beamInfoDir);
		    setBeamlet();
		    int totalbmlt = totalBeamlet();
		    //....Generacion DDM
		    DDM M; 
		    M = new DDM(o, initialBACs, pathFile,beamletAngles,totalbmlt);
		    
		    //Obtener dosis de cada organo
		    setDD();
		    //Asignacion de pesos para la funcion objetivo
		    this.w = w;
		    
		    InitialSolution newInitial=new InitialSolution(w,beamletSet,initialBACs,pathFile,M,o,dd,beamletAngles,init_intensity, 20, max_delta, max_iter, max_time, seed, 2, numAngles,numOrgans);

		    
		    TreatmentPlan solution = newInitial.test4gEUD();
		    solution.evaluationMethod=1;
		    solution.evaluationFunction=3;
		    solution.evaluateSolution();
			Algorithm heuristica= new Algorithm(beamletSet,initialBACs,M,o,dd,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
			try {
				heuristica.LSPen(M, solution);
			} catch (GRBException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			sol=solution;
			//solution.evaluateSolutionPen(M, o, w, solution.intensity);

			//System.out.println(solution.singleObjectiveValue);
	    	
			



		}	
		
	public double[] returnFluenceMap() {
		return sol.intensity;
		
	}	
	public double returnobjectiveFunction() {
		return sol.singleObjectiveValue;
		
	}
	public void doseprint(ArrayList<ArrayList<Double>> Doses,String init,int bac, int iter){
		
		String fullName=init+"-"+bac+"-"+iter+"_";
		System.out.println("Working Directory = " + System.getProperty("user.dir"));
		PrintWriter writer;
		try {
			writer = new PrintWriter((fullName+"0.txt"), "UTF-8");
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
			writer = new PrintWriter((fullName+"1.txt"), "UTF-8");
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
			writer = new PrintWriter((fullName+"2.txt"), "UTF-8");
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
	public DDM getDDM() {
		return null;
	}
	 public static int[][] loadNumBixels(String beamsInfoDir) throws IOException{
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
	            for (int i=0;i<360;i++){
	                line=fileIn.readLine();
	                auxReader = line.split("\t");
	                beamlets[i][0]=(int) Double.parseDouble(auxReader[0]); //beamIndex
	                beamlets[i][1]=(int) Double.parseDouble(auxReader[1]); //numBeamlets
	                beamlets[i][2]= auxBeamletIndex;                       //firstBeamletIndex
	                beamlets[i][3]=(int) Double.parseDouble(auxReader[2]) - 1; //lastBeamletIndex
	                auxBeamletIndex = auxBeamletIndex + (int) Double.parseDouble(auxReader[1]);
	            }
	            fileIn.close();
	        }
	        return (beamlets);
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
	
	
	//----------------------------------------------------
	
	public  void setBeams(TreatmentPlan solution, String dir) throws IOException{
		totalbmlt=0;
		for(int i=0; i<numAngles; i++){
			solution.selAngles[i] = new Beam(beamletSet[initialBACs[i]][1],initialBACs[i], pathFile);
			totalbmlt=totalbmlt+solution.selAngles[i].beamlets;
		}
		solution.beamlets=totalbmlt;
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
			solution.beamlets=beamlets;
		}
	}

	public void readInputFile(String dir) throws IOException{
	        
    	String sp="\\s+";
        //String dir = "./inputFile.txt";
        File f = new File(dir);
        //System.out.println(f.getAbsolutePath());
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
        
        //Mximo delta para la variacin de intensidad en las aperturas
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
			
			//Solucin inicial completamente abierta
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
			
			//Solucin inicial completamente abierta
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
			
			//Solucin inicial completamente abierta
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
			
			//Solucin inicial completamente abierta
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
			
			//Soluci�n inicial con limites random
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
	
	
	//---------------- Funcion de evaluacion--------------------------------------
	public double evaluateSolution(TreatmentPlan solution, DDM M, Organs[] o, double w[]){
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
		//double[] solutionVector = new double[]{4,4,4,21,21,15,10,4,4,21,21,21,21,15,15,10,10,4,21,21,15,15,15,15,10,8,4,21,15,15,15,15,10,10,10,4,21,15,15,10,10,10,10,10,4,21,15,15,15,15,10,10,10,4,15,21,21,15,10,10,10,4,10,21,10,8,4,4,4,4,10,10,10,11,4,4,4,10,10,10,17,17,14,11,4,10,10,10,17,14,14,14,4,11,11,11,11,11,17,10,4,4,10,10,14,14,14,11,10,4,4,14,14,14,14,14,14,4,4,10,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,9,4,4,4,4,4,4,9,9,9,9,9,4,4,4,9,13,13,13,14,6,4,4,4,4,4,6,6,9,9,4,4,4,4,4,4,4,9,9,4,4,4,4,4,4,4,14,4,4,4,4,9,9,9,9,9,4,4,6,9,6,4,4,4,4,4,4,4,4,7,5,4,4,5,9,9,9,4,4,7,7,13,13,9,5,4,4,5,13,13,13,13,9,9,5,4,9,9,9,9,13,9,9,9,4,7,9,9,9,9,9,13,5,4,7,7,7,7,7,7,7,9,4,5,7,7,7,7,7,9,13,4,4,4,5,7,4,4,4,4,4,4,4,4,4,4,4,4,7,4,4,4,12,11,11,11,11,11,7,4,4,7,21,12,12,12,12,12,4,4,4,21,12,12,12,7,7,7,4,4,21,12,12,12,12,12,4,4,11,11,11,11,11,12,12,4,7,11,11,11,12};
		//double[] solutionVector = new double[]{10,20,20,14,20,19,11,10,20,5,6,7,16,19,18,14,9,20,10,4,5,14,17,16,11,10,20,12,4,5,13,17,15,9,10,20,13,5,4,13,18,12,8,11,19,14,7,4,14,17,9,10,12,20,9,4,17,18,8,14,10,20,20,20,5,14,13,20,20,19,11,18,20,5,10,7,10,15,20,15,3,15,6,3,12,17,20,9,2,16,4,3,9,20,20,6,3,19,2,4,1,5,17,17,0,0,17,0,4,0,0,17,20,1,1,19,1,0,1,0,13,0,1,19,0,4,0,7,11,4,0,20,15,3,0,0,0,9,10,6,20,0,1,3,8,10,19,7,20,6,3,6,8,9,20,5,9,20,8,1,4,5,6,17,6,7,13,8,1,2,3,2,10,5,8,16,10,2,3,3,1,5,4,4,20,9,7,9,3,4,3,5,5,3,2,4,7,4,6,2,7,4,0,0,20,9,20,20,20,20,12,1,20,20,20,20,19,18,9,10,17,20,18,15,14,18,17,7,13,16,19,14,15,14,20,20,6,17,15,20,14,16,12,16,19,5,18,19,20,20,19,13,14,19,3,19,20,20,20,20,12,7,20,4,17,19,20,6,0,17,1,7,18,2,0,0,3,4,0,0,0,0,0,0,0,2,0,0,11,17,19,10,11,15,0,10,10,8,17,18,10,13,14,11,2,8,14,20,20,13,13,10,12,0,4,14,19,14,11,8,8,12,0,17,19,9,8,7,7,11,20,19,4,5,6,6};
		//double[] solutionVector=new double[] {10.31,20.00,20.00,14.49,20.00,19.37,11.44,10.08,20.00,5.23,6.49,7.40,16.28,18.54,17.56,13.60,9.32,20.00,10.41,3.67,4.82,13.81,16.93,15.60,10.51,9.81,20.00,11.63,4.02,4.74,12.96,16.99,14.88,8.83,9.79,20.00,13.07,4.90,3.71,13.46,17.62,12.28,8.30,11.12,19.34,13.87,7.48,4.13,13.53,16.78,9.34,9.68,12.06,20.00,8.51,4.18,17.12,18.09,7.79,14.35,10.17,20.00,20.00,20.00,4.55,14.37,13.27,20.00,20.00,19.11,10.89,18.46,20.00,4.95,10.04,6.66,9.82,15.16,19.87,14.52,3.23,14.87,6.21,3.34,11.99,16.70,20.00,8.77,2.48,16.18,3.84,2.63,9.30,20.00,20.00,5.58,2.86,19.32,2.23,3.68,0.68,5.38,17.18,17.28,0.00,0.00,16.50,0.00,4.28,0.00,0.00,17.06,20.00,1.04,0.58,19.19,0.92,0.00,0.54,0.45,13.30,0.00,1.39,18.92,0.00,4.37,0.00,7.01,11.15,4.37,0.00,20.00,15.45,2.62,0.00,0.00,0.00,8.54,9.78,5.60,20.00,0.00,1.43,3.35,8.33,10.28,19.01,6.61,20.00,5.81,2.98,6.48,8.28,8.88,20.00,5.25,8.81,20.00,8.45,1.12,3.79,5.25,5.60,16.52,6.01,7.16,12.88,7.73,0.68,1.61,3.07,1.61,10.29,5.49,7.75,16.13,9.73,1.64,3.13,2.94,1.12,5.49,3.85,3.76,20.00,8.51,6.98,9.28,2.81,3.76,3.33,5.10,5.30,2.95,2.44,3.95,7.35,3.54,5.58,1.85,6.66,3.53,0.00,0.00,20.00,9.14,19.85,20.00,20.00,20.00,12.22,1.47,20.00,20.00,20.00,20.00,18.61,17.75,8.77,9.74,16.53,19.93,18.08,15.31,13.80,18.36,16.93,6.52,13.03,16.03,19.12,13.84,14.91,13.61,20.00,19.84,6.39,16.74,15.23,19.66,13.69,15.69,12.40,15.65,18.58,4.84,18.19,18.92,20.00,20.00,18.85,12.91,13.85,19.30,3.44,19.49,20.00,20.00,20.00,20.00,12.08,7.44,19.75,3.63,17.48,18.65,20.00,6.44,0.00,17.11,1.27,7.16,18.34,1.94,0.00,0.00,2.97,3.64,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.83,0.00,0.00,10.85,17.47,18.61,9.50,11.06,14.53,0.00,9.56,9.64,8.43,17.10,17.73,10.13,12.52,13.76,10.97,1.70,7.53,13.55,20.00,20.00,13.21,12.79,9.77,12.10,0.00,3.59,14.10,19.19,13.58,10.81,8.01,8.19,11.79,0.00,17.35,19.34,9.14,8.15,6.98,7.11,11.09,20.00,18.97,3.56,5.08,5.88,5.87};
		//double[] solutionVector=new double[] {20.0,20.0,0.0,12.0,8.0,12.0,12.0,4.0,0.0,12.0,8.0,4.0,8.0,4.0,8.0,8.0,0.0,20.0,16.0,4.0,0.0,4.0,0.0,4.0,4.0,0.0,20.0,16.0,4.0,0.0,4.0,0.0,4.0,4.0,0.0,20.0,20.0,4.0,0.0,8.0,4.0,4.0,4.0,0.0,20.0,20.0,8.0,0.0,8.0,4.0,4.0,4.0,4.0,20.0,12.0,4.0,12.0,8.0,8.0,12.0,4.0,20.0,16.0,16.0,16.0,20.0,20.0,0.0,20.0,20.0,8.0,20.0,16.0,20.0,20.0,16.0,0.0,16.0,20.0,12.0,16.0,20.0,16.0,12.0,4.0,16.0,20.0,8.0,16.0,16.0,12.0,8.0,4.0,20.0,20.0,8.0,12.0,16.0,8.0,8.0,20.0,4.0,12.0,16.0,4.0,8.0,8.0,0.0,4.0,0.0,0.0,12.0,20.0,4.0,4.0,4.0,0.0,0.0,0.0,0.0,16.0,4.0,4.0,4.0,0.0,4.0,0.0,0.0,4.0,4.0,0.0,4.0,4.0,8.0,0.0,0.0,0.0,4.0,0.0,0.0,20.0,0.0,4.0,8.0,12.0,8.0,16.0,20.0,20.0,12.0,8.0,16.0,20.0,16.0,20.0,20.0,20.0,20.0,12.0,4.0,12.0,20.0,12.0,16.0,20.0,20.0,16.0,12.0,8.0,12.0,20.0,12.0,16.0,12.0,12.0,16.0,12.0,8.0,16.0,20.0,12.0,12.0,12.0,8.0,20.0,12.0,12.0,20.0,20.0,12.0,12.0,16.0,8.0,20.0,12.0,12.0,16.0,8.0,16.0,8.0,16.0,12.0,0.0,0.0,20.0,16.0,20.0,20.0,20.0,20.0,16.0,8.0,20.0,20.0,20.0,20.0,16.0,16.0,16.0,16.0,20.0,20.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,16.0,20.0,12.0,16.0,16.0,20.0,20.0,16.0,20.0,20.0,20.0,12.0,16.0,12.0,16.0,20.0,16.0,20.0,20.0,20.0,20.0,20.0,16.0,20.0,20.0,16.0,20.0,20.0,20.0,20.0,20.0,16.0,12.0,20.0,16.0,16.0,16.0,20.0,0.0,0.0,12.0,4.0,12.0,20.0,0.0,0.0,0.0,4.0,4.0,0.0,4.0,4.0,0.0,0.0,0.0,0.0,8.0,0.0,0.0,8.0,12.0,8.0,4.0,4.0,4.0,4.0,8.0,4.0,4.0,12.0,8.0,4.0,4.0,8.0,20.0,0.0,4.0,12.0,16.0,16.0,12.0,8.0,8.0,20.0,20.0,0.0,12.0,20.0,16.0,8.0,8.0,4.0,16.0,0.0,16.0,20.0,16.0,8.0,4.0,4.0,20.0,20.0,20.0,16.0,8.0,8.0,8.0,};
		//double[] solutionVector=new double[] {6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814094235, 6.809247814096439, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814094235, 6.809247814094235, 6.809247814094235, 6.809247814094235, 6.809247814096439, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814094235, 12.76375812908914, 6.809247814094235, 6.809247814094235, 6.809247814096439, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814094235, 7.663985313185287, 7.663985313185287, 7.663985313185287, 7.663985313187491, 7.663985313179215, 7.663985313179215, 6.809247814088163, 6.809247814088163, 6.807326388063126, 6.807326388063126, 6.807326388063126, 6.807326388063126, 6.807326388065331, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388063126, 7.662063887154178, 7.662063887154178, 7.662063887154178, 7.662063887156383, 7.662063887148107, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388063126, 12.759915277026924, 6.807326388063126, 6.807326388065331, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388063126, 6.807326388065331, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322238302, 13.616766322224539, 13.616766322224539, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322238302, 13.616766322224539, 13.616766322224539, 13.616766322224539, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322238302, 13.616766322224539, 13.616766322224539, 13.616766322224539, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322238302, 13.616766322224539, 13.616766322224539, 13.616766322224539, 13.616766322224539, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290232553, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290232553, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290232553, 13.86658529021879, 13.86658529021879, 9.049069919423301, 9.049069919423301, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 16.523844075238834, 16.523844075238834, 16.52384407524183, 12.707081704052866, 8.516896861836218, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.52384407524183, 12.707081704052866, 12.707081704052866, 12.707081704052866, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.52384407524183, 16.52384407524183, 16.52384407524183, 12.707081704052866, 12.707081704052866, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.52384407524183, 12.707081704052866, 12.707081704052866, 12.707081704052866, 12.707081704052866, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119484485, 11.523065748295522, 11.523065748295522, 11.523065748295522, 11.523065748295522, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119484485, 11.523065748295522, 11.523065748295522, 7.3328809060788736, 7.3328809060788736, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119484485, 11.523065748295522, 11.523065748295522, 11.523065748295522, 11.523065748295522, 15.339828119484485, 11.523065748295522, 4.190184842219642, 4.190184842219642, 4.190184842219642, 11.523065748295522, 11.523065748295522, 11.523065748295522, 11.523065748295522, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.315450570383323, 8.315450570385886, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.315450570383323, 8.315450570383323, 8.315450570383323, 8.315450570385886, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.315450570383323, 8.315450570383323, 8.315450570383323, 8.315450570383323, 8.315450570385886, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.87868622158971, 11.87868622158971, 11.87868622158971, 7.064866331967789, 7.064866331947883, 7.064866331947883, 7.064866331947883, 7.064866331947883, 11.87868622158971, 11.87868622158971, 11.87868622158971, 11.87868622158971, 11.878686221592748, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.878686221572842, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345};
		//double[] solutionVector=new double[] {[7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 14.0, 14.0, 10.0, 10.0, 10.0, 10.0, 14.0, 14.0, 14.0, 14.0, 10.0, 10.0, 10.0, 10.0, 14.0, 14.0, 14.0, 14.0, 14.0, 10.0, 10.0, 10.0, 10.0, 16.0, 16.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 16.0, 16.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 16.0, 16.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 16.0, 16.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 16.0, 16.0, 16.0, 12.0, 12.0, 12.0, 12.0, 12.0, 5.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
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
			//System.out.println("organo: "+w[i]);
		}
		solution.singleObjectiveValue=F;
		//System.out.println("Solucion: "+F);
		return F;
	}
	
	/**La funcion evaluada pero penalizada
	 * 
	 * @param M
	 * @param o
	 * @param w
	 * @param solutionVector
	 * @return
	 */
			
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
			//solution.scores[i]=score;
			F+=pen/aux_index.size();
			//System.out.println("organo: "+solution.scores[i]);
		}
		//solution.singleObjectiveValue=F;
		System.out.println("Solucion: "+F);
		return F;
	}
	
	
	
	
	/**
	 * 
	 * Todo:
	 * agregar variable de nombre de archivo con bac
	 */
	public static void serializeDataOut(TreatmentPlan solutionGenerated)throws IOException{
	    String fileName= "treatment511.tpg";
	    FileOutputStream fos = new FileOutputStream(fileName);
	    ObjectOutputStream oos = new ObjectOutputStream(fos);
	    oos.writeObject(solutionGenerated);
	    oos.close();
	}
	/**
	 * 
	 * Funcion que guarda el objeto solution para ser cargado despues
	 */
	public static void serializeDataOut(TreatmentPlan solutionGenerated, String fileName)throws IOException{
	    
	    FileOutputStream fos = new FileOutputStream(fileName);
	    ObjectOutputStream oos = new ObjectOutputStream(fos);
	    oos.writeObject(solutionGenerated);
	    oos.close();
	}
	public static TreatmentPlan serializeDataIn() throws IOException, ClassNotFoundException{
	   String fileName= "treatment.tpg";
	   FileInputStream fin = new FileInputStream(fileName);
	   ObjectInputStream ois = new ObjectInputStream(fin);
	   TreatmentPlan saveTreatmentPlan= (TreatmentPlan) ois.readObject();
	   ois.close();
	   return saveTreatmentPlan;
	}
	public static TreatmentPlan serializeDataIn(String fileName) throws IOException, ClassNotFoundException{
		   
	   FileInputStream fin = new FileInputStream(fileName);
	   ObjectInputStream ois = new ObjectInputStream(fin);
	   TreatmentPlan saveTreatmentPlan= (TreatmentPlan) ois.readObject();
	   ois.close();
	   return saveTreatmentPlan;
	}
	
	
	public ArrayList<ArrayList<Double>> DoseVolume(TreatmentPlan solution, DDM M, Organs[] o){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		
		ArrayList<ArrayList<Double>> organ=new ArrayList<ArrayList<Double>>();
		double[] solutionVector = solution.intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		//double[] solutionVector=new double[] {4.603306596809836E-10, 4.1522815865467795, 4.180678677244826, 4.531710564128731E-11, 3.983994950386257E-12, 10.284709871426173, 12.820877644391064, 8.03881222582465, 9.92179932065358E-12, 2.6179339841827396, 13.119283914517439, 16.586771860443562, 10.09920393986466, 25.567088012609336, 25.587036791035846, 0.7572689661837781, 2.4802033377687382, 12.406384310009138, 16.484953209912852, 4.725601544672685, 1.8734994277032416E-10, 19.80195329382713, 20.619885463206717, 0.8280479389577354, 0.41990102286193376, 12.40468752700521, 12.747076861731564, 1.2815497388152461, 2.6194154307349034E-11, 18.403813154770877, 17.77017882719856, 4.5382239965917835E-12, 5.0432802541544965E-12, 13.325668724626066, 9.746877584118822, 0.6532332817197885, 3.362748376436782E-11, 18.479107000161715, 13.282513669990516, 0.993800888658923, 2.272797716123517E-10, 14.371434933024625, 9.243131254418175, 2.143611297767546, 0.6146522494139087, 18.554076188795456, 9.767191814464367, 3.141741342315267, 0.5225721398106938, 14.689668519661586, 11.336857260548596, 2.591928686824542, 12.702077619109915, 1.505354124372825E-10, 3.654794364915715, 0.6995811382254522, 13.87562522859049, 13.010114218406848, 1.4267922505403765, 1.3445092287980867E-10, 0.5992413713984279, 13.886352072902211, 12.746641433325735, 3.820773660235344E-12, 3.1555788169473194E-11, 13.728966518613705, 9.077678587960634, 1.9860596169275874E-11, 11.777472311303692, 1.2003422650688997E-11, 4.4378665830939274E-10, 8.147005940018076, 18.86766414773367, 2.42834714447056, 22.992195263256516, 17.822468424911648, 14.0789697986716, 4.286479660823125, 7.093611067408727E-12, 6.918446272992865, 19.883906108769118, 18.978797499960354, 4.456010630712492E-10, 5.0077399298800954E-11, 0.4570056718388383, 10.648576005563065, 20.053018414206548, 18.761044403388933, 1.0906474305990947E-11, 42.34609365798736, 4.644198536446167, 2.7295013650621427E-11, 3.0213114341335645, 15.106648652044598, 20.46903328574452, 18.361282795329668, 7.789868897546472, 43.80135714588408, 4.5492545883296214E-11, 1.7589440670096994, 7.473552202774985, 19.441342763670548, 20.794030118186704, 16.722570206228934, 15.305192659306215, 41.28729762887966, 3.083175711172887E-11, 1.025604188791615E-11, 10.792537048956328, 23.627904782274573, 21.350469793366887, 16.058188986523678, 21.39244645752887, 38.90245135176751, 25.17189537760626, 22.106430530667275, 14.09087028032827, 24.34241410486692, 41.075480854187305, 40.31664845802923, 20.49473080862755, 10.685986343943453, 24.35542961972121, 39.06581728064339, 23.233797169797167, 6.758232707896598, 25.48488503102567, 39.64391192740875, 48.06070368723197, 1.4075753476130723E-10, 24.268323187565297, 37.3153652817524, 7.880110986335836, 11.894439255818552, 70.27820883227845, 8.112045725745546, 13.722491922589947, 1.4047392579004107E-11, 19.376512258811992, 53.74696151878255, 53.28908074780916, 5.821141373312123, 10.653743600891367, 3.958950261692384, 4.48586453443563, 0.39892767280863495, 46.309215309496125, 22.718738252224473, 35.17268026725934, 4.576264286488712, 9.294090519364786, 5.178394456893013, 0.15773635389379204, 0.44878301950184146, 29.29366740548964, 20.64086885266853, 26.115334869272047, 3.9225504992244598, 8.619732335948278, 5.147494085885313, 9.325078614859089E-13, 9.387503571847502E-12, 17.71187802927547, 18.91024561780447, 20.161971593614656, 2.3904847021334805, 8.325370994860174, 5.178170511008548, 0.8330055974022106, 0.0387377940397586, 4.620180096988544, 18.925180474158722, 17.184370916123633, 1.6568795912970186E-11, 6.456214139868457, 5.4310958585770095, 2.9805566662674803, 2.4761288156856703E-11, 21.848496190358517, 10.556408473706695, 7.33415051640792E-12, 5.007386436049916, 3.7447657418516207, 4.282969183557636, 2.462001371582888, 3.8009212814011786, 3.6520818112540985, 2.5767720752904064, 1.1226362138695845, 1.0872688199434506, 5.053186246752819, 0.9426932891485297, 0.19469093978161253, 1.1608230517596583E-11, 11.917400290553239, 1.0851728568829095E-10, 1.3431266461782982, 4.316622428887439E-4, 1.866056850542748, 4.66389020588248, 6.355914181898651, 35.22706604201244, 7.164244095746165E-11, 11.789583605065706, 3.863559529224403, 4.4725891113590455E-11, 9.332871541041382E-12, 7.828763882435012, 0.738805066745458, 1.3521104971418592, 3.267662976144205, 5.306705167616092, 5.81118704810603E-12, 4.083533586988921, 5.200188107144168, 3.493678122665157, 2.532207086319474, 0.44807797587347653, 2.356938368865201, 3.906850469368954, 1.638349926396288E-11, 8.580038958608831, 2.7814433676809527, 1.3272717518394072E-11, 1.7780316539615943, 3.727061634547518, 22.079556453717267, 2.1299348032732188E-11, 0.23882561007422548, 11.167404777927375, 3.2104810352356394, 6.328781400745852E-12, 1.831560763091865, 4.809805742550826, 28.06511242213768, 1.0853901036587765, 1.4400813075243544, 15.695215466021954, 3.9808770914376654, 0.3991905287527887, 1.0297843410462117, 6.856104218405361, 35.14032805190652, 2.0495808076800137, 2.419231802657311, 19.44709887738912, 5.067173985431579, 0.9721827398260882, 9.487178396967178E-13, 5.768438502778807, 53.70417571964857, 1.5977206574095795E-11, 3.037576983418051E-11, 20.556660574379542, 4.073387577135681, 0.3822902064659835, 0.20246712358648225, 5.473463880504645, 26.37537258535767, 19.51764618696966, 3.622638739647977, 1.3156641620498037, 4.005729985202039, 6.87813457969773, 6.659379698942805, 5.6791151193025495, 11.206735773581935, 24.442033866708286, 12.050700037069506, 13.33801182613702, 5.255328332511672, 10.315062375826383, 29.374593741435792, 20.74896780674686, 11.559023858762712, 2.7382871719113882E-11, 22.14114692808186, 29.39529987929727, 20.490351881055794, 13.470096846661798, 3.028543995297505E-11, 28.185656773217218, 24.59418016361163, 21.290132978540782, 12.193940143516276, 2.0402563635764312E-11, 1.0601569839385626E-11, 1.2322480147936392E-11, 16.563293809266515, 28.670403228023147, 22.095389656310108, 22.150882470222676, 10.729056479088984, 3.5754588935436564E-10, 10.571596565179895, 7.3279349793428565, 19.70311945630245, 26.50443057166728, 18.408477785370305, 20.1786832517126, 6.92756724493922, 16.9323959989572, 9.768498128515617, 6.564080650679913, 21.207008881835883, 25.91978083839274, 15.339068766556297, 9.928986122389496E-12, 8.886119920526898, 18.261986231497986, 7.286940506674815, 4.688046106070827, 23.556793829923304, 26.053119443941075, 9.129285199393365, 47.08862422353665, 2.5407542395374867E-11, 1.0868513797171119E-11, 26.7258696938149, 24.842977146171418, 4.475074549696293, 27.548927844587627, 27.75185379249043, 26.450846443344503, 0.11422923155235745, 1.5165050636951056, 9.090487395259379E-11};
	
		System.out.println("Largo soluci�n :"+solutionVector.length);
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
			organ.add(Doses);

		}

		return organ;
	}
	
	public double gradientSolution(TreatmentPlan solution, DDM M, Organs[] o, double w[]){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams,intensityVoxel2BeamsIntensidad;
		Integer key, beam;
		
		double[] solutionVector = solution.intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		//double[] vectorAperturas = solution.intesities;
		double[] vectorAperturas = new double[]{4,4,4,21,21,15,10,4,4,21,21,21,21,15,15,10,10,4,21,21,15,15,15,15,10,8,4,21,15,15,15,15,10,10,10,4,21,15,15,10,10,10,10,10,4,21,15,15,15,15,10,10,10,4,15,21,21,15,10,10,10,4,10,21,10,8,4,4,4,4,10,10,10,11,4,4,4,10,10,10,17,17,14,11,4,10,10,10,17,14,14,14,4,11,11,11,11,11,17,10,4,4,10,10,14,14,14,11,10,4,4,14,14,14,14,14,14,4,4,10,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,9,4,4,4,4,4,4,9,9,9,9,9,4,4,4,9,13,13,13,14,6,4,4,4,4,4,6,6,9,9,4,4,4,4,4,4,4,9,9,4,4,4,4,4,4,4,14,4,4,4,4,9,9,9,9,9,4,4,6,9,6,4,4,4,4,4,4,4,4,7,5,4,4,5,9,9,9,4,4,7,7,13,13,9,5,4,4,5,13,13,13,13,9,9,5,4,9,9,9,9,13,9,9,9,4,7,9,9,9,9,9,13,5,4,7,7,7,7,7,7,7,9,4,5,7,7,7,7,7,9,13,4,4,4,5,7,4,4,4,4,4,4,4,4,4,4,4,4,7,4,4,4,12,11,11,11,11,11,7,4,4,7,21,12,12,12,12,12,4,4,4,21,12,12,12,7,7,7,4,4,21,12,12,12,12,12,4,4,11,11,11,11,11,12,12,4,7,11,11,11,12};
		//double[] solutionVector = new double[]{10,20,20,14,20,19,11,10,20,5,6,7,16,19,18,14,9,20,10,4,5,14,17,16,11,10,20,12,4,5,13,17,15,9,10,20,13,5,4,13,18,12,8,11,19,14,7,4,14,17,9,10,12,20,9,4,17,18,8,14,10,20,20,20,5,14,13,20,20,19,11,18,20,5,10,7,10,15,20,15,3,15,6,3,12,17,20,9,2,16,4,3,9,20,20,6,3,19,2,4,1,5,17,17,0,0,17,0,4,0,0,17,20,1,1,19,1,0,1,0,13,0,1,19,0,4,0,7,11,4,0,20,15,3,0,0,0,9,10,6,20,0,1,3,8,10,19,7,20,6,3,6,8,9,20,5,9,20,8,1,4,5,6,17,6,7,13,8,1,2,3,2,10,5,8,16,10,2,3,3,1,5,4,4,20,9,7,9,3,4,3,5,5,3,2,4,7,4,6,2,7,4,0,0,20,9,20,20,20,20,12,1,20,20,20,20,19,18,9,10,17,20,18,15,14,18,17,7,13,16,19,14,15,14,20,20,6,17,15,20,14,16,12,16,19,5,18,19,20,20,19,13,14,19,3,19,20,20,20,20,12,7,20,4,17,19,20,6,0,17,1,7,18,2,0,0,3,4,0,0,0,0,0,0,0,2,0,0,11,17,19,10,11,15,0,10,10,8,17,18,10,13,14,11,2,8,14,20,20,13,13,10,12,0,4,14,19,14,11,8,8,12,0,17,19,9,8,7,7,11,20,19,4,5,6,6};
		//double[] solutionVector=new double[] {10.31,20.00,20.00,14.49,20.00,19.37,11.44,10.08,20.00,5.23,6.49,7.40,16.28,18.54,17.56,13.60,9.32,20.00,10.41,3.67,4.82,13.81,16.93,15.60,10.51,9.81,20.00,11.63,4.02,4.74,12.96,16.99,14.88,8.83,9.79,20.00,13.07,4.90,3.71,13.46,17.62,12.28,8.30,11.12,19.34,13.87,7.48,4.13,13.53,16.78,9.34,9.68,12.06,20.00,8.51,4.18,17.12,18.09,7.79,14.35,10.17,20.00,20.00,20.00,4.55,14.37,13.27,20.00,20.00,19.11,10.89,18.46,20.00,4.95,10.04,6.66,9.82,15.16,19.87,14.52,3.23,14.87,6.21,3.34,11.99,16.70,20.00,8.77,2.48,16.18,3.84,2.63,9.30,20.00,20.00,5.58,2.86,19.32,2.23,3.68,0.68,5.38,17.18,17.28,0.00,0.00,16.50,0.00,4.28,0.00,0.00,17.06,20.00,1.04,0.58,19.19,0.92,0.00,0.54,0.45,13.30,0.00,1.39,18.92,0.00,4.37,0.00,7.01,11.15,4.37,0.00,20.00,15.45,2.62,0.00,0.00,0.00,8.54,9.78,5.60,20.00,0.00,1.43,3.35,8.33,10.28,19.01,6.61,20.00,5.81,2.98,6.48,8.28,8.88,20.00,5.25,8.81,20.00,8.45,1.12,3.79,5.25,5.60,16.52,6.01,7.16,12.88,7.73,0.68,1.61,3.07,1.61,10.29,5.49,7.75,16.13,9.73,1.64,3.13,2.94,1.12,5.49,3.85,3.76,20.00,8.51,6.98,9.28,2.81,3.76,3.33,5.10,5.30,2.95,2.44,3.95,7.35,3.54,5.58,1.85,6.66,3.53,0.00,0.00,20.00,9.14,19.85,20.00,20.00,20.00,12.22,1.47,20.00,20.00,20.00,20.00,18.61,17.75,8.77,9.74,16.53,19.93,18.08,15.31,13.80,18.36,16.93,6.52,13.03,16.03,19.12,13.84,14.91,13.61,20.00,19.84,6.39,16.74,15.23,19.66,13.69,15.69,12.40,15.65,18.58,4.84,18.19,18.92,20.00,20.00,18.85,12.91,13.85,19.30,3.44,19.49,20.00,20.00,20.00,20.00,12.08,7.44,19.75,3.63,17.48,18.65,20.00,6.44,0.00,17.11,1.27,7.16,18.34,1.94,0.00,0.00,2.97,3.64,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.83,0.00,0.00,10.85,17.47,18.61,9.50,11.06,14.53,0.00,9.56,9.64,8.43,17.10,17.73,10.13,12.52,13.76,10.97,1.70,7.53,13.55,20.00,20.00,13.21,12.79,9.77,12.10,0.00,3.59,14.10,19.19,13.58,10.81,8.01,8.19,11.79,0.00,17.35,19.34,9.14,8.15,6.98,7.11,11.09,20.00,18.97,3.56,5.08,5.88,5.87};
		//double[] solutionVector=new double[] {10.3097,20,20,14.4913,20,19.3684,11.4378,10.075,20,5.2253,6.49239,7.39881,16.2819,18.5441,17.563,13.6014,9.32086,20,10.4118,3.67387,4.81576,13.8092,16.9314,15.6042,10.5102,9.80699,20,11.6295,4.02062,4.74281,12.9592,16.9919,14.8824,8.83135,9.78841,20,13.0725,4.89888,3.70851,13.4564,17.6188,12.2771,8.30026,11.1225,19.3446,13.8665,7.48044,4.1291,13.5259,16.7837,9.34001,9.68003,12.0576,20,8.5066,4.17597,17.1152,18.0933,7.79289,14.3488,10.1703,20,20,20,4.55037,14.3681,13.2736,20,20,19.1062,10.8932,18.4637,20,4.95187,10.0398,6.66412,9.82459,15.1614,19.8728,14.5203,3.22511,14.8658,6.20638,3.34013,11.9883,16.6987,19.9997,8.77295,2.47631,16.1838,3.8422,2.63254,9.29702,20,20,5.58084,2.86308,19.3228,2.23349,3.67659,0.677008,5.37564,17.1813,17.2778,5.34E-10,1.04E-09,16.5044,1.38E-10,4.27637,6.97E-10,1.66E-10,17.0639,20,1.04173,0.576239,19.1878,0.924726,2.55E-10,0.541289,0.450083,13.3042,4.80E-10,1.39492,18.915,2.08E-10,4.37446,5.83E-10,7.01266,11.153,4.37109,8.41E-10,20,15.4508,2.61777,6.29E-10,7.45E-10,4.03E-10,8.54297,9.77775,5.60274,20,6.33E-10,1.43363,3.35208,8.32635,10.2826,19.0077,6.61247,20,5.80767,2.97914,6.47571,8.28413,8.88179,20,5.2468,8.80924,20,8.44605,1.12493,3.78676,5.25264,5.59573,16.5185,6.01393,7.15919,12.8792,7.73481,0.680547,1.60642,3.06524,1.60948,10.2937,5.48907,7.74604,16.1253,9.73036,1.64011,3.12747,2.93724,1.11502,5.48732,3.8463,3.75988,20,8.51367,6.97555,9.27777,2.80921,3.75663,3.32548,5.10275,5.30177,2.95004,2.43942,3.94816,7.35041,3.54059,5.58316,1.85172,6.65664,3.52539,3.64E-09,3.36E-09,20,9.14021,19.854,20,20,20,12.2175,1.4684,20,20,20,20,18.6082,17.7476,8.7686,9.73554,16.5349,19.9298,18.0823,15.3061,13.8026,18.3586,16.931,6.52291,13.0256,16.0252,19.1165,13.8403,14.9063,13.6087,20,19.8418,6.39398,16.7429,15.234,19.6629,13.6932,15.6869,12.3977,15.6517,18.5804,4.84421,18.1889,18.918,20,20,18.8455,12.9068,13.8491,19.3014,3.44182,19.4928,20,20,20,20,12.0805,7.44039,19.7482,3.63102,17.4829,18.6455,20,6.43969,3.12E-09,17.1054,1.26997,7.15592,18.3381,1.9401,5.58E-10,2.14E-10,2.96576,3.64303,1.14E-09,1.55E-09,2.63E-09,9.97E-11,1.23E-10,1.79E-10,5.91E-11,1.82695,4.63E-10,5.50E-10,10.8482,17.4672,18.6106,9.50046,11.0649,14.5312,1.64E-10,9.56407,9.63549,8.42675,17.102,17.7268,10.1277,12.5151,13.7564,10.9661,1.69777,7.53269,13.5452,20,20,13.2125,12.7909,9.77117,12.097,5.81E-10,3.5894,14.1014,19.1897,13.5777,10.805,8.00704,8.19408,11.7871,5.21E-09,17.3491,19.3405,9.136,8.14758,6.98342,7.10701,11.0933,20,18.9689,3.55873,5.08056,5.87941,5.86829	};
		//double[] solutionVector=new double[] {6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814094235, 6.809247814096439, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814094235, 6.809247814094235, 6.809247814094235, 6.809247814094235, 6.809247814096439, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814094235, 12.76375812908914, 6.809247814094235, 6.809247814094235, 6.809247814096439, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814088163, 6.809247814094235, 7.663985313185287, 7.663985313185287, 7.663985313185287, 7.663985313187491, 7.663985313179215, 7.663985313179215, 6.809247814088163, 6.809247814088163, 6.807326388063126, 6.807326388063126, 6.807326388063126, 6.807326388063126, 6.807326388065331, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388063126, 7.662063887154178, 7.662063887154178, 7.662063887154178, 7.662063887156383, 7.662063887148107, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388063126, 12.759915277026924, 6.807326388063126, 6.807326388065331, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388063126, 6.807326388065331, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 6.807326388057055, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322238302, 13.616766322224539, 13.616766322224539, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322238302, 13.616766322224539, 13.616766322224539, 13.616766322224539, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322238302, 13.616766322224539, 13.616766322224539, 13.616766322224539, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322235602, 13.616766322238302, 13.616766322224539, 13.616766322224539, 13.616766322224539, 13.616766322224539, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290232553, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290232553, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.866585290229853, 13.866585290229853, 13.866585290229853, 13.866585290232553, 13.86658529021879, 13.86658529021879, 9.049069919423301, 9.049069919423301, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 13.86658529021879, 16.523844075238834, 16.523844075238834, 16.52384407524183, 12.707081704052866, 8.516896861836218, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.52384407524183, 12.707081704052866, 12.707081704052866, 12.707081704052866, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.52384407524183, 16.52384407524183, 16.52384407524183, 12.707081704052866, 12.707081704052866, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.523844075238834, 16.52384407524183, 12.707081704052866, 12.707081704052866, 12.707081704052866, 12.707081704052866, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119484485, 11.523065748295522, 11.523065748295522, 11.523065748295522, 11.523065748295522, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119484485, 11.523065748295522, 11.523065748295522, 7.3328809060788736, 7.3328809060788736, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119481492, 15.339828119484485, 11.523065748295522, 11.523065748295522, 11.523065748295522, 11.523065748295522, 15.339828119484485, 11.523065748295522, 4.190184842219642, 4.190184842219642, 4.190184842219642, 11.523065748295522, 11.523065748295522, 11.523065748295522, 11.523065748295522, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.315450570383323, 8.315450570385886, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.315450570383323, 8.315450570383323, 8.315450570383323, 8.315450570385886, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.315450570383323, 8.315450570383323, 8.315450570383323, 8.315450570383323, 8.315450570385886, 8.31545057031774, 8.31545057031774, 8.31545057031774, 8.31545057031774, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008783044, 9.192655008848627, 9.192655008848627, 9.19265500885119, 9.192655008783044, 9.192655008783044, 9.192655008783044, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.87868622158971, 11.87868622158971, 11.87868622158971, 7.064866331967789, 7.064866331947883, 7.064866331947883, 7.064866331947883, 7.064866331947883, 11.87868622158971, 11.87868622158971, 11.87868622158971, 11.87868622158971, 11.878686221592748, 11.878686221572842, 11.878686221572842, 11.878686221572842, 11.878686221572842, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345, 13.538909175635345, 13.538909175652211, 13.538909175652211, 13.538909175652211, 13.538909175655249, 13.538909175635345, 13.538909175635345};
		//double[] solutionVector=new double[] {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,67.0,67.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,67.0,67.0,0.0,0.0,0.0,0.0,0.0,0.0,67.0,67.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,67.0,67.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,54.0,54.0,0.0,0.0,0.0,0.0,0.0,0.0,54.0,54.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,13.0,13.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,13.0,13.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,13.0,13.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,13.0,13.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,13.0,13.0,0.0,0.0,54.0,54.0,0.0,0.0,0.0,13.0,13.0,0.0,0.0,0.0,53.0,30.0,30.0,0.0,0.0,0.0,30.0,83.0,83.0,83.0,83.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,80.0,0.0,0.0,0.0,0.0,80.0,83.0,83.0,80.0,0.0,0.0,0.0,0.0,0.0,80.0,83.0,83.0,80.0,0.0,0.0,0.0,0.0,0.0,80.0,83.0,83.0,80.0,0.0,0.0,0.0,0.0,0.0,80.0,83.0,83.0,80.0,0.0,0.0,0.0,0.0,0.0,80.0,83.0,83.0,80.0,0.0,0.0,3.0,83.0,83.0,83.0,83.0,0.0,75.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,75.0,0.0,0.0,0.0,0.0,75.0,75.0,75.0,75.0,0.0,0.0,0.0,0.0,0.0,75.0,75.0,75.0,75.0,0.0,0.0,0.0,0.0,0.0,75.0,75.0,75.0,75.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,75.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,75.0,0.0,0.0,75.0,75.0,75.0,75.0};
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
	            intensityVoxel2BeamsIntensidad = 0.0;
	            for(int b=0;b<beams.size(); b++){
	            	beam = beams.get(b);
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	intensityVoxel2Beams+= solutionVector[beam] * radiation;
	            	intensityVoxel2BeamsIntensidad+= vectorAperturas[beam] * radiation;
	         
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
            		pen += w[i]* 2* (dd[i] - intensityVoxel2Beams)*-(intensityVoxel2BeamsIntensidad);
            	}else{
   					pen += w[i] *2 * Math.max((intensityVoxel2Beams - dd[i]),0) ;
   				}
            	
	        }
			//System.out.println(i+" "+Zmax[i]+" ");
			score=pen/aux_index.size();
			solution.scores[i]=score;
			F+=pen/aux_index.size();
			//System.out.println("organo: "+solution.scores[i]);
			//System.out.println("organo: "+w[i]);
		}
		solution.singleObjectiveValue=F;
		//System.out.println("Solucion: "+F);
		return F;
	}
	
	public  void printBAC() {
		for(int i =0;i<initialBACs.length;i++) {
			System.out.print(initialBACs[i]+" ");
		}
		System.out.println();
		
	}
	public static void printAperture(TreatmentPlan sol) {
		
		for(int i=0;i<sol.apertures.size();i++) {

			for(int x=0;x<sol.apertures.get(i).size();x++) {
				
				sol.apertures.get(i).get(x).printAperture();		
			}
			System.out.println();			
		}
		System.out.println();
	}
	public static void printFirstSolution(TreatmentPlan sol){
		Vector<List<Aperture>> stations = sol.apertures;
		
		//recorrido de angulos
		for(int i=0; i<stations.size(); i++){
			System.out.println("Angulo: "+i);
			List<Aperture> list_apertures = stations.get(i);
			//recorrido de aperturas por angulo
			for(int j=0; j<list_apertures.size(); j++){
				int intensity = (int) list_apertures.get(j).intensity;
				int[][] matrix = list_apertures.get(j).aperture;
				System.out.println("Apertura N�: "+j);
				System.out.println("Intensidad: "+intensity+"\n");
				
				for(int x=0;x<matrix.length; x++){
					System.out.println(matrix[x][0]+", "+matrix[x][1]);
				}
				
			}
			System.out.println("\n");
		}
		
	}
	
}
