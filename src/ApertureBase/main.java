package ApertureBase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import IMRT_DAO.*;
import gurobi.GRBException;



public class main {
    public static int[][] beamletSet;
    public static int numAngles;  
    public static int[] initialBACs;
    public static int numOrgans;
    public static Organs[] o;
    public static String pathFile = "";
    public static int[][] Vx; 				       //Vx index. % of the volume receiving more than x% of the prescribed dose
	  //One row per organ (usually only for target regions)
	public static int[][] Dx;
	public static int init_intensity;   
	public static Vector<int[][]> apertureShapes;
    public static int option;
    public static int selectionCriterion;
    public static int[] beamletAngles;
    public static int max_intensity;                // # DAO : Maximum intensity for apertures
    public static int [] zmax,zmin;
    public static double [] dd;
    public static double [] fmo;
    public static double [] fmround;
    public static double [] fmround2;
    public static double [] fmround4;
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		//Crear DDM y cargar amgulos y pesos
		String inputFile = args[0];
		readInputFile(inputFile);
		beamletSet = new int[numAngles][4];
		beamletAngles = new int[numAngles];
		
		for(int i = 0;i<1;i++) {
			printBAC() ;
	    
		    //Constructor DAO for DDM
		    //long time=System.currentTimeMillis();
		    String beamInfoDir = pathFile+"beamsInfo.txt";
		    
		    beamletSet = loadNumBixels(beamInfoDir);
		    setBeamlet();
		    int totalbmlt = totalBeamlet();
		    DDM M; 
		    
		    ///int totalbmlt =206;
			
			//Organs[] or=new Organs[1];
			/*
			int[] initialBAC=new int[1];
			int[] beamletAngle=new int[1];
			initialBAC[0]=initialBACs[0];
	
			beamletAngle[0]=beamletAngles[0];
	
			or[0]=o[0];
			*/
			//M = new DDM(or, initialBAC, pathFile,beamletAngle,totalbmlt);
				
		    M = new DDM(o, initialBACs, pathFile,beamletAngles,totalbmlt);
			//Resolver FMO por cada Angulo
		    
			System.out.println("");
			//Se genera 4 FMO por cada conjunto de angulo (Original - Redondeado - MOD2 - MOD4)
			//Se calcula los valores para los nuevos FMO
			double[]wT=new double[]{0,1,1};
			double[]wOAR=new double[]{0,1,1};
		    setDD();
		    //M.writeDDM(totalbmlt, beams, org, fileName, processName)
		    //TreatmentPlan solution = null;
			//evaluateSolution(solution,M,o,w);
		    //Solver
		    try {
				FMO_Solver newModel=new FMO_Solver(totalbmlt,M,o,beamletAngles,dd,wT);
				//Gurobi_Solver newModel=new Gurobi_Solver(totalbmlt,M,or,beamletAngle,dd,w);
				System.out.print("OBJ: "+newModel.objVal+" ");
				fmo=newModel.newIntensity;
			} catch (GRBException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		    evaluateSolution(M,o,wT,fmo);
		    System.out.println();
		    //printBMLTS();
		    int total=0;
		    int bot=0;
		    double [] x1;
		    double [] x2;
		    double [] x3;
		    TreatmentPlan solution = new TreatmentPlan(init_intensity, max_intensity, 1, 1, 1, 1, 1, numAngles,numOrgans);
		    setBeams(solution, pathFile);
		    MLC obt=new MLC();
		    x1=roundBMLTS();
		    setInitialApertureShape( solution, x1);
		    for(int j=0;j<apertureShapes.size();j++) {
		    	obt.getAperture(apertureShapes.get(j));
		    	total=obt.matrix+total;
		    	bot=bot+obt.beam_on_time;
		    }
		    
		    System.out.print(" "+total);
		    System.out.println(" "+bot);
	
		    total=0;
		    bot=0;
		    evaluateSolution(M,o,wT,x1);
		    x2=roundBMLTS(2);
		    setInitialApertureShape( solution, x2);
		    for(int j=0;j<apertureShapes.size();j++) {
		    	obt.getAperture(apertureShapes.get(j));
		    	total=obt.matrix+total;
		    	bot=bot+obt.beam_on_time;
		    }
		    System.out.println(total);
		    System.out.println(bot);
		    total=0;
		    bot=0;
		   evaluateSolution(M,o,wT,x2);
		    x3=roundBMLTS(4);
		    setInitialApertureShape( solution, x3);
		    for(int j=0;j<apertureShapes.size();j++) {
		    	obt.getAperture(apertureShapes.get(j));
		    	bot=bot+obt.beam_on_time;
		    	total=obt.matrix+total;
		    }
		    System.out.println(total);
		    System.out.println(bot);
		    total=0;
		    evaluateSolution(M,o,wT,x3);
	
		    //System.out.println("param bmlt :=	"+totalBeamlet()+"; param R1 :="+M.OrganVoxels[0]+"; param R2 :="+M.OrganVoxels[1]+"; param R3 :="+M.OrganVoxels[2]+";");
		    for(int h =0;h<initialBACs.length;h++) {
				initialBACs[h]=initialBACs[h]+5;
			}
		}
	}
	
	
	
	//******* FUNCIONES    ******
	public static void printBMLTS() {
		for(int i =0;i<fmo.length;i++) {
			System.out.println(fmo[i]);
		}
		
	}
	public static void printBAC() {
		for(int i =0;i<initialBACs.length;i++) {
			System.out.print(initialBACs[i]+" ");
		}
		System.out.println();
		
	}
	
	public static double [] roundBMLTS() {
		double [] x;
		x=new double[fmo.length];
		for(int i =0;i<fmo.length;i++) {
			x[i]=Math.round(fmo[i]);
		}
		return x;
		
	}
	public static double [] roundBMLTS(int mod) {
		double [] x;
		x=new double[fmo.length];
		for(int i =0;i<fmo.length;i++) {
			double aux=fmo[i]/mod;
			aux=Math.round(aux);
			x[i]=aux*mod;
		}
		return x;
	}
	
	public static void setDD() {
		dd = new double[o.length];
		for(int i=0;i<o.length;i++) {
			if(o[i].isTarget) {
				dd[i]=o[i].doseLB;
			}else {
				dd[i]=o[i].doseUB;
			}
		}
	}
	public static int totalBeamlet() {
		int sum=0;
		for(int i=0 ;i<numAngles;i++) {
			sum+=beamletAngles[i];
		}
		return sum;
	}
	public static void setBeamlet() {
		for(int i=0 ;i<numAngles;i++) {
			beamletAngles[i]=beamletSet[initialBACs[i]][1];
		}
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
	
	public static void readInputFile(String dir) throws IOException{
        
    	String sp="\\s+";
        //String dir = "./inputFile.txt";
        File f = new File(dir);
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
                int max_iter = Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get max time LS
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                int max_time = Integer.parseInt(auxReader[0]);
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
        fileIn.close();
	 } 
	
	public static double evaluateSolution( DDM M, Organs[] o, double w[],double[] solutionVector){
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Double radiation, intensityVoxel2Beams;
		Integer key, beam;
		
		//double[] solutionVector = solution.intensity; //solucion a evaluar en forma de vector (agregar como variable del Treatment Plan y crer metodo que la cree)
		//double[] solutionVector = new double[]{4,4,4,21,21,15,10,4,4,21,21,21,21,15,15,10,10,4,21,21,15,15,15,15,10,8,4,21,15,15,15,15,10,10,10,4,21,15,15,10,10,10,10,10,4,21,15,15,15,15,10,10,10,4,15,21,21,15,10,10,10,4,10,21,10,8,4,4,4,4,10,10,10,11,4,4,4,10,10,10,17,17,14,11,4,10,10,10,17,14,14,14,4,11,11,11,11,11,17,10,4,4,10,10,14,14,14,11,10,4,4,14,14,14,14,14,14,4,4,10,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,9,4,4,4,4,4,4,9,9,9,9,9,4,4,4,9,13,13,13,14,6,4,4,4,4,4,6,6,9,9,4,4,4,4,4,4,4,9,9,4,4,4,4,4,4,4,14,4,4,4,4,9,9,9,9,9,4,4,6,9,6,4,4,4,4,4,4,4,4,7,5,4,4,5,9,9,9,4,4,7,7,13,13,9,5,4,4,5,13,13,13,13,9,9,5,4,9,9,9,9,13,9,9,9,4,7,9,9,9,9,9,13,5,4,7,7,7,7,7,7,7,9,4,5,7,7,7,7,7,9,13,4,4,4,5,7,4,4,4,4,4,4,4,4,4,4,4,4,7,4,4,4,12,11,11,11,11,11,7,4,4,7,21,12,12,12,12,12,4,4,4,21,12,12,12,7,7,7,4,4,21,12,12,12,12,12,4,4,11,11,11,11,11,12,12,4,7,11,11,11,12};
		//double[] solutionVector = new double[]{10,20,20,14,20,19,11,10,20,5,6,7,16,19,18,14,9,20,10,4,5,14,17,16,11,10,20,12,4,5,13,17,15,9,10,20,13,5,4,13,18,12,8,11,19,14,7,4,14,17,9,10,12,20,9,4,17,18,8,14,10,20,20,20,5,14,13,20,20,19,11,18,20,5,10,7,10,15,20,15,3,15,6,3,12,17,20,9,2,16,4,3,9,20,20,6,3,19,2,4,1,5,17,17,0,0,17,0,4,0,0,17,20,1,1,19,1,0,1,0,13,0,1,19,0,4,0,7,11,4,0,20,15,3,0,0,0,9,10,6,20,0,1,3,8,10,19,7,20,6,3,6,8,9,20,5,9,20,8,1,4,5,6,17,6,7,13,8,1,2,3,2,10,5,8,16,10,2,3,3,1,5,4,4,20,9,7,9,3,4,3,5,5,3,2,4,7,4,6,2,7,4,0,0,20,9,20,20,20,20,12,1,20,20,20,20,19,18,9,10,17,20,18,15,14,18,17,7,13,16,19,14,15,14,20,20,6,17,15,20,14,16,12,16,19,5,18,19,20,20,19,13,14,19,3,19,20,20,20,20,12,7,20,4,17,19,20,6,0,17,1,7,18,2,0,0,3,4,0,0,0,0,0,0,0,2,0,0,11,17,19,10,11,15,0,10,10,8,17,18,10,13,14,11,2,8,14,20,20,13,13,10,12,0,4,14,19,14,11,8,8,12,0,17,19,9,8,7,7,11,20,19,4,5,6,6};
		//double[] solutionVector=new double[] {10.31,20.00,20.00,14.49,20.00,19.37,11.44,10.08,20.00,5.23,6.49,7.40,16.28,18.54,17.56,13.60,9.32,20.00,10.41,3.67,4.82,13.81,16.93,15.60,10.51,9.81,20.00,11.63,4.02,4.74,12.96,16.99,14.88,8.83,9.79,20.00,13.07,4.90,3.71,13.46,17.62,12.28,8.30,11.12,19.34,13.87,7.48,4.13,13.53,16.78,9.34,9.68,12.06,20.00,8.51,4.18,17.12,18.09,7.79,14.35,10.17,20.00,20.00,20.00,4.55,14.37,13.27,20.00,20.00,19.11,10.89,18.46,20.00,4.95,10.04,6.66,9.82,15.16,19.87,14.52,3.23,14.87,6.21,3.34,11.99,16.70,20.00,8.77,2.48,16.18,3.84,2.63,9.30,20.00,20.00,5.58,2.86,19.32,2.23,3.68,0.68,5.38,17.18,17.28,0.00,0.00,16.50,0.00,4.28,0.00,0.00,17.06,20.00,1.04,0.58,19.19,0.92,0.00,0.54,0.45,13.30,0.00,1.39,18.92,0.00,4.37,0.00,7.01,11.15,4.37,0.00,20.00,15.45,2.62,0.00,0.00,0.00,8.54,9.78,5.60,20.00,0.00,1.43,3.35,8.33,10.28,19.01,6.61,20.00,5.81,2.98,6.48,8.28,8.88,20.00,5.25,8.81,20.00,8.45,1.12,3.79,5.25,5.60,16.52,6.01,7.16,12.88,7.73,0.68,1.61,3.07,1.61,10.29,5.49,7.75,16.13,9.73,1.64,3.13,2.94,1.12,5.49,3.85,3.76,20.00,8.51,6.98,9.28,2.81,3.76,3.33,5.10,5.30,2.95,2.44,3.95,7.35,3.54,5.58,1.85,6.66,3.53,0.00,0.00,20.00,9.14,19.85,20.00,20.00,20.00,12.22,1.47,20.00,20.00,20.00,20.00,18.61,17.75,8.77,9.74,16.53,19.93,18.08,15.31,13.80,18.36,16.93,6.52,13.03,16.03,19.12,13.84,14.91,13.61,20.00,19.84,6.39,16.74,15.23,19.66,13.69,15.69,12.40,15.65,18.58,4.84,18.19,18.92,20.00,20.00,18.85,12.91,13.85,19.30,3.44,19.49,20.00,20.00,20.00,20.00,12.08,7.44,19.75,3.63,17.48,18.65,20.00,6.44,0.00,17.11,1.27,7.16,18.34,1.94,0.00,0.00,2.97,3.64,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.83,0.00,0.00,10.85,17.47,18.61,9.50,11.06,14.53,0.00,9.56,9.64,8.43,17.10,17.73,10.13,12.52,13.76,10.97,1.70,7.53,13.55,20.00,20.00,13.21,12.79,9.77,12.10,0.00,3.59,14.10,19.19,13.58,10.81,8.01,8.19,11.79,0.00,17.35,19.34,9.14,8.15,6.98,7.11,11.09,20.00,18.97,3.56,5.08,5.88,5.87};
		//solutionVector=new double[] {12.40409206207636, 12.40409206207636, 12.40409206207636, 12.4040920620766, 13.404092062076922, 12.40409206207636, 12.40409206207636, 12.40409206207636, 12.40409206207636, 12.4040920620766, 12.4040920620766, 12.4040920620766, 12.4040920620766, 12.4040920620766, 11.404092062076039, 11.404092062076039, 11.404092062076039, 11.404092062076039, 12.4040920620766, 12.4040920620766, 12.4040920620766, 12.4040920620766, 13.404092062076922, 12.40409206207636, 12.40409206207636, 12.40409206207636, 12.40409206207636, 12.4040920620766, 12.4040920620766, 12.4040920620766, 12.4040920620766, 13.404092062076922, 12.40409206207636, 12.40409206207636, 12.40409206207636, 12.40409206207636, 12.835024060238142, 12.835024060238142, 12.835024060238142, 12.835024060238142, 13.835024060238464, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060238142, 11.835024060237581, 11.835024060237581, 11.835024060237581, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060237902, 11.835024060237581, 11.835024060237581, 11.835024060237581, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060238142, 13.835024060238464, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060237902, 12.835024060237902, 6.201184891891641, 6.201184891891641, 6.201184891891641, 7.201184891891915, 7.201184891891915, 7.201184891891915, 13.304122113484263, 6.201184891891641, 6.201184891891641, 6.201184891891641, 7.201184891891915, 7.201184891891915, 7.201184891891915, 7.201184891891915, 13.304122113484263, 13.304122113484263, 13.304122113484263, 13.304122113484263, 14.304122113484537, 7.201184891891915, 7.201184891891915, 7.201184891891915, 13.304122113484263, 13.304122113484263, 13.304122113484263, 13.304122113484263, 14.304122113484537, 7.201184891891915, 7.201184891891915, 7.201184891891915, 7.201184891891915, 13.479585680813999, 13.479585680813999, 13.479585680813999, 13.479585680813999, 14.479585680814273, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 13.479585680813999, 13.479585680813999, 13.479585680813999, 13.479585680813999, 14.479585680814273, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 13.479585680813999, 13.479585680813999, 13.479585680813999, 14.479585680814273, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 7.37664845922165, 10.944008890957328, 10.944008890957328, 11.944008890957623, 10.944008890953835, 10.944008890953835, 10.944008890957328, 10.944008890957328, 10.944008890957328, 10.944008890957328, 11.944008890957623, 10.944008890953835, 10.944008890953835, 10.944008890953835, 10.944008890957328, 10.944008890957328, 10.944008890957328, 10.944008890957328, 11.944008890957623, 10.944008890953835, 10.944008890953835, 10.944008890953835, 10.944008890953835, 10.944008890957328, 10.944008890957328, 10.944008890957328, 10.944008890957328, 10.944008890957328, 9.944008890953537, 9.944008890953537, 9.944008890953537, 9.944008890953537, 7.515583649927876, 7.515583649927876, 7.515583649927876, 7.515583649927876, 8.515583649928175, 7.515583649924384, 7.515583649924384, 7.515583649924384, 7.515583649924384, 7.515583649927876, 7.515583649927876, 7.515583649927876, 7.515583649927876, 8.515583649928175, 7.515583649924384, 7.515583649924384, 7.515583649924384, 7.515583649924384, 7.515583649927876, 4.72427788792661, 4.72427788792661, 4.72427788792661, 4.72427788792661, 3.724277887922821, 3.724277887922821, 3.724277887922821, 3.724277887922821, 8.515583649928175, 7.515583649924384, 7.515583649924384, 7.515583649924384, 7.515583649924384, 7.515583649924384, 7.515583649924384, 7.515583649924384, 7.515583649924384, 9.379681954441573, 9.379681954441573, 9.379681954441573, 9.379681954441573, 9.379681954441573, 9.379681954441573, 9.379681954442336, 10.379681954442644, 9.379681954441573, 9.379681954441573, 9.379681954441573, 9.379681954441573, 9.379681954442336, 9.379681954442336, 9.379681954442336, 10.379681954442644, 9.379681954441573, 9.379681954441573, 9.379681954441573, 9.379681954441573, 9.379681954442336, 9.379681954442336, 9.379681954442336, 9.379681954442336, 10.379681954442644, 9.379681954441573, 9.379681954441573, 9.379681954441573, 9.379681954441573, 10.010592249165972, 10.010592249165972, 10.010592249165972, 10.010592249165972, 11.01059224916628, 10.010592249165208, 10.010592249165208, 10.010592249165208, 10.010592249165208, 10.010592249165972, 10.010592249165972, 10.010592249165972, 10.010592249165972, 11.01059224916628, 10.010592249165208, 10.010592249165208, 10.010592249165208, 10.010592249165208, 10.010592249165972, 10.010592249165972, 10.010592249165972, 10.010592249165972, 11.01059224916628, 10.010592249165208, 10.010592249165208, 10.010592249165208, 10.010592249165208, 10.010592249165972, 10.010592249165972, 10.010592249165972, 10.010592249165972, 11.01059224916628, 10.010592249165208, 10.010592249165208, 10.010592249165208, 10.010592249165208, 10.010592249165972, 10.010592249165972, 11.01059224916628, 10.010592249165208, 10.010592249165208, 10.010592249165208, 17.387417269640824, 17.387417269640824, 17.387417269640824, 17.387417269640824, 17.387417269640824, 17.387417269640824, 18.65931412010837, 18.65931412010837, 18.65931412010837, 19.659314120108665, 17.387417269640824, 17.387417269640824, 17.387417269640824, 17.387417269640824, 18.65931412010837, 18.65931412010837, 18.65931412010837, 18.65931412010837, 18.65931412010837, 16.38741726964053, 16.38741726964053, 16.38741726964053, 16.38741726964053, 19.84135811131747, 8.120820075129703, 8.120820075129703, 8.120820075129703, 9.12082007513, 6.84892322466216, 6.84892322466216, 6.84892322466216, 6.84892322466216, 19.84135811131747, 19.84135811131747, 19.84135811131747, 19.84135811131747, 20.841358111317767, 18.569461260849927, 18.569461260849927, 18.569461260849927, 18.569461260849927, 19.84135811131747, 8.120820075129703, 8.120820075129703, 8.120820075129703, 9.12082007513, 6.84892322466216, 6.84892322466216, 6.84892322466216, 19.84135811131747, 8.120820075129703, 8.120820075129703, 8.120820075129703, 9.12082007513, 6.84892322466216, 6.84892322466216, 6.84892322466216, 19.84135811131747, 19.84135811131747, 19.84135811131747, 20.841358111317767, 18.569461260849927, 18.569461260849927};
		
		
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
            		//pen += w[i] * Math.pow((dd[i] - intensityVoxel2Beams),2);
            	}else{
   					pen += w[i] * Math.pow(Math.max((intensityVoxel2Beams - dd[i]),0),2);
   				}
            	
	        }
			//System.out.println(i+" "+Zmax[i]+" ");
			score=pen/aux_index.size();
			//solution.scores[i]=score;
			F+=pen/aux_index.size();
			if(i == 0){
				
			}else {
				System.out.print("organo "+i+": "+score);
			}
			
		}
		//solution.singleObjectiveValue=F;
		System.out.print(" Solucion: "+F);
		return F;
	}

	public static void setInitialApertureShape(TreatmentPlan solution,double[] solutionVector){
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
	public static void setBeams(TreatmentPlan solution, String dir) throws IOException{
		int totalbmlt = totalBeamlet();;
		for(int i=0; i<numAngles; i++){
			solution.selAngles[i] = new Beam(beamletSet[initialBACs[i]][1],initialBACs[i], pathFile);
			totalbmlt=totalbmlt+solution.selAngles[i].beamlets;
		}
		solution.beamlets=totalbmlt;
	}

}
