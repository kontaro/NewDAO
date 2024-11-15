package MultiObj;

import java.io.IOException;
import java.util.ArrayList;

import IMRT_DAO.DDM;
import IMRT_DAO.TreatmentPlan;
import SingleObj.Algorithm;
import SingleObj.DAO;
import SingleObj.InitialSolution;
import gurobi.GRBException;

public class MODAO extends DAO {
	public MODAO(String inputFile) throws IOException {
		super();
		readInputFile(inputFile);
	    beamletSet = new int[numAngles][4];
		beamletAngles = new int[numAngles];
	    String beamInfoDir = pathFile+"beamsInfo.txt";
	    beamletSet = loadNumBixels(beamInfoDir);
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
			    
			    w = new double[]{5,1,1};
			    //-----Generar solucion Inicial------/
			    
			    InitialSolution newInitial=new InitialSolution(w,beamletSet,initialBACs,pathFile,M,o,dd,beamletAngles,init_intensity, 20, max_delta, max_iter, max_time, seed, 2, numAngles,numOrgans);
			    System.out.println("Run: "+ j);
			    ArrayList<TreatmentPlan> front= new ArrayList<TreatmentPlan>();
			    TreatmentPlan solution = newInitial.test4();
			    //evaluateSolution(solution,M,o,w);
			    MultiObjectiveAlgorithm heuristica= new MultiObjectiveAlgorithm(beamletSet,initialBACs,M,o,dd,beamletAngles,init_intensity, max_intensity, max_delta, max_iter, max_time, seed, step_intensity, numAngles,numOrgans);
				double initialvalue=solution.singleObjectiveValue;

				
				try {
					int flag=1;
					if(flag==0) {
					//heuristica.Test2(M, solution);
					//heuristica.WeightedSum(M, solution);
					//heuristica.initialFront(M, solution);
						front=heuristica.PLS2(M, solution);
					}
					else {
						//heuristica.ParetoSurfaceNavigation(M,serializeDataIn("tpPen-05-bac4.tpg"));
						front=heuristica.AttachedParetoSurfaceNavigationTest(M);
					}
				} catch (GRBException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}catch (ClassNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				//solution.evaluateSolution();
				//solution.scorePrint();
				
				System.out.println(solution.singleObjectiveValue+" "+solution.getLowIntensityperBeam()+" "+solution.getLowIntensityAngle()+" "+(System.currentTimeMillis()-time)+" "+solution.beamOnTime()+" "+heuristica.iterations);
				
				if(j==0) {
					for(int x=0;x<front.size();x++) {

						ArrayList<ArrayList<Double>> Doses=DoseVolume(front.get(x), front.get(x).M, o);
						String palabra="MOSurfaceIntersection_bac-"+x;
						doseprint(Doses,palabra,i,j);
							
					}	
				}
				//solution.numIntensity();
				;
				
				
	    	}
			for(int h =0;h<initialBACs.length;h++) {
				initialBACs[h]=initialBACs[h]+5;
				System.out.print(initialBACs[h]+" ");
			}
	    }
	    

	    
	    
	    
	}
}
