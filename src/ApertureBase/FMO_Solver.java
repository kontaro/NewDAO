package ApertureBase;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import IMRT_DAO.DDM;
import IMRT_DAO.Organs;
import IMRT_DAO.TreatmentPlan;
import gurobi.*;


/**
*
* @author Moyano
*/

public class FMO_Solver {
	
    public int organs;//number of organ
    public int beams;//number of beams
    public int []R;//number of voxels by region
    public int aperture;//number of aperture
    
    public int[] bmlts;// #number of beamlets by angles
    public double[] weight; //weights of objective
    public double[] EUD0;
    
    
    public int totalBmlts; //total beamlets
    public int[] angles;
    
    public double[] newIntensity; 
    
    public double[] LB;
    public double[] UB;
    public Boolean[] isTarget;
    public double epsilon;
    //public double t;
    public double[] x;
    public String jobThreadID;
    public String solver;
    public double maxIntensity;
    public double minIntensity;
    public double objVal;
    
    public double[] eud;
	GRBEnv env;
	GRBModel model;
	DDM M;
	TreatmentPlan sol;
	public FMO_Solver(int beams,DDM M,Organs[] o,int[] selAngles,double[] dd,double[] weight) throws GRBException {
		this.beams=beams;
		this.eud=dd;
		R=new int[o.length];
		organs=o.length;
		bmlts=new int[beams];
		this.weight=weight;
		bmlts=selAngles;


		this.M=M;

		for(int i=0;i<R.length;i++) {
			R[i]=M.OrganVoxels[i];
		}
		setEnv();
		//setModelFMO();
		setModelFMO_OARsPen();
		//writeModel();
		model.dispose();
        env.dispose();
		
	}
	
	public void setNewModel(TreatmentPlan sol) {
		aperture=sol.apertures.get(0).size();
		this.sol=sol;
		boolean error=true;
		do {
		
			try {
				setEnv();
				setModelFMO();
				//setModel();
				//model.dispose();
		        //env.dispose();
		        error=false;
			} catch (GRBException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}while(error);
		
		
		
	}
	
	public void reset() {
		try {
			model.reset();
		} catch (GRBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public void setEnv() throws GRBException {
	      this.env = new GRBEnv(true);
	      env.set("OutputFlag", "0");
	     // env.set("logFile", "mip1.log");
	      
	      //env.set(GRB.IntParam.Threads, 4);
	      env.set("Aggregate","0");
	      env.set("Method", "2");
	      env.set("NormAdjust", "2");
	      env.start();
	}
	/**
	 * Setea el modelo que minimiza la penalizacion de todos los organos,
	 * no tiene restricciones
	 * @throws GRBException
	 */
	
	
	
	
	/**
	 * utiliza el modelo minimiza la penalizacion de los organos en riesgo y tiene como restricion el eud 0 para el tumor
	 * @throws GRBException
	 */
	public void setModelFMO() throws GRBException {
		this.model = new GRBModel(env);
		model.set(GRB.StringAttr.ModelName, "Direct Aperture Optimization");
		
		
		//set Variables
		
		// set variables for intensity an V
		GRBVar[] bmlt= new GRBVar[beams];//intensity for beam,for aperture
		GRBVar[][]dOrgan= new GRBVar[(R.length)][];
		GRBVar[][] voxel=new GRBVar[(R.length)][];
		int indexi=0,indexj=0;
		
		for (int i = 0; i < this.beams; i++) {

             bmlt[i] = model.addVar(0,20,0.0, GRB.CONTINUOUS ,"bmlt" +"["+(i+1) +"]");
            
        }
		
		for (int i = 0; i < (this.R.length); i++) {
			dOrgan[i]=new GRBVar[R[i]];
            for (int j = 0; j < R[i]; j++) {
            	indexi=i+1;
            	indexj=j+1;
            	dOrgan[i][j] = model.addVar((GRB.INFINITY*-1), GRB.INFINITY,0.0, GRB.CONTINUOUS ,"dOrgan" + indexi + "[" + indexj+"]");
                
            }
        }
		for (int i = 0; i < (this.R.length); i++) {
			voxel[i]=new GRBVar[R[i]];
            for (int j = 0; j < R[i]; j++) {
            	indexi=i+1;
            	indexj=j+1;
            	voxel[i][j] = model.addVar((GRB.INFINITY*-1), GRB.INFINITY,0.0, GRB.CONTINUOUS ,"v" + indexi + "[" + indexj+"]");
                
            }
        }

		//set constraints
		//for 
	
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Integer key, beamblet, totalBeamblets,beamIndex, count_voxel;
		Double radiation, coefficent;
		int diffBeamblets=0;
		
		for (int o = 0;  o< organs; o++) {
			aux_index = index_dao_ddm.get(o);
			aux_values = value_dao_ddm.get(o);
			keys = aux_index.keys();
			
			//Recorremos claves de voxel por organo para su evaluación
			
			count_voxel = 0;
			while(keys.hasMoreElements()){
				GRBLinExpr voxelRadiation= new GRBLinExpr();
				
				GRBLinExpr voxelRadiationOAR= new GRBLinExpr();
				if(o!=0) {
					voxelRadiationOAR.clear();
				}
				key = keys.nextElement();
	            beams = aux_index.get(key);
	            //repeat(beams);
	            for(int b=0;b<beams.size(); b++){ // de aqui vamos a sacar el beam (indice del angulo)
	            	value_index_key = key+"-"+beams.get(b);
	            
	            	radiation = aux_values.get(value_index_key);
	            	beamblet = beams.get(b);
	            	totalBeamblets = 0;
	            	beamIndex = 0;
	            	diffBeamblets=0;
	            	for(int z=0; z<bmlts.length;z++) {
	            		totalBeamblets+= bmlts[z];
	            		if(beamblet < totalBeamblets){
	            			beamIndex = z;
	            			break;
	            		}
	            		diffBeamblets+=bmlts[z];
	            	}
	            	
            		coefficent = radiation;
            		//coefficent = (double)aper[beamIndex][a][(beamblet-diffBeamblets)]*radiation;
            		if(o==0) {
            			
            			coefficent=coefficent*-1;
						voxelRadiation.addTerm(coefficent,bmlt[beamblet] );
						//System.out.println(coefficent);
					}
					else {
						voxelRadiation.addTerm(coefficent,bmlt[beamblet] );
					}
	            		
	            	
	            	
	            }
	            GRBLinExpr dm=new GRBLinExpr();
				dm.addTerm(1,dOrgan[o][count_voxel]);
			
				
				
	            if(o==0) {
	          
	            	double constEud=eud[o];
					voxelRadiation.addConstant(constEud);
					model.addConstr(dm,GRB.EQUAL,voxelRadiation,"dOrgan"+o+"["+(count_voxel+1)+"]");
					
				}
				else{
					
					double constEud=eud[o]*-1;
					voxelRadiation.addConstant(constEud);
					model.addConstr(dm,GRB.EQUAL,voxelRadiation,"dOrgan"+o+"["+(count_voxel+1)+"]");
				}
				
	            GRBLinExpr V=new GRBLinExpr();
				V.addTerm(1,voxel[o][count_voxel]);
				
				if(o==0) {
					
					model.addConstr(V, GRB.EQUAL,dm, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
					
					
				}else {

					model.addConstr(V , GRB.GREATER_EQUAL, dm, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
				}
				
				count_voxel++;
			}
			
			//System.out.println("");
		}
		
		//set model
		GRBQuadExpr objFunc= new GRBQuadExpr();
		for(int o = 0;  o< organs; o++) {
			double coef=(double)((weight[o]/R[o]));
			
			for (int j = 0; j < R[o]; ++j) {

				objFunc.addTerm(coef, voxel[o][j],voxel[o][j]);
				
			}
		}

		model.setObjective(objFunc, GRB.MINIMIZE);
		model.update();
        //writeModel();
		model.optimize();
        model.update();
        //model.computeIIS();
        //model.write("mod.ilp");
		double[]getIntensity=new double[this.beams];
        for (int i = 0; i < this.beams; ++i) {

        	getIntensity[i]=bmlt[i].get(GRB.DoubleAttr.X);
        	//String varName="intensity"+i+"."+j;
        	//getIntensity[i][j]=model.getVarByName(varName).get(GRB.DoubleAttr.X);
        	//intensity[i][j] = model.addVar(minIntensity,maxIntensity,0.0, GRB.CONTINUOUS ,"Intensity" + i + "." + j);
            
            
        }
		//GRBVar[] vars = model.getVars();
        newIntensity=getIntensity;
        objVal=model.get(GRB.DoubleAttr.ObjVal);
        
		
		
	}

	
	/**
	 * utiliza el modelo minimiza la penalizacion de los organos en riesgo y tiene como restricion el eud 0 para el tumor
	 * @throws GRBException
	 */
	public void setModelFMO_OARsPen() throws GRBException {
		this.model = new GRBModel(env);
		model.set(GRB.StringAttr.ModelName, "Direct Aperture Optimization");
		
		
		//set Variables
		
		// set variables for intensity an V
		GRBVar[] bmlt= new GRBVar[beams];//intensity for beam,for aperture
		GRBVar[][]dOrgan= new GRBVar[(R.length)][];
		GRBVar[][] voxel=new GRBVar[(R.length)][];
		int indexi=0,indexj=0;
		
		for (int i = 0; i < this.beams; i++) {

             bmlt[i] = model.addVar(0,20,0.0, GRB.CONTINUOUS ,"bmlt" +"["+(i+1) +"]");
            
        }
		
		for (int i = 0; i < (this.R.length); i++) {
			dOrgan[i]=new GRBVar[R[i]];
            for (int j = 0; j < R[i]; j++) {
            	indexi=i+1;
            	indexj=j+1;
            	dOrgan[i][j] = model.addVar((GRB.INFINITY*-1), GRB.INFINITY,0.0, GRB.CONTINUOUS ,"dOrgan" + indexi + "[" + indexj+"]");
                
            }
        }
		for (int i = 0; i < (this.R.length); i++) {
			voxel[i]=new GRBVar[R[i]];
            for (int j = 0; j < R[i]; j++) {
            	indexi=i+1;
            	indexj=j+1;
            	voxel[i][j] = model.addVar((GRB.INFINITY*-1), GRB.INFINITY,0.0, GRB.CONTINUOUS ,"v" + indexi + "[" + indexj+"]");
                
            }
        }

		//set constraints
		//for 
	
		ArrayList<Hashtable<Integer, ArrayList<Integer>>> index_dao_ddm = M.index_dao_ddm;
		ArrayList<Hashtable<String, Double>> value_dao_ddm = M.value_dao_ddm;
		Hashtable<Integer, ArrayList<Integer>> aux_index;
		Hashtable<String, Double> aux_values;
		Enumeration<Integer> keys;
		ArrayList<Integer> beams;
		String value_index_key;
		Integer key, beamblet, totalBeamblets,beamIndex, count_voxel;
		Double radiation, coefficent;
		int diffBeamblets=0;
		
		for (int o = 0;  o< organs; o++) {
			aux_index = index_dao_ddm.get(o);
			aux_values = value_dao_ddm.get(o);
			keys = aux_index.keys();
			
			//Recorremos claves de voxel por organo para su evaluación
			
			count_voxel = 0;
			while(keys.hasMoreElements()){
				GRBLinExpr voxelRadiation= new GRBLinExpr();
				
				GRBLinExpr voxelRadiationOAR= new GRBLinExpr();
				if(o!=0) {
					voxelRadiationOAR.clear();
				}
				key = keys.nextElement();
	            beams = aux_index.get(key);
	            //repeat(beams);
	            for(int b=0;b<beams.size(); b++){ // de aqui vamos a sacar el beam (indice del angulo)
	            	value_index_key = key+"-"+beams.get(b);
	            
	            	radiation = aux_values.get(value_index_key);
	            	beamblet = beams.get(b);
	            	totalBeamblets = 0;
	            	beamIndex = 0;
	            	diffBeamblets=0;
	            	for(int z=0; z<bmlts.length;z++) {
	            		totalBeamblets+= bmlts[z];
	            		if(beamblet < totalBeamblets){
	            			beamIndex = z;
	            			break;
	            		}
	            		diffBeamblets+=bmlts[z];
	            	}
	            	
            		coefficent = radiation;
            		//coefficent = (double)aper[beamIndex][a][(beamblet-diffBeamblets)]*radiation;
            		if(o==0) {
            			
            			coefficent=coefficent;
						voxelRadiation.addTerm(coefficent,bmlt[beamblet] );
						//System.out.println(coefficent);
					}
					else {
						voxelRadiation.addTerm(coefficent,bmlt[beamblet] );
					}
	            		
	            	
	            	
	            }
	            GRBLinExpr dm=new GRBLinExpr();
				
			
				
				
	            if(o==0) {
	            	dm.addTerm(1,dOrgan[o][count_voxel]);
	            	double constEud=eud[o];
	            	//voxelRadiation.addConstant(constEud);
	            	model.addConstr(dm,GRB.EQUAL,voxelRadiation,"dOrgan"+o+"["+(count_voxel+1)+"]");
					model.addConstr(dm,GRB.GREATER_EQUAL,constEud,"dOrganP"+o+"["+(count_voxel+1)+"]");
				}
				else{
					dm.addTerm(1,dOrgan[o][count_voxel]);
					double constEud=eud[o]*-1;
					voxelRadiation.addConstant(constEud);
					model.addConstr(dm,GRB.EQUAL,voxelRadiation,"dOrgan"+o+"["+(count_voxel+1)+"]");
				}
				
	            GRBLinExpr V=new GRBLinExpr();
				V.addTerm(1,voxel[o][count_voxel]);
				
				if(o==0) {
				//	model.addConstr(V , GRB.GREATER_EQUAL, dm, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
					
				}else {

					model.addConstr(V , GRB.GREATER_EQUAL, dm, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
				}
				
				count_voxel++;
			}
			
			//System.out.println("");
		}
		
		//set model
		GRBQuadExpr objFunc= new GRBQuadExpr();
		for(int o = 1;  o< organs; o++) {
			double coef=(double)((weight[o]/R[o]));
			
			for (int j = 0; j < R[o]; ++j) {

				objFunc.addTerm(coef, voxel[o][j],voxel[o][j]);
				
			}
		}

		model.setObjective(objFunc, GRB.MINIMIZE);
		model.update();
        //writeModel();
		model.optimize();
        model.update();
        //model.computeIIS();
        model.write("mod4.lp");
		double[]getIntensity=new double[this.beams];
        for (int i = 0; i < this.beams; ++i) {

        	getIntensity[i]=bmlt[i].get(GRB.DoubleAttr.X);
        	//String varName="intensity"+i+"."+j;
        	//getIntensity[i][j]=model.getVarByName(varName).get(GRB.DoubleAttr.X);
        	//intensity[i][j] = model.addVar(minIntensity,maxIntensity,0.0, GRB.CONTINUOUS ,"Intensity" + i + "." + j);
            
            
        }
		//GRBVar[] vars = model.getVars();
        newIntensity=getIntensity;
        objVal=model.get(GRB.DoubleAttr.ObjVal);
        System.out.print("LLegue");
        
		
		
	}

	
	
	
	public void writeModel() throws GRBException {
		model.write("out1.lp");
		model.write("out1.mst");
		model.write("out1.mps");
		model.write("out1.sol");
		
	}
	
	
	public void cuadraticSum(TreatmentPlan sol) {
		
	}
	
	public void returnIntensity() {
		
	}
	public void repeat(ArrayList<Integer> lista) {
		for(int i=0;i<lista.size();i++) {
			for(int j=0;j<lista.size();j++) {
				if(lista.get(i)==lista.get(j) && i!=j) {
					System.out.println("se repite");
				}
				
			}
		}
	}
}
