	package SingleObj;
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

public class Gurobi_Solver {
	
    public int organs;//number of organ
    public int beams;//number of organ
    public int []R;//number of voxels by region
    public int aperture;//number of aperture
    
    public int[] bmlts;// #number of beamlets by angles
    public double[] weight; //weights of objective
    public double[] EUD0;
    public double[] ar=new double[] {-10,8,2};
    
    
    public int totalBmlts; //total beamlets
    public int[] angles;
    
    public double[][] newIntensity; 
    
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
    
    public int[] eud;
	GRBEnv env;
	GRBModel model;
	DDM M;
	TreatmentPlan sol;
	public Gurobi_Solver(TreatmentPlan sol,DDM M,Organs[] o,int[] selAngles,double[] dd,double[] weight) throws GRBException {
		beams=sol.beams;
		this.eud=doubletoint(dd);
		R=new int[o.length];
		organs=o.length;
		bmlts=new int[beams];
		this.weight=weight;
		bmlts=selAngles;
		aperture=sol.apertures.get(0).size();
		minIntensity=0;
		maxIntensity=sol.max_intesity;
		this.M=M;
		this.sol=sol;
		for(int i=0;i<R.length;i++) {
			R[i]=M.OrganVoxels[i];
		}
		setEnv();
		setModelPen();
		//setModelORPen();
		writeModel();
		model.dispose();
        env.dispose();
		
	}
	public Gurobi_Solver(TreatmentPlan sol,DDM M,Organs[] o,int[] selAngles,double[] dd,double[] weight,int a,int n) throws GRBException {
		beams=sol.beams;
		this.eud=doubletoint(dd);
		R=new int[o.length];
		organs=o.length;
		bmlts=new int[beams];
		this.weight=weight;
		bmlts=selAngles;
		aperture=sol.apertures.get(0).size();
		minIntensity=n;
		maxIntensity=sol.max_intesity;
		this.M=M;
		this.sol=sol;
		for(int i=0;i<R.length;i++) {
			R[i]=M.OrganVoxels[i];
		}
		setEnv();
		//setModelgEUD();
		setModelORPen();
		//writeModel();
		model.dispose();
        env.dispose();
		
	}
	public Gurobi_Solver(TreatmentPlan sol,DDM M,Organs[] o,int[] selAngles,double[] dd,double[] weight,int min) throws GRBException {
		beams=sol.beams;
		this.eud=doubletoint(dd);
		R=new int[o.length];
		organs=o.length;
		bmlts=new int[beams];
		this.weight=weight;
		bmlts=selAngles;
		aperture=sol.apertures.get(0).size();
		minIntensity=min;
		maxIntensity=sol.max_intesity;
		this.M=M;
		this.sol=sol;
		for(int i=0;i<R.length;i++) {
			R[i]=M.OrganVoxels[i];
		}
		setEnv();
		setModelPen();
		//setModelORPen();
		//writeModel();
		model.dispose();
        env.dispose();
		
	}
	public Gurobi_Solver(TreatmentPlan sol,DDM M,Organs[] o,int[] selAngles,int[] dd,double[] weight) throws GRBException {
		beams=sol.beams;
		this.eud=dd;
		R=new int[o.length];
		organs=o.length;
		bmlts=new int[beams];
		this.weight=weight;
		bmlts=selAngles;
		aperture=sol.apertures.get(0).size();
		minIntensity=0;
		maxIntensity=sol.max_intesity;
		this.M=M;
		this.sol=sol;
		for(int i=0;i<R.length;i++) {
			R[i]=M.OrganVoxels[i];
		}
		setEnv();
		setModelPen();
		//setModelORPen();
		//writeModel();
		model.dispose();
        env.dispose();
		
	}
	public int[] doubletoint(double[]doubleArray) {
		 int[] intArray = new int[doubleArray.length];
		for (int i=0; i<intArray.length; ++i)
		    intArray[i] = (int) doubleArray[i];
		return intArray;
	}
	public void setNewModel(TreatmentPlan sol) {
		aperture=sol.apertures.get(0).size();
		this.sol=sol;
		boolean error=true;
		do {
		
			try {
				setEnv();
				setModelORPen();
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
	public void setModelPen() throws GRBException {
		this.model = new GRBModel(env);
		model.set(GRB.StringAttr.ModelName, "Direct Aperture Optimization");
		
		
		//set Variables
		
		// set variables for intensity an V
		GRBVar[][] intensity= new GRBVar[beams][aperture];//intensity for beam,for aperture
		GRBVar[][] voxel=new GRBVar[R.length][];
		int indexi=0,indexj=0;
		
		for (int i = 0; i < this.beams; ++i) {
            for (int j = 0; j < this.aperture; ++j) {
            	indexi=i+1;
            	indexj=j+1;
            	intensity[i][j] = model.addVar(minIntensity,20,0.0, GRB.CONTINUOUS ,"Intensity" +"["+ indexi + "." + indexj+"]");
            	intensity[i][j].set(GRB.DoubleAttr.Start, sol.apertures.get(i).get(j).intensity); 
            }
        }
		
		
		for (int i = 0; i < this.R.length; ++i) {
			voxel[i]=new GRBVar[R[i]];
            for (int j = 0; j < R[i]; ++j) {
            	indexi=i+1;
            	indexj=j+1;
            	if(i==0) {
            		voxel[i][j] = model.addVar(-GRB.INFINITY,GRB.INFINITY,0.0, GRB.CONTINUOUS ,"v" + indexi + "[" + indexj+"]");
            		
            	}
            	else{
            		voxel[i][j] = model.addVar(-GRB.INFINITY,GRB.INFINITY,0.0, GRB.CONTINUOUS ,"v" + indexi + "[" + indexj+"]");
            	}
            	
                
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
				key = keys.nextElement();
	            beams = aux_index.get(key);
	            
	            for(int b=0;b<beams.size(); b++){ // de aqui vamos a sacar el beam (indice del angulo)
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	beamblet = beams.get(b);
	            	totalBeamblets = 0;
	            	beamIndex = 0;
	            	diffBeamblets=0;
	            	for(int z=0; z<bmlts.length;z++) {
	            		totalBeamblets+= bmlts[z];//70 +69
	            		if(beamblet < totalBeamblets){//
	            			beamIndex = z;
	            			break;
	            		}
	            		diffBeamblets+=bmlts[z];
	            	}

	            	for(int a=0;a<aperture;a++) {
	            		coefficent = (double)sol.aperturesBmlts.get(beamIndex)[a][beamblet-diffBeamblets]*radiation;
	            		//coefficent = (double)aper[beamIndex][a][(beamblet-diffBeamblets)]*radiation;
	            		if(coefficent!=0 && o==0) {
	            			
							coefficent=coefficent*-1;
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][a] );
							//System.out.println(coefficent);
						}
						else {
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][a] );
						}
	            		
	            		
	            	}
	            }
	            if(o==0) {
					voxelRadiation.addConstant(eud[o]);
				}
				else{
					int constEud=eud[o]*-1;
					voxelRadiation.addConstant(constEud);
				}
				GRBLinExpr V=new GRBLinExpr();
				V.addTerm(1,voxel[o][count_voxel]);
				
				if(o==0) {
					model.addConstr(V, GRB.EQUAL, voxelRadiation, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
				}else {
					model.addConstr(V, GRB.GREATER_EQUAL, voxelRadiation, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
				}
				
				count_voxel++;
			}
		}
		
		//set model
		GRBQuadExpr objFunc= new GRBQuadExpr();
		for(int o = 0;  o< organs; o++) {
			double coef=(double)((weight[o]/R[o]));
			
			for (int j = 0; j < R[o]; ++j) {
				objFunc.addTerm(coef, voxel[o][j],voxel[o][j]);
				System.out.print("");
			}
		}

		model.setObjective(objFunc, GRB.MINIMIZE);
		model.update();
		writeModel();
        model.optimize();
        model.update();
        //model.computeIIS();
        //model.write("mod.ilp");
		double[][]getIntensity=new double[this.beams][this.aperture];
        for (int i = 0; i < this.beams; ++i) {
            for (int j = 0; j < this.aperture; ++j) {
            	getIntensity[i][j]=intensity[i][j].get(GRB.DoubleAttr.X);
            	//String varName="intensity"+i+"."+j;
            	//getIntensity[i][j]=model.getVarByName(varName).get(GRB.DoubleAttr.X);
            	//intensity[i][j] = model.addVar(minIntensity,maxIntensity,0.0, GRB.CONTINUOUS ,"Intensity" + i + "." + j);
                
            }
        }
		//GRBVar[] vars = model.getVars();
        newIntensity=getIntensity;
        objVal=model.get(GRB.DoubleAttr.ObjVal);
        
		
		
	}
	
	
	
	/**
	 * Setea el modelo que minimiza la penalizacion de todos los organos,
	 * no tiene restricciones
	 * @throws GRBException
	 */
	public void setModelgEUD() throws GRBException {
		this.model = new GRBModel(env);
		model.set(GRB.StringAttr.ModelName, "Direct Aperture Optimization");
		
		
		//set Variables
		
		// set variables for intensity an V
		GRBVar[][] intensity= new GRBVar[beams][aperture];//intensity for beam,for aperture
		GRBVar[][] voxel=new GRBVar[R.length][];
		GRBVar[][] a=new GRBVar[R.length][];
		int indexi=0,indexj=0;
		
		for (int i = 0; i < this.beams; ++i) {
            for (int j = 0; j < this.aperture; ++j) {
            	indexi=i+1;
            	indexj=j+1;
            	intensity[i][j] = model.addVar(0,20,0.0, GRB.CONTINUOUS ,"Intensity" +"["+ indexi + "." + indexj+"]");
            	intensity[i][j].set(GRB.DoubleAttr.Start, sol.apertures.get(i).get(j).intensity); 
            }
        }
		
		
		for (int i = 0; i < this.R.length; ++i) {
			voxel[i]=new GRBVar[R[i]];
            for (int j = 0; j < R[i]; ++j) {
            	indexi=i+1;
            	indexj=j+1;
            		voxel[i][j] = model.addVar(0,GRB.INFINITY,0.0, GRB.CONTINUOUS ,"v" + indexi + "[" + indexj+"]");

            }
        }
		for (int i = 0; i < this.R.length; ++i) {
			a[i]=new GRBVar[R[i]];
            for (int j = 0; j < R[i]; ++j) {
            	indexi=i+1;
            	indexj=j+1;

            		a[i][j] = model.addVar(0,GRB.INFINITY,0.0, GRB.CONTINUOUS ,"aux1" + indexi + "[" + indexj+"]");
 
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
				key = keys.nextElement();
	            beams = aux_index.get(key);
	            
	            for(int b=0;b<beams.size(); b++){ // de aqui vamos a sacar el beam (indice del angulo)
	            	value_index_key = key+"-"+beams.get(b);
	            	radiation = aux_values.get(value_index_key);
	            	beamblet = beams.get(b);
	            	totalBeamblets = 0;
	            	beamIndex = 0;
	            	diffBeamblets=0;
	            	for(int z=0; z<bmlts.length;z++) {
	            		totalBeamblets+= bmlts[z];//70 +69
	            		if(beamblet < totalBeamblets){//
	            			beamIndex = z;
	            			break;
	            		}
	            		diffBeamblets+=bmlts[z];
	            	}

	            	for(int j=0;j<aperture;j++) {
	            		coefficent = (double)sol.aperturesBmlts.get(beamIndex)[j][beamblet-diffBeamblets]*radiation;
	            		//coefficent = (double)aper[beamIndex][a][(beamblet-diffBeamblets)]*radiation;
	            		if(coefficent!=0 && o==0) {
	            			
							
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][j] );
							
						}
						else {
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][j] );
						}
	            		
	            		
	            	}
	            }
	            
				GRBLinExpr V=new GRBLinExpr();
				V.addTerm(1,voxel[o][count_voxel]);
				
				if(o==0) {
					model.addConstr(V, GRB.EQUAL, voxelRadiation, "voxelRadiation"+o+"["+(count_voxel+1)+"]");

				}else {
					model.addConstr(V, GRB.EQUAL, voxelRadiation, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
				}
				
				
				model.addGenConstrPow(voxel[o][count_voxel], a[o][count_voxel], ar[o], "dosisElevada"+o+"["+(count_voxel+1)+"]", null);
				count_voxel++;
			}
		}
		
		//set model
		
		
		GRBVar[] organ=new GRBVar[R.length];

		
		for (int i = 0; i < this.organs; ++i) {
            
        	organ[i] = model.addVar(0,GRB.INFINITY,0.0, GRB.CONTINUOUS ,"organ" + "[" + i+"]");
            
        }
		GRBVar[] b=new GRBVar[R.length];

		
		for (int i = 0; i < this.organs; ++i) {

        	b[i] = model.addVar(0,GRB.INFINITY,0.0, GRB.CONTINUOUS ,"b" + "[" + i+"]");

            
        }
		
		
		
		GRBQuadExpr objFunc= new GRBQuadExpr();
		
		
		
		
		for(int o = 0;  o< organs; o++) {
			double coef=(double)((R[o]));
			GRBLinExpr auxb=new GRBLinExpr();
			auxb.addTerm(coef,b[o]);
			GRBLinExpr first=new GRBLinExpr();
			
			for (int j = 0; j < R[o]; ++j) {
				first.addTerm(1, a[o][j]);
				System.out.print("");
			}
			model.addConstr(auxb, GRB.EQUAL, first, "AuxB"+"["+(o+1)+"]");
			model.addGenConstrPow(b[o], organ[o], ((double)(1/ar[o])), "organ", null);
			objFunc.addTerm(weight[o], organ[o]);
		}

		model.setObjective(objFunc, GRB.MINIMIZE);
		model.update();
		//a[1][3111].get(GRB.DoubleAttr.X);
		//writeModel();
        model.optimize();
        model.update();
        model.computeIIS();
        model.write("mod.ilp");
        model.write("out1.mps");

		System.out.println(voxel[0][0].get(GRB.DoubleAttr.X));
		double[][]getIntensity=new double[this.beams][this.aperture];
        for (int i = 0; i < this.beams; ++i) {
            for (int j = 0; j < this.aperture; ++j) {
            	getIntensity[i][j]=intensity[i][j].get(GRB.DoubleAttr.X);
            	//String varName="intensity"+i+"."+j;
            	//getIntensity[i][j]=model.getVarByName(varName).get(GRB.DoubleAttr.X);
            	//intensity[i][j] = model.addVar(minIntensity,maxIntensity,0.0, GRB.CONTINUOUS ,"Intensity" + i + "." + j);
                
            }
        }
		//GRBVar[] vars = model.getVars();
        newIntensity=getIntensity;
        objVal=model.get(GRB.DoubleAttr.ObjVal);
        
		
		
	}
	
	
	
	/**
	 * Setea el modelo que minimiza la penalizacion de todos los organos, esta funcion considera que un angulo puede tener mas aperturas que otro
	 * no tiene restricciones
	 * @throws GRBException
	 */
	public void setModelPenAllApertures() throws GRBException {
		this.model = new GRBModel(env);
		model.set(GRB.StringAttr.ModelName, "Direct Aperture Optimization");
		
		
		//set Variables
		
		// set variables for intensity an V
		GRBVar[][] intensity= new GRBVar[beams][];//intensity for beam,for aperture
		GRBVar[][] voxel=new GRBVar[R.length][];
		int indexi=0,indexj=0;
		
		for (int i = 0; i < this.beams; ++i) {
			intensity[i]=new GRBVar[sol.apertures.get(i).size()];
            for (int j = 0; j < sol.apertures.get(i).size(); ++j) {
            	indexi=i+1;
            	indexj=j+1;
            	intensity[i][j] = model.addVar(minIntensity,20,0.0, GRB.CONTINUOUS ,"Intensity" +"["+ indexi + "." + indexj+"]");
            	intensity[i][j].set(GRB.DoubleAttr.Start, sol.apertures.get(i).get(j).intensity); 
            }
        }
		
		
		for (int i = 0; i < this.R.length; ++i) {
			voxel[i]=new GRBVar[R[i]];
            for (int j = 0; j < R[i]; ++j) {
            	indexi=i+1;
            	indexj=j+1;
            	voxel[i][j] = model.addVar(0.0,1000.0,0.0, GRB.CONTINUOUS ,"v" + indexi + "[" + indexj+"]");
                
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
				key = keys.nextElement();
	            beams = aux_index.get(key);
	            
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

	            	for(int a=0;a<sol.apertures.get(beamIndex).size();a++) {
	            		coefficent = (double)sol.aperturesBmlts.get(beamIndex)[a][beamblet-diffBeamblets]*radiation;
	            		//coefficent = (double)aper[beamIndex][a][(beamblet-diffBeamblets)]*radiation;
	            		if(coefficent!=0 && o==0) {
							coefficent=coefficent*-1;
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][a] );
							//System.out.println(coefficent);
						}
						else {
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][a] );
						}
	            		
	            		
	            	}
	            }
	            if(o==0) {
					voxelRadiation.addConstant(eud[o]);
				}
				else{
					int constEud=eud[o]*-1;
					voxelRadiation.addConstant(constEud);
				}
				GRBLinExpr V=new GRBLinExpr();
				V.addTerm(1,voxel[o][count_voxel]);
				
				if(o==0) {
					model.addConstr(V, GRB.EQUAL, voxelRadiation, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
				}else {
					model.addConstr(V, GRB.GREATER_EQUAL, voxelRadiation, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
				}
				
				count_voxel++;
			}
		}
		
		//set model
		GRBQuadExpr objFunc= new GRBQuadExpr();
		for(int o = 0;  o< organs; o++) {
			double coef=(double)((weight[o]/R[o]));
			
			for (int j = 0; j < R[o]; ++j) {
				objFunc.addTerm(coef, voxel[o][j],voxel[o][j]);
				System.out.print("");
			}
		}

		model.setObjective(objFunc, GRB.MINIMIZE);
		model.update();
		writeModel();
        model.optimize();
        model.update();
       // model.computeIIS();
        //model.write("mod.ilp");
		double[][]getIntensity=new double[this.beams][];
        for (int i = 0; i < this.beams; ++i) {
        	getIntensity[i]=new double[sol.apertures.get(i).size()];
            for (int j = 0; j < sol.apertures.get(i).size(); ++j) {
            	getIntensity[i][j]=intensity[i][j].get(GRB.DoubleAttr.X);
            	//String varName="intensity"+i+"."+j;
            	//getIntensity[i][j]=model.getVarByName(varName).get(GRB.DoubleAttr.X);
            	//intensity[i][j] = model.addVar(minIntensity,maxIntensity,0.0, GRB.CONTINUOUS ,"Intensity" + i + "." + j);
                
            }
        }
		//GRBVar[] vars = model.getVars();
        newIntensity=getIntensity;
        objVal=model.get(GRB.DoubleAttr.ObjVal);
        
		
		
	}
	
	
	/**
	 * utiliza el modelo minimiza la penalizacion de los organos en riesgo y tiene como restricion el eud 0 para el PTV
	 * @throws GRBException
	 */
	public void setModelORPen() throws GRBException {
		this.model = new GRBModel(env);
		model.set(GRB.StringAttr.ModelName, "Direct Aperture Optimization");
		
		
		//set Variables
		
		// set variables for intensity an V
		GRBVar[][] intensity= new GRBVar[beams][aperture];//intensity for beam,for aperture
		GRBVar[][] voxel=new GRBVar[(R.length-1)][];
		int indexi=0,indexj=0;
		
		for (int i = 0; i < this.beams; ++i) {
            for (int j = 0; j < this.aperture; ++j) {
            	indexi=i+1;
            	indexj=j+1;
            	intensity[i][j] = model.addVar(minIntensity,20,0.0, GRB.CONTINUOUS ,"Intensity" +"["+ indexi + "." + indexj+"]");
            	intensity[i][j].set(GRB.DoubleAttr.Start, sol.apertures.get(i).get(j).intensity); 
            }
        }
		
		
		for (int i = 0; i < (this.R.length-1); ++i) {
			voxel[i]=new GRBVar[R[(i+1)]];
            for (int j = 0; j < R[(i+1)]; ++j) {
            	indexi=i+1;
            	indexj=j+1;
            	voxel[i][j] = model.addVar(0.0,80.0,0.0, GRB.CONTINUOUS ,"v" + indexi + "[" + indexj+"]");
                
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
				key = keys.nextElement();
	            beams = aux_index.get(key);
	            
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
	            	for(int a=0;a<aperture;a++) {
	            		coefficent = (double)sol.aperturesBmlts.get(beamIndex)[a][beamblet-diffBeamblets]*radiation;
	            		//coefficent = (double)aper[beamIndex][a][(beamblet-diffBeamblets)]*radiation;
	            		if(coefficent!=0 && o==0) {
							
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][a] );
							//System.out.println(coefficent);
						}
						else {
							voxelRadiation.addTerm(coefficent,intensity[beamIndex][a] );
						}
	            		
	            		
	            	}
	            }
	            if(o==0) {
				//	voxelRadiation.addConstant(eud[o]);
				}
				else{
					int constEud=eud[o]*-1;
					voxelRadiation.addConstant(constEud);
				}
				
				
				if(o==0) {
					model.addConstr(voxelRadiation, GRB.GREATER_EQUAL,eud[o] , "gEUD"+o+"["+(count_voxel+1)+"]");
				}else {
					GRBLinExpr V=new GRBLinExpr();
					V.addTerm(1,voxel[o-1][count_voxel]);
					model.addConstr(V, GRB.GREATER_EQUAL, voxelRadiation, "voxelRadiation"+o+"["+(count_voxel+1)+"]");
				}
				
				count_voxel++;
			}
		}
		
		//set model
		GRBQuadExpr objFunc= new GRBQuadExpr();
		for(int o = 1;  o< organs; o++) {
			double coef=(double)((weight[o]/R[o]));
			
			for (int j = 0; j < R[o]; ++j) {

				objFunc.addTerm(coef, voxel[o-1][j],voxel[o-1][j]);
				
				System.out.print("");
			}
		}

		model.setObjective(objFunc, GRB.MINIMIZE);
		model.update();
        //
		model.optimize();
        model.update();
       // writeModel();
        //model.computeIIS();
        model.write("modPen.lp");
		double[][]getIntensity=new double[this.beams][this.aperture];
        for (int i = 0; i < this.beams; ++i) {
            for (int j = 0; j < this.aperture; ++j) {
            	getIntensity[i][j]=intensity[i][j].get(GRB.DoubleAttr.X);
            	//String varName="intensity"+i+"."+j;
            	//getIntensity[i][j]=model.getVarByName(varName).get(GRB.DoubleAttr.X);
            	//intensity[i][j] = model.addVar(minIntensity,maxIntensity,0.0, GRB.CONTINUOUS ,"Intensity" + i + "." + j);
                
            }
        }
		//GRBVar[] vars = model.getVars();
        newIntensity=getIntensity;
        objVal=model.get(GRB.DoubleAttr.ObjVal);
        
		
		
	}
	
	public void writeModel() throws GRBException {
		model.write("out1.lp");
//		model.write("out1.mst");
		//model.write("out1.mps");
		//model.write("out1.sol");
		
	}
	
	
	public void cuadraticSum(TreatmentPlan sol) {
		
	}
	
	public void returnIntensity() {
		
	}
}
