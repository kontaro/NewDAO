package IMRT_BAO;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import IMRT_DAO.TreatmentPlan;
import SingleObj.DAO;

public class algorithm {
	
	String inputFile;
	public algorithm(String inputFile){
		this.inputFile=inputFile;
	}
	
	
	public ArrayList<TreatmentPlan> ParetoLocalSearch(int initialBACs[]) throws IOException {
		
		ArrayList<TreatmentPlan> spTPlan=new ArrayList<TreatmentPlan>();
		ArrayList<TreatmentPlan> finalTP=new ArrayList<TreatmentPlan>();
		spTPlan.addAll(FirstPhase(initialBACs));
		for (int i = 0; i < spTPlan.size(); i++) {
			spTPlan.get(i).printValues();
		}
		finalTP.addAll(SecondPhase(spTPlan));
		for (int i = 0; i < spTPlan.size(); i++) {
			finalTP.get(i).printValues();
		}
		
		
		
		return finalTP;
		
	}
	
	public ArrayList<TreatmentPlan> test(int initialBACs[]) throws IOException {
		ArrayList<TreatmentPlan> Frontnew=new ArrayList<TreatmentPlan>();
		int initialBAC4[]= {0,70,145,210,280};
		int initialBAC6[]= {0,70,140,210,285};
		int initialBAC7[]= {355,75,140,210,280};
		int initialBAC8[]= {355,70,145,210,280};
		ArrayList<double[]> w=new ArrayList<double[]>();
		w.addAll(weightGeneratorTestMO());
		
		
		Frontnew.addAll(WeightedSum(initialBAC4, w));
		
		for (int i = 0; i < Frontnew.size(); i++) {
			Frontnew.get(i).printValues();
		}
		
		Frontnew.clear();
		Frontnew.addAll(WeightedSum(initialBAC6, w));
		
		for (int i = 0; i < Frontnew.size(); i++) {
			Frontnew.get(i).printValues();
		}
		Frontnew.clear();
		Frontnew.addAll(WeightedSum(initialBAC7, w));
		
		for (int i = 0; i < Frontnew.size(); i++) {
			Frontnew.get(i).printValues();
		}
		Frontnew.clear();
		Frontnew.addAll(WeightedSum(initialBAC8, w));
		
		for (int i = 0; i < Frontnew.size(); i++) {
			Frontnew.get(i).printValues();
		}
		return Frontnew;
	}
	
	/**
	 * This function made the N1 movement of the work of maikol
	 * @param solution
	 * @return
	 */
	public ArrayList<int[]> generateBACNeigbourhood(int[] solution) {
		ArrayList<int[]> neighbours=new ArrayList<int[]>();
		for (int i = 0; i < solution.length; i++) {
			int[] newNeigbourPlus= new int[solution.length];
			int[] newNeigbourLess= new int[solution.length];
			System.arraycopy(solution, 0, newNeigbourPlus, 0, solution.length);
			System.arraycopy(solution, 0, newNeigbourLess, 0, solution.length);
			
			
			if(newNeigbourPlus[i]==355) {
				newNeigbourPlus[i]=0;
			}else {
				newNeigbourPlus[i]=newNeigbourPlus[i]+5;
			}
			if(newNeigbourLess[i]==0) {
				newNeigbourLess[i]=355;
			}else {
				newNeigbourLess[i]=newNeigbourLess[i]-5;
			}
				
			
			if(areAllElementsDistinct(newNeigbourPlus))neighbours.add(newNeigbourPlus);
			if(areAllElementsDistinct(newNeigbourLess))neighbours.add(newNeigbourLess);
		}
		
		return neighbours;
		
	}
	// Método para verificar si todos los elementos son distintos
    public  boolean areAllElementsDistinct(int[] array) {
        Set<Integer> seen = new HashSet<>();  // Crear un conjunto para almacenar los números

        for (int num : array) {
            if (!seen.add(num)) {  // Intenta agregar el número al conjunto
                return false;  // Si el número ya estaba, devuelve falso
            }
        }
        return true;  // Si no se encontraron duplicados, devuelve verdadero
    }
 // Método para verificar si dos arrays son iguales
    private boolean arraysAreEqual(int[] array1, int[] array2) {
        return Arrays.equals(array1, array2);
    }
    
    
 // Método para eliminar los arrays de 'vecinos' que están en 'visitados'
    public void removeVisited(ArrayList<int[]> vecinos, ArrayList<int[]> visitados) {
        Iterator<int[]> iterator = vecinos.iterator();
 
        while (iterator.hasNext()) {
            int[] vecino = iterator.next();
            boolean isVisited = false;
            
            for (int[] visitado : visitados) {
                if (arraysAreEqual(vecino, visitado)) {
                    isVisited = true;
                    break;
                }
            }

            if (isVisited) {
                iterator.remove();
            }
        }
    }
    
	public TreatmentPlan samplePoint(int[] initialBACs) throws IOException {
		double w[]={0,1,1};
		DAO sPDAO=new DAO(inputFile,initialBACs,w);// call the dao proces to objatin a sample point to the bac
		// save the vlues of the sample point
		return sPDAO.sol;
		
		
	}
	
	public ArrayList<TreatmentPlan> WeightedSum(TreatmentPlan Tp, ArrayList<double[]> w) throws IOException {
		int initialBACs[]=Tp.getBAC();
		ArrayList<TreatmentPlan> newList=new ArrayList<TreatmentPlan>();
		for (double[] i : w) {
			DAO sPDAO=new DAO(inputFile,initialBACs,i);
			newList.add(sPDAO.sol);
			
		}
		// call the dao proces to objatin a sample point to the bac
		// save the vlues of the sample point
		return newList;
		
		
	}
	public ArrayList<TreatmentPlan> WeightedSum(int[] initialBACs, ArrayList<double[]> w) throws IOException {
		
		ArrayList<TreatmentPlan> newList=new ArrayList<TreatmentPlan>();
		for (double[] i : w) {
			DAO sPDAO=new DAO(inputFile,initialBACs,i);
			newList.add(sPDAO.sol);
			
		}
		// call the dao proces to objatin a sample point to the bac
		// save the vlues of the sample point
		return newList;
		
		
	}
	
	
    // Método para verificar si 'this' domina a otro TreatmentPlan
    public boolean dominates(TreatmentPlan thisTP,TreatmentPlan otherTP) {
        return (thisTP.scores[0] <= otherTP.scores[0] && thisTP.scores[1] <= otherTP.scores[1]) &&
               (thisTP.scores[0] < otherTP.scores[0] || thisTP.scores[1] < otherTP.scores[1]);
    }
	
	
	// Método para actualizar la frontera de dominancia
    public  void updateDominanceFrontier(ArrayList<TreatmentPlan> frontier, TreatmentPlan newPlan) {
        boolean isDominated = false;
        Iterator<TreatmentPlan> iterator = frontier.iterator();

        while (iterator.hasNext()) {
            TreatmentPlan plan = iterator.next();
            // Si el nuevo plan es dominado por uno ya en la frontera, no necesita ser añadido
            if (dominates(plan,newPlan)) {
                isDominated = true;
                break;
            }
            // Si el nuevo plan domina a alguno en la frontera, ese plan debe ser eliminado
            if (dominates(newPlan,plan)) {
                iterator.remove();
            }
        }

        // Si el nuevo plan no es dominado por ninguno en la frontera, debe ser añadido
        if (!isDominated) {
            frontier.add(newPlan);
        }
    }
	
	public ArrayList<double[]> weightGeneratorTestMO(){
		ArrayList<double[]> weights=new ArrayList<double[]>();
		double a[]={0,5,1};// +PTV
		double c[]={0,3,1};//
		double g[]={0,1,3};//
		double i[]={0,1,5};// +Rectum
		weights.add(a);
		weights.add(c);
		weights.add(g);
		weights.add(i);

		return weights;
	}
	public ArrayList<TreatmentPlan> FirstPhase(int[] initialBACs) throws IOException {
		ArrayList<int[]> visitedBACs=new ArrayList<int[]>();
		ArrayList<int[]> toVisit=new ArrayList<int[]>();
		ArrayList<TreatmentPlan> fronList=new ArrayList<TreatmentPlan>();
		toVisit.add(initialBACs);
		fronList.add(samplePoint(initialBACs));
		int iter=0;
		while(!toVisit.isEmpty()&&iter<=7) {
			System.out.println(toVisit.size());
			ArrayList<int[]> newNeigbours=generateBACNeigbourhood(toVisit.remove(0));
			removeVisited(newNeigbours, visitedBACs);
			removeVisited(newNeigbours, toVisit);
			toVisit.addAll(newNeigbours);
			
			
			fronList.add(samplePoint(toVisit.remove(0)));
			iter++;
			
		}
		

		
		return fronList;
		
	}
	
	public ArrayList<TreatmentPlan> SecondPhase(ArrayList<TreatmentPlan> Front) throws IOException {
		ArrayList<TreatmentPlan> Frontnew=new ArrayList<TreatmentPlan>();
		ArrayList<double[]> w=new ArrayList<double[]>();
		w.addAll(weightGeneratorTestMO());
		for (TreatmentPlan treatmentPlan : Front) {
			Frontnew.addAll(WeightedSum(treatmentPlan, w));
		}
		
		return Frontnew;
		
	}
	

}
