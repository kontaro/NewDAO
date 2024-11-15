package IMRT_DAO;

import java.io.Serializable;
import java.util.Comparator;

public class Aperture implements Serializable {

	private static final long serialVersionUID = 1L;
	public int[][] aperture;
	public double intensity;
	int flag;

	public Aperture(int[][] a, int i){
		this.aperture = a;
		this.intensity = i;
		flag=-1;
	}
	public Aperture(Aperture original){
		this.aperture = new int[original.aperture.length][original.aperture[0].length];
		this.intensity = original.intensity;
		
		for(int i=0;i<original.aperture.length;i++) {
			for(int j=0;j<original.aperture[0].length;j++) {
				aperture[i][j]= original.aperture[i][j];
			}
		}
		
		
		
	}
	
	public Aperture(int[][] a, float i) {
		this.aperture = a;
		this.intensity = i;
		flag=-1;
	}
	/**
	 * Se cambia la apertura por una nueva
	 * @param newAperture
	 */
	public void changeAperture(int[][] newAperture) {
		for(int i=0;i<aperture.length;i++) {
			for(int j=0;j<aperture[0].length;j++) {
				aperture[i][j]= newAperture[i][j];
			}
		}
	}
	/**
	 * Cambia una fila completa por un vector ingresado
	 * @param newRow tiene el vectopr que reemplazara la fila 
	 * @param index tiene la posicion de la fila
	 */
	public void changeRowAperture(int index,int[] newRow) {
		
		for(int j=0;j<aperture[0].length;j++) {
			aperture[index][j]= newRow[j];
		}
		
	}
	/**
	 * funcion que abre o cierra la apertura(no considera si existen movimientos infactibles)
	 * @param row 
	 * @param state si es 1 se abre la apertura, en caso contrario se cierra
	 * @param direction indica en que direccion se hace la operacion, si es true se hace en la izquierda, si es false en la derecha
	 */
	public void changeRow(int row, int state,boolean direction) {
		int aux=-1;
		if(state==1) aux=1;
		
		if(aperture[row][0]==-1 && aperture[row][1]==-1) {
			
		}else {
			if(direction==true) {
				aperture[row][0]=aperture[row][0]+aux;
			}
			else {
				aperture[row][1]=aperture[row][1]+aux;
			}
		}
	}
	
	/**
	 * funcion que abre o cierra la apertura(no considera si existen movimientos infactibles)
	 * @param row 
	 * @param state si es 1 se abre la apertura, en caso contrario se cierra
	 * @param direction indica en que direccion se hace la operacion, si es true se hace en la izquierda, si es false en la derecha
	 *
	 * TODO:
	 * 
	 * - modificar para que cierre hoja
	 * - que pasa si cierra algo ya cerrado?
	 * - hacer un if para casos de aperture y para casos de cerrado, es mas facil asi ver si 
	 */
	public int changeRow(int row, int state,boolean direction,int leftLimit,int rightLimit) {
		int aux=-1;
		int lastPositionChanged=-1;
		if(state==1) aux=1;
		int middle=((rightLimit+1)/(leftLimit+1))-1;
		
		//Apertura de hoja
		if(aperture[row][0]==-1 && aperture[row][1]==-1) {
			aperture[row][0]=middle;
			aperture[row][1]=middle;
			lastPositionChanged=aperture[row][0];

			return lastPositionChanged;
		}else {
			//Apertura de hoja por lado izquierdo
			if(direction==true) {
				aperture[row][0]=aperture[row][0]+aux;
				lastPositionChanged=aperture[row][0];
				if(aux==1)
					return lastPositionChanged;
				else
					return lastPositionChanged+1;
			}
			//Apertura de hoja por lado derecho
			else {
				aperture[row][1]=aperture[row][1]+aux;
				lastPositionChanged=aperture[row][1];
				if(aux==1)
					return lastPositionChanged;
				else
					return lastPositionChanged+1;
			}
		}
	}
	
	
	public int changeRow2(int row, int state,boolean direction,int leftLimit,int rightLimit) {
		int aux=-1;
		int lastPositionChanged=-1;
		
		int middle=((rightLimit+1)/(leftLimit+1))-1;
		
		//Separamos los casos para la apertura y cierre de hojas
		if(state==1) {
			if(aperture[row][0]==-1 && aperture[row][1]==-1) {
				//-1,-1  -> X,X
				aperture[row][0]=middle;
				aperture[row][1]=middle;
				lastPositionChanged=aperture[row][0];
				flag=1;
				return lastPositionChanged;
			}else {
				//Apertura de hoja por lado izquierdo
				if(direction==true) {
					//5,7  -> 4,7
					aux=-1;
					aperture[row][0]=aperture[row][0]+aux;
					lastPositionChanged=aperture[row][0];
					flag=2;
					return lastPositionChanged-1;
				}
				//Apertura de hoja por lado derecho
				else {
					//5,7  -> 5,8
					aux=1;
					aperture[row][1]=aperture[row][1]+aux;
					lastPositionChanged=aperture[row][1];
					flag=3;
					return lastPositionChanged;

				}
			}
		}
		
		//en el caso de que se cierre se hace lo sigueinte
		else {
			if((aperture[row][0]==-1)  && (-1 ==aperture[row][1])) {
				flag=11;
				return -3;
			}
			
			
			if(aperture[row][0]==aperture[row][1]) {
				//6,6 -> 7,6 -> -1,-1
				//6,6  -> 6,5 -> -1,-1
				lastPositionChanged=aperture[row][0];
				aperture[row][0]=-1;
				aperture[row][1]=-1;
				
				flag=4;
				return lastPositionChanged;
			}else {
				if(direction==true) {
					// lado izquierdo
					//5,7  -> 6,7
					
					aux=1;
					lastPositionChanged=aperture[row][0];
					aperture[row][0]=aperture[row][0]+aux;
					//lastPositionChanged=aperture[row][0];
					flag=5;
					return lastPositionChanged;
					
				}
				//Apertura de hoja por lado derecho
				else {
					//5,7  -> 5,6
					
					aux=-1;
					aperture[row][1]=aperture[row][1]+aux;
					lastPositionChanged=aperture[row][1];
					flag=6;
					return lastPositionChanged+1;
				}
			}
		}
	}
	
	
	
	public void printAperture() {
		for(int i=0;i<aperture.length;i++) {
			for(int j=0;j<aperture[0].length;j++) {
				System.out.print(aperture[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	public int getTotalApertureOpen() {
		int total=0;
		for(int i=0;i<aperture.length;i++) {
			for(int j=0;j<aperture[0].length;j++) {
				if(aperture[i][j]==1)
					total=total+1;
			}
		}
		return total;
	}
	
	/**
	 * La funcion combrueba que las filas sean factibles.
	 * Esto debido a que como no consideramos ninguna restriccion del mlc solo aqui se puede producir infactibilidad
	 * Funcion pensada para los movimientos de cruzamiento (no deberian llegar casos de -1 en un lado ya que es cruzamiento)
	 * @return
	 */
	
	public int isFeasible() {
		
		for(int i=0;i<aperture.length;i++) {
			if(aperture[i][0]>aperture[i][1] || aperture[i][1]<aperture[i][0])
				return 1;
		}
		return 0;
	}
	public static Comparator<Aperture> intensityComparator = new Comparator<Aperture>() {
        
	      public int compare(Aperture s1, Aperture s2){
	  


	                  // ascending order
	                  return Double.valueOf(s2.intensity).compareTo(Double.valueOf(s1.intensity));
	  
	              }

	
	          };
	
}
