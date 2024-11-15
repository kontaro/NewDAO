package ApertureBase;

public class Decomposition {
	int [][] decomposition;;
	int alpha;
	
	public Decomposition(int [][] decomposition, int alpha){
		
		this.decomposition = copyMatrix(decomposition);
		this.alpha = alpha;
	}

	public int[][] getDecomposition() {
		return decomposition;
	}
	public int[][] copyMatrix(int [][] decomposition){
		int[][] dupA = new int [decomposition.length][decomposition[0].length];
		
		//int[][] auxMatrix = new int[decomposition.length][decomposition[0].length];
		for(int i=0 ; i<decomposition.length; i++){
			dupA[i] = decomposition[i].clone();
		}
		return dupA;
	}

	public void setDecomposition(int[][] decomposition) {
		this.decomposition = decomposition;
	}

	public int getAlpha() {
		return alpha;
	}

	public void setAlpha(int alpha) {
		this.alpha = alpha;
	}
}
