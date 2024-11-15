package ApertureBase;


import java.io.*;
import java.util.Arrays;
import java.util.Scanner;

import IMRT_DAO.Aperture;

import java.io.IOException;


import java.util.ArrayList;
import java.util.List;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;


public class MLC {
	int matrix=0;
	int beam_on_time=0;
	ArrayList<Aperture> decompositionAperture=new ArrayList<Aperture>();
	int maxIntensity=0;
	public 	MLC () {
		
	}

	void getAperture (int [][] A) {	
		
//		int [][] A = readFile("mat.txt");
		
		//System.out.println(Arrays.toString(matriz));
		

		int [][] a = new int [A.length][A[0].length+1];

		
		//Generate negative matrix
		for (int i = 0; i < A.length; i ++ ){
	      for(int j=0 ; j < A[0].length+1 ; j++ ){
	    	if(j == 0){
	    		a[i][j] = A[i][j]; 
	    	}else if(j == A[0].length){
	    		a[i][j] = -A[i][j-1];
	    	}else{
	    		a[i][j] = A[i][j] - A[i][j-1];
	    	}
	      }
		}
		ArrayList<Interval> [] intervals = generateIntervals(A, a);	
		ArrayList<Decomposition> decompositionMatrix = genetareDecompositionMatrix(A,intervals);
		getApertures(decompositionMatrix);
		matrix=decompositionMatrix.size();
		beam_on_time=0;
		for(int b=0;b<decompositionMatrix.size();b++) {
			
			beam_on_time=decompositionMatrix.get(b).alpha+beam_on_time;
			if(decompositionMatrix.get(b).alpha>maxIntensity)maxIntensity=decompositionMatrix.get(b).alpha;
		}
		
		
	}
	
	public void getApertures(ArrayList<Decomposition> decompositionMatrix) {
		decompositionAperture.clear();
		for(int i=0;i<decompositionMatrix.size();i++) {
			Aperture newAperture=new Aperture(convertToAperture(decompositionMatrix.get(i)),decompositionMatrix.get(i).getAlpha());
			decompositionAperture.add(newAperture);
		}
		
	}

	
	public int[][] convertToAperture(Decomposition decompositionMatrix){
		int[][] Aper=decompositionMatrix.getDecomposition();
		int[][] nAper= new int [Aper.length][2];
		int flag=0;
		for(int i=0;i<Aper.length;i++) {
			flag=0;
			nAper[i][0]=-1;
			nAper[i][1]=-1;
			for(int j=0;j<Aper[0].length;j++) {
				if(flag==0 && Aper[i][j]==1) {
					nAper[i][0]=j;
					flag=1;
				}else {
					if(flag==1 && Aper[i][j]==0) {
						nAper[i][1]=(j-1);
						flag=2;
					}
				}
				
			}
			if(flag==1) {nAper[i][1]=(Aper[0].length-1);}
		}
				
				
		
		return nAper;
		
	}
	
	public void imprimirMatriz(int[][] x){
		
		for(int i=0 ; i < x.length ; i++){
			for(int j=0 ; j< x[0].length ; j++){
				System.out.print(x[i][j] + " ");
			}
			System.out.println();
		}
	}
	
	//A = original matrix ; a = negative matrix 
	public ArrayList<Interval>[] generateIntervals(int[][] A, int[][] a ){
		
		int[][] dupA = new int [A.length][A[0].length];
		for(int i=0 ; i<A.length; i++){
			dupA[i] = A[i].clone();
		}
		
		ArrayList<Interval> [] intervals = new ArrayList [A.length];
		List<Integer> lb = new ArrayList<Integer>(); //lower bounds
		List<Integer> ub = new ArrayList<Integer>(); //upper bounds
		
		//Iteration for each row of matrix 
		for(int i= 0;  i < A.length ; i++){
			//Generation of lower and upper bounds by row  
			for(int j=0 ; j< a[0].length ; j++){
				if(a[i][j] > 0){
					lb.add(j);
				}else if(a[i][j] < 0){
					ub.add(j);
				}
			}
			
			int alpha;
			ArrayList<Interval> arrIntervals = new ArrayList<Interval>();
			//Generation of intervals for each row
			while(!emptyRow(dupA[i])){
				int [] pr ={lb.get(0), ub.get(0)};
				
				//Calculate alpha for interval
				if(a[i][lb.get(0)] < -(a[i][ub.get(0)])){
					alpha = a[i][lb.get(0)];
				}else{
					alpha = -a[i][ub.get(0)];
				}
				
				//Reduce intensities for dupA matrix 
				for(int z=lb.get(0); z<ub.get(0) ; z++){
					dupA[i][z] = dupA[i][z] - alpha;
				}
				
				//Remove lm or rm of lower bounds or upper bounds if they are equal to the value of a (negative) matrix
				a[i][lb.get(0)] -= alpha;
				a[i][ub.get(0)] += alpha;
				if(a[i][lb.get(0)] == 0){
					lb.remove(0);
				}
				if((-a[i][ub.get(0)]) == 0){
					ub.remove(0);
				}
				
				//Insert interval and alpha
				Interval interv = new Interval(pr , alpha);
				arrIntervals.add(interv);
				
			}	
			intervals[i] = arrIntervals;

		}
		return intervals;
	}
	
	//Generate decomposition of original matrix
	public ArrayList<Decomposition> genetareDecompositionMatrix(int [][] A, ArrayList<Interval> [] intervals){
		int[][] dupA = new int [A.length][A[0].length];
		int minAlpha;
		int[][] auxMatrix = new int[A.length][A[0].length];
		for(int i=0 ; i<A.length; i++){
			dupA[i] = A[i].clone();
		}	
		ArrayList<Decomposition> decompositionMatrix = new ArrayList<>();
		
		while(!emptyMatrix(dupA)){
			minAlpha = 10000;
			for(int i=0; i< intervals.length; i++){

				//Choose minimum alpha of intervals 
				if((intervals[i].size()>0) && (intervals[i].get(0).getAlpha() < minAlpha)){
					minAlpha = intervals[i].get(0).getAlpha();
				}
				
			}
			for(int i=0; i< intervals.length; i++){
				
				//Generate a decomposition matrix
				for(int j=0; j<dupA[0].length; j++){
					if((intervals[i].size() > 0) && (j >= intervals[i].get(0).getInterval()[0]) && (j < intervals[i].get(0).getInterval()[1])){
						auxMatrix[i][j] = 1;
					}else{
						auxMatrix[i][j] = 0;
					}
				}
				if((intervals[i].size()>0) && intervals[i].get(0).getAlpha() == minAlpha){
					intervals[i].remove(0);
				}else if((intervals[i].size()>0)){ 
					intervals[i].get(0).setAlpha(intervals[i].get(0).getAlpha() - minAlpha);
				}	
			}
			Decomposition matrix = new Decomposition(auxMatrix, minAlpha);
			decompositionMatrix.add(matrix);
			
			
			for(int i=0; i<dupA.length; i++){
				for(int j=0; j<dupA[0].length; j++ ){
					dupA[i][j] = dupA[i][j] - auxMatrix[i][j] * minAlpha;  
				}
			}
			
		}
		return decompositionMatrix;
	}
	
	//Validate if matrix is empty
	public boolean emptyMatrix(int[][] A){
		boolean empty = true;
		for(int i=0; i<A.length ; i++){
			if(!emptyRow(A[i])){
				empty = false;
			}
		}
		return empty;
	}
	
	//Validate if row is empty
	public boolean emptyRow(int[] A){
		boolean empty = true; 
		for(int i=0; i <A.length ; i++){
			if(A[i] != 0){
				empty = false;
			}
		}
		return empty;
	}
	
	//Leer fila
	public static int[][] readFile(String file) {
		try {
			File f = new File(file);
			Scanner s = new Scanner(f);//Lee numero de columnas
			FileReader fr = new FileReader(file);
			LineNumberReader lnr = new LineNumberReader(fr);//Lee numero de lineas
			int largo = 0;
			int alto = 0;
			while (lnr.readLine() != null){
	        	alto++;
	        	while (s.hasNextInt() && s.next()!="\r\n") {
					largo++;
					s.nextInt();
				}
	        }
			System.out.println(alto);
			System.out.println(largo);
			
			int[][] arr = new int[alto][largo];
			
			Scanner s1 = new Scanner(f);
			for (int i=0; i < alto; i++) {
				for(int j=0; j < largo; j++) {
					//System.out.println(s1.nextInt());
					arr[i][j] = s1.nextInt();
				}
			}
			System.out.println((arr[0][0]));
			return arr;
		}
		catch(Exception e) {
			return null;
		}
	}
	

}
