package IMRT_BAO;

import java.io.IOException;

import SingleObj.DAO;


public class main {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String cerr="inputFile.txt";
		try {
			algorithm alg=new algorithm(cerr);
			int initialBac[]={0,70,140,210,280};
			alg.test(initialBac);
			//MODAO Prueba=new MODAO(cerr);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
