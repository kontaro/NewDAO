package SingleObj;

import java.io.IOException;


import MultiObj.MODAO;

public class main {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String cerr="inputFile.txt";
		try {

			DAO Prueba=new DAO(cerr);
			//MODAO Prueba=new MODAO(cerr);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
