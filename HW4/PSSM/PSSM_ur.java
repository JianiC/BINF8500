import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.stream.Collectors;

public class PSSM_ur {

	public static void main(String[] args) throws FileNotFoundException{
	// Scanner for user input
	
	System.out.println("Please type in motif file name, sequence file name, scure cut off");
	Scanner user =new Scanner (System.in);
	
	String motif_f = user.next();
	String seq_f =user.next();
	double S= user.nextDouble();
	
		
		
		
		
		
	File input1= new File(motif_f);
	Scanner scnr = new Scanner(input1);	
	ArrayList<String> motif  = new ArrayList<String>();
	
	while(scnr.hasNext()){
		String n=scnr.nextLine();
		motif.add(n);		
	}
	
// frequency matrix	
	int l= motif.get(1).length();
	
	int matrix[][]=new int [l][4];
	
	for (int i=0;i<l;i++) {
		for(int n=0;n<motif.size();n++) {
			if (motif.get(n).charAt(i)=='a'){
					matrix[i][0]++;	
		}
			else if (motif.get(n).charAt(i)=='c'){
				matrix[i][1]++;	
			}
			else if (motif.get(n).charAt(i)=='g'){
				matrix[i][2]++;	
			}
			else if (motif.get(n).charAt(i)=='t'){
				matrix[i][3]++;	
			}		
		}		
	}
	
	System.out.println("frequency matrix:");
	System.out.println("pos.\tA\tC\tG\tT\t");
	for(int r=0;r<matrix.length;r++) {
		System.out.print(r+"\t");
		for(int c=0;c<matrix[0].length;c++) {
			System.out.print(matrix[r][c]+"\t");
		}
		System.out.println();
	}

//change to psedo-counts matrix and covert to score matrix
// covert to probabity matrix: sum for each position = 12 +0.25*4
	double p_matrix[][]=new double[l][4];
	for(int r=0; r<l; r++) {
		for(int c = 0; c < 4; c++) {
			p_matrix[r][c]=(matrix[r][c]+0.25)/13;
		}
	}
	
	System.out.println("\nprobablity matrix:");
	System.out.println("pos.\tA\tC\tG\tT\t");

	for(int r=0;r<matrix.length;r++) {
		System.out.print(r+"\t");
		for(int c=0;c<matrix[0].length;c++) {
			System.out.print(df3.format(p_matrix[r][c])+"\t");
		}
		System.out.println();
	}
	


	//score assigned
	// read sequence into array
	File seq = new File(seq_f ); // input file
	Scanner scnr2 = new Scanner(seq);	
	List<String> fasta = new ArrayList<String>();
	String line="";		
	while(scnr2.hasNext()) {
		line=scnr2.nextLine();
		fasta.add(line);					 
	}
	
	String name=fasta.get(0);
	fasta.remove(0);

	char[] sequence = fasta.stream().collect(Collectors.joining()).toCharArray();
	
	// get the complementary string
	
	char [] sequence_re = new char [sequence.length];
	for(int i=0; i<sequence.length;i++) {
		if (sequence[i]=='A') {
			sequence_re[i]='T';
		}
		else if (sequence[i]=='C') {
			sequence_re[i]='G';
		}
		else if (sequence[i]=='G') {
			sequence_re[i]='C';
		}
		else if (sequence[i]=='T') {
			sequence_re[i]='A';
		}
		
	}

	//get the sequence background
	
	int count_A=0;
	int count_C=0;
	int count_G=0;
	int count_T=0;
	for(int i=0;i<sequence.length;i++) {
		if (sequence[i]=='A') {
			count_A++;
		}
		else if (sequence[i]=='C') {
			count_C++;
		}
		else if (sequence[i]=='G') {
			count_G++;
		}
		else if (sequence[i]=='T') {
			count_T++;
		}
	}
	
	double q_A=(double)(count_A+count_T)/(2*sequence.length);
	double q_C=(double)(count_C+count_G)/(2*sequence.length);
	double q_G=(double)(count_G+count_C)/(2*sequence.length);
	double q_T=(double)(count_T+count_A)/(2*sequence.length);

	double s_matrix[][]=new double[l][4];
	for(int r=0; r<l; r++) {		
		s_matrix[r][0]=Math.log(p_matrix[r][0]/q_A);
		s_matrix[r][1]=Math.log(p_matrix[r][1]/q_C);
		s_matrix[r][2]=Math.log(p_matrix[r][2]/q_G);
		s_matrix[r][3]=Math.log(p_matrix[r][3]/q_T);
		
	}

	System.out.println("\nPSSM score:");
	System.out.println("pos.\tA\tC\tG\tT\t");
	for(int r=0;r<matrix.length;r++) {
		System.out.print(r+"\t");
		for(int c=0;c<matrix[0].length;c++) {
			System.out.print(df3.format(s_matrix[r][c])+"\t");
		}
		System.out.println();
	}
	
	System.out.print("\nStart\tEnd\tStrand\tScore\tSequence");
	 	
	for(int i=0;i<sequence.length-matrix.length;i++) {
		double sum=0;
		double sum_re=0;
		for(int k=0;k<matrix.length;k++) {
			if(sequence[i+k]=='A') {
				sum+=s_matrix[k][0];
				// for complementary chain
				sum_re+=s_matrix[k][3];
			}
			else if(sequence[i+k]=='C') {
				sum+=s_matrix[k][1];
				sum_re+=s_matrix[k][2];
			}
			else if(sequence[i+k]=='G') {
				sum+=s_matrix[k][2];
				sum_re+=s_matrix[k][1];
			}
			else if(sequence[i+k]=='T') {
				sum+=s_matrix[k][3];
				sum_re+=s_matrix[k][0];
			}
			
		
		}
		if(sum>S) {
		
			int end = i + matrix.length;
			System.out.print("\n"+i+"\t"+ end +"\t+\t"+df3.format(sum)+"\t");
			for (int n=0;n<l;n++) {
				System.out.print(sequence[i+n]);
			}	
		}
			
		if(sum_re>S) {
		
			int end = i + matrix.length;
			System.out.print("\n"+i+"\t"+ end +"\t-\t"+df3.format(sum_re)+"\t");
			for (int n=0;n<l;n++) {
				System.out.print(sequence_re[i+n]);
			}
			
		}
				
	}
	
	scnr.close();
	scnr2.close();

	}
	private static DecimalFormat df3 = new DecimalFormat("#.###");

}

