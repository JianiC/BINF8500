// BINF 8500 HW 3
// input file and scores are hardcodes
//sequence , alignment score and alignment sequence will be print out at console




import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.stream.Collectors;

public class alignment2 {
	
	public static void main(String[] args) throws FileNotFoundException{
		File text1 = new File("HIV1a.fasta");	// input file			
		Scanner scnr1 = new Scanner(text1);	
		List<String> fasta1 = new ArrayList<String>();
		String line1="";		
		while(scnr1.hasNext()) {
			line1=scnr1.nextLine();
			fasta1.add(line1);					 
		}
		
		for(int i=0;i<fasta1.size();i++) {
			System.out.print(fasta1.get(i));
			System.out.println();
		}
		String name1=fasta1.get(0);
		fasta1.remove(0);
		fasta1.add(0," ");
		char[] f1 = fasta1.stream().collect(Collectors.joining()).toCharArray();
				
		File text2 = new File("HIV1b.fasta"); // input file
		Scanner scnr2 = new Scanner(text2);	
		List<String> fasta2 = new ArrayList<String>();
		String line2="";		
		while(scnr2.hasNext()) {
			line2=scnr2.nextLine();
			fasta2.add(line2);					 
		}
		
		System.out.println();
		
		for(int i=0;i<fasta2.size();i++) {
			System.out.print(fasta2.get(i));
			System.out.println();
		}
		
		String name2=fasta2.get(0);
		fasta2.remove(0);
		fasta2.add(0," ");
		char[] f2 = fasta2.stream().collect(Collectors.joining()).toCharArray();
		
		
		
		
		System.out.println();

		
		
		
		// initialize match, mismatch and gap score
		int gap = -1;
		int mismatch=0;
		int match=1;
		int dia,hori,vert,best;
		
		
		int [][]matrix=new int[f1.length][f2.length];
		// fill first row and column
		for(int r=0;r<matrix.length;r++) {
			matrix[r][0]=gap*r;
		}
	
		for(int c=0;c<matrix[0].length;c++) {
			matrix[0][c]=gap*c;
		}
		
		// fil the rest of the matrix
		for(int r=1;r<matrix.length;r++) {
			for(int c=1;c<matrix[0].length;c++) {
				if(f1[r]==f2[c]) {
					dia=matrix[r-1][c-1]+match;					
				}
				else {
					dia=matrix[r-1][c-1]+mismatch;
				}
				
				hori=matrix[r][c-1]+gap;
				vert=matrix[r-1][c]+gap;
				
				best=dia;
				if (hori>=best) {
					best=hori;
				}
				else if(vert>=best) {
					best=vert;
				}
				matrix[r][c]=best;
			}
		}
		
		System.out.println();
		

		System.out.println("Best score: "+matrix[matrix.length-1][matrix[0].length-1]);
		
		
		
		//trace back to find the right alignment
		
		ArrayList<Character>seq1 = new ArrayList<Character>();
		ArrayList<Character>seq2 = new ArrayList<Character>();
		ArrayList<String>align = new ArrayList<String>();
		
		int r=matrix.length-1;
		int c=matrix[0].length-1;
		
		while( r>0 || c > 0) {
			if(r >0 && c >0) {					
				if(matrix[r][c]==matrix[r-1][c-1]+match) {
					seq1.add(0,f1[r]);
					seq2.add(0,f2[c]);
					align.add(0,"*");
					
				
					r--;
					c--;
				}
				
				
				
				else if (matrix[r][c]==matrix[r-1][c]+gap) {
					seq1.add(0,f1[r]);
					seq2.add(0,'-');
					align.add(0," ");
					r--;
				}
				
				else if(matrix[r][c]==matrix[r][c-1]+gap) {
					seq1.add(0,'-');
					seq2.add(0,f2[c]);
					align.add(0," ");
					c--;
				}
				else if(matrix[r][c]==matrix[r-1][c-1]+mismatch) {
					seq1.add(0,f1[r]);
					seq2.add(0,f2[c]);
					align.add(0," ");
					r--;
					c--;
				}
			}
			
			else if(r==0 && c>0) {
				seq1.add(0,'-');
				seq2.add(0,f2[c]);
				align.add(0," ");
				c--;
			}
			
			else {
				seq1.add(0,f1[r]);
				seq2.add(0,'-');
				align.add(0," ");
				r--;
			}				
		}
		
		

		for(int n=1;n<seq1.size()/60;n++) {
			System.out.println("\n");
			int index=((n-1)*60)+1;
			String sindex=String.format("%4d", index);
			System.out.print("\t  "+sindex+": ");
			for(int i=(n-1)*60;i<n*60;i++) {
				System.out.print(seq1.get(i));				
			}
			System.out.println();
			System.out.print("\t\t");
			for(int i=(n-1)*60;i<n*60;i++) {
				System.out.print(align.get(i));				
			}
			System.out.println();
			System.out.print("\t  "+sindex+": ");
			for(int i=(n-1)*60;i<n*60;i++) {
				System.out.print(seq2.get(i));				
			}
			
		}
		
		int count=seq1.size()/60;
		int index2=60*(count-1)+1;
		String sindex2=String.format("%4d", index2);
		System.out.println("\n");
		System.out.print("\t  "+sindex2+": ");
		for(int i=60*count;i<seq1.size();i++ ) {
			System.out.print(seq1.get(i));
		}
		System.out.println();
		System.out.print("\t\t");
		for(int i=60*count;i<seq1.size();i++ ) {
			System.out.print(align.get(i));
		}
		System.out.println();
		System.out.print("\t  "+sindex2+": ");
		for(int i=60*count;i<seq1.size();i++ ) {
			System.out.print(seq2.get(i));
		}

		scnr1.close();
		scnr2.close();


	}

}
