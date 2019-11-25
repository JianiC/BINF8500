/* Name: Gibbs
 * Author: Jiani Chen
 * Date: 11/18/19
 * Usage" gibbs InputFile(sequences in Fasta format) MotifLength(estimated length of motif)
 * this codes was set as 50 seeds, 200 times update, adjust motif length every 5 times update
 * 
 */
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

public class gibbs {
	public static void main(String[] args) throws FileNotFoundException{
		File seq = new File(args[0]);	// input file			
		Scanner scnr = new Scanner(seq);	
		List<String> seq1 = new ArrayList<String>();
		List<String> taxa = new ArrayList<String>();
		
		while(scnr.hasNext()) {
			String line = scnr.nextLine().trim();
			
			if (line.charAt(0)=='>') {
				taxa.add(line);
			}
			else {
				seq1.add(line);
			}		
		}
		
		int seq_size = seq1.size()/taxa.size(); // get the line num for single sequence				
		String[] seq2 = new String[taxa.size()];
		
		for(int i=0;i<seq2.length;i++) {
			StringBuilder sb= new StringBuilder();
			int n= seq_size*i;
			for(int j=0; j<seq_size;j++) {
				sb.append(seq1.get(n+j));
			}
			String sequence= sb.toString().toUpperCase();
			seq2[i]= sequence;
		
		}
		
		/*
		for(int i=0;i<seq2.length;i++) {
			System.out.println(taxa.get(i));
			System.out.println(seq2[i]);
		}		
		*/
		
	
		// initialize the length m
		int m = Integer.valueOf(args[1]);
		System.out.println("Initial motif length: "+ m);
		
	
		/*
		

			// function test
		  
		  int[] test = new int[]{31,32,42,36,53,3,52,48,32,48,43,45,72,11,39,29};
		  
		  // get motif
		  List<String> test2 = new ArrayList<String>();
		  for(int i=0;i<taxa.size();i++){
			  test2.add(get_motif2(test[i],m,seq2[i]));
			  int start =test[i]+1;
			  int end =start+m;
			  
			  
			  System.out.println(test2.get(i)+ "\t"+start+"-"+end+"\t"+taxa.get(i));
		  }
		  
		  // print the global pssm score for this motif
		  //temp_pssm=total_pssm(ind_new, m ,seq2);
			  
		  double a = total_pssm(test,m,seq2);
		  
		  
		  
		  System.out.println("score for test is : "+ a);
		  
		*/
		
		int t=50;// how many seeds you want to try
		int seed=0;
		// to store the parameters from different seeds
		List<Double> pssm_max=new ArrayList<Double>();
		int[][]global_ind=new int[t+1][seq2.length];
		List<Integer>length=new ArrayList<Integer>();
				
		while(seed<t) {
			seed++;
		
		//seed
		
		//for all of the sequence get the  index and put them into a list		
		int ind[]=new int[taxa.size()];
		
		for(int i=0;i<seq2.length;i++) {
			Random random=new Random();
			ind[i] = random.nextInt(seq2[i].length()-m);			
		}
		
		//Update
		double temp_pssm;
		double temp_pssmadd;
		double temp_pssmminus;
		double max_pssm=-10000;
		int max_length=m;
		int ind_max[]=new int[taxa.size()];
		int count=0;
		while(count<100) {
			count++;
			int ind_new[]=new int[taxa.size()];
			 int ind_add[]=new int[taxa.size()];
			  int ind_minus[]=new int[taxa.size()];
			
			// from the current index, you want to generate a new index based on previous propotion motif 
			List<String> motif_update = new ArrayList<String>();
			for(int i=0; i< taxa.size();i++) {				
				  motif_update.add(get_motif2(ind[i],m,seq2[i]));  			
			}
			
			for(int i=0;i<ind.length;i++) { //loop through all seq, i is the seq number
				
				List<String> motif_update2 = new ArrayList<String>();		  			
				 for(int n=0; n< motif_update.size();n++) {
						if(n!=i) {
							motif_update2.add(motif_update.get(n));
						}
				  }				 				
				//convert it to array
				String [] motif_u= new String [taxa.size()-1];
				motif_u= motif_update2.toArray(motif_u);
				
				ind_new[i]= index(motif_u,seq2[i], m);

			}
			
			// get the pssm score for the new motif we generate	
			
				temp_pssm=total_pssm(ind_new, m ,seq2);
				//System.out.println("\nthe new pssm is "+temp_pssm);
				
				ind=ind_new.clone();
				/*	
				if (temp_pssm>max_pssm) {
					max_length=m;
					max_pssm=temp_pssm;
					ind_max=ind_new.clone();
					 //System.out.println(" pssm is " +max_pssm+" length "+ max_length);
					}
				
				*/
				// adjust the length of the motif with every 5 times update			
				if(count % 5==0) {
					// increase motif length 1
					temp_pssmadd=-1000;
					
					// covert ind into a list to test if the motif length could be increase by index -1				
					List<Integer> ind_test=new ArrayList<Integer>();
					for (int r=0;r<ind_new.length;r++) {
						ind_test.add(ind_new[r]);
					}
					
					if (ind_test.contains(0)==false) {
						// get the pssm by motif length +1
						for (int r=0;r<ind_new.length;r++) {
							ind_add[r]=ind_new[r]-1;
						}
						temp_pssmadd=total_pssm(ind_add,m+1,seq2);					
					}
					
					// decrease motif length by 1			
					temp_pssmminus=-1000;
					// check the motif length could be dicrease need m-1>0
					if(m>1) {
						for (int r=0;r<ind_new.length;r++) {
							ind_minus[r]=ind_new[r]+1;
						}
											
						temp_pssmminus=total_pssm(ind_minus,m-1,seq2);
						
					}
					
								
					// compare with these three situation to adjust motif length
					
					//situation 1
					if (temp_pssmadd>temp_pssm) {
						m=m+1;
						ind=ind_add.clone(); // the index of motif for next update
					

						if (temp_pssmadd>max_pssm) { // get the max for the all updating process
							max_length=m;
							max_pssm=temp_pssmadd;
							ind_max=ind_add.clone();
							//System.out.println(" pssm is " +max_pssm+" length "+ max_length);
						}
						
					}
				
					
					// situation2
					else if (temp_pssmminus>temp_pssm) {
						m=m-1;
						ind=ind_minus.clone();

						if (temp_pssmadd>max_pssm) {
							max_length=m;
							max_pssm=temp_pssmminus;
							ind_max=ind_minus.clone();
							//System.out.println(" pssm is " +max_pssm+" length "+ max_length);
						}
						
					}
					
					
					// the same as when count % 5 !=0
					else if (temp_pssm>=temp_pssmadd || temp_pssm>=temp_pssmminus) {
						ind=ind_new.clone();// motif length does not change
						
						if (temp_pssm>max_pssm) {
							max_length=m;
							max_pssm=temp_pssm;
							ind_max=ind_new.clone();
							//System.out.println(" pssm is " +max_pssm+" length "+ max_length);
						}					
					}			
				}
				
				
				
				else if (count %5!=0) { //skip the adjustment of motif length
				
					
					ind=ind_new.clone();
					if (temp_pssm>max_pssm) {
						max_length=m;
						max_pssm=temp_pssm;
						ind_max=ind_new.clone();
						//System.out.println(" pssm is " +max_pssm+" length "+ max_length);
					}
					
				}
				
				

			// loop through the update process and get the maximum value from updat
				
				
				
				
				
				
				
				
				
				
				
				
				
	
		}
		//store the maximal from each seed
		 pssm_max.add(max_pssm);
		 length.add(max_length);
		 for(int j=0;j<ind_max.length;j++) {
			 global_ind[seed][j]=ind_max[j];
		 }
		
		
			}
			 
	//  get global maximal from different seeds
		for(int i=0;i<t;i++) {
			//System.out.println(i+"\t"+pssm_max.get(i));
		}
		
		
	
				int max_seed=0;  
				double global_max=pssm_max.get(0);
				for(int i=0;i<t;i++) {
					if(pssm_max.get(i)>global_max) {
						global_max=pssm_max.get(i);
						max_seed=i;
					}
				}
				
				int motif_length=length.get(max_seed);
				
				System.out.println(max_seed);
				System.out.println("Final motif length: "+ motif_length);
				System.out.println("Final score:\t "+ global_max);
				System.out.println("\nMotif sequences and locations:");
				// get motif
				List<String> global_motif = new ArrayList<String>();
				  for(int i=0; i< taxa.size();i++) {				
					  global_motif.add(get_motif2(global_ind[max_seed][i],length.get(max_seed),seq2[i])); 
					  int start = global_ind[max_seed][i];
					  int end= start + motif_length;
					  System.out.println(global_motif.get(i)+"\t"+start +"-"+end+"\t"+taxa.get(i));

		}
		  
		  
		  
		  
		
		scnr.close();
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// calculate the log2
	public static double log2(double x) {
		return(double)(Math.log(x)/Math.log(2));
	}
	
	
	
	
	//extract random sequence helper	
	public static String get_motif(int m, String seq) {
		// generate start index randomly, m is the length of motif
		Random random=new Random();
		int start = random.nextInt(seq.length()-m);
		String motif = seq.substring(start, start+m);
		
		return motif;
	}
	
	// extract sequence helper2 for updating step
	public static String get_motif2(int start, int m, String seq) {
		// generate start index randomly, m is the length of motif

		String motif = seq.substring(start, start+m);		
		return motif;
	}
	
	
	
	
	
	
	
	
	
	// PSSM function
	public static double[][] pssm  (int m, String[] motif, String seq) {
		
		// get frequency matrix
		int f_matrix[][]=new int[m][4];

		for(int i=0; i< m;i++) {
			for(int j=0;j<motif.length;j++) {
				if(motif[j].charAt(i)=='A') {
					f_matrix[i][0]++;
				}
				else if (motif[j].charAt(i)=='C') {
					f_matrix[i][1]++;
				}
				else if (motif[j].charAt(i)=='G') {
					f_matrix[i][2]++;
				}
				else if (motif[j].charAt(i)=='T') {
					f_matrix[i][3]++;
				}
			}
		}
		
		//add pseudo-counts to probability matrix
		double p_matrix[][]=new double[m][4];
		double total = 0.25*4 + motif.length ; // each column add 0.25 to the total nucleotide at each position in motif.txt
		for(int i=0; i< m; i++) {
			for(int j=0; j<4; j++) {
				p_matrix[i][j]=(f_matrix[i][j]+0.25)/ total;
			}
		}

		
		//get GC background from the sequence
		int GC =0;
		for(int i=0;i<seq.length();i++) {
			if (seq.charAt(i)=='G'||seq.charAt(i)=='C') {
				GC++;
			}
		}
		
		double q_GC = (double)GC/(2*seq.length());
		double q_AT = (double)(1 - 2* q_GC)/2;
		
		//pssm matrix
		double pssm[][]=new double [m][4];
		
		for(int i=0; i< m;i++) {
			pssm[i][0]=log2(p_matrix[i][0]/q_AT);
			pssm[i][1]=log2(p_matrix[i][1]/q_GC);
			pssm[i][2]=log2(p_matrix[i][2]/q_GC);
			pssm[i][3]=log2(p_matrix[i][3]/q_AT);
			
		}
	
		return pssm;
		
	
	}
	

	public static int index(String[]motif,String seq, int m) {
		double pssm[][]=pssm(m,motif,seq);
		
		//assign score to each position in seq and return the position 
		int max=0; // max is the start index with the motif has the largest pssm
		double sum []=new double[seq.length()-m];
		for(int i=0;i<seq.length()-m;i++) {	
			double temp_sum=0;
			for(int j=0;j<m;j++) {
				if(seq.charAt(i+j)=='A') {
					temp_sum+=pssm[j][0];
				}
				if(seq.charAt(i+j)=='C') {
					temp_sum+=pssm[j][1];
				}
				if(seq.charAt(i+j)=='G') {
					temp_sum+=pssm[j][2];
				}
				if(seq.charAt(i+j)=='T') {
					temp_sum+=pssm[j][3];
				}			
			}
			sum[i]=temp_sum;		
		}
		
		//covert this sum list into a propotion list
		
		
		double sum_p[]=new double[sum.length];
		for(int i=0; i< sum_p.length;i++) {
				
			sum_p[i]=Math.pow(2,sum[i]);
		}
		
		// convert this sum positive list into a propotion score list
		
		double p_score[]=new double[sum_p.length];
		for(int i=0;i<p_score.length;i++) {
			for(int j=0;j<=i;j++) {
				p_score[i]+=sum_p[j];
			}
		}
		
		
	
		//Get the real propotion list
		// get the total value of p_score list
		double total=0;
		for (int i=0; i<p_score.length;i++) {
			total+=sum_p[i];
		}
		
		double p[]=new double[p_score.length];
		for(int i=0;i<p.length;i++) {
			p[i]=p_score[i]/total;
		}
	
	
		//generate a random seed to get the start index of motif that are used for the next round
		// generate a random number from 0 to 1
		Random random=new Random();
		double n= random.nextDouble();
		
		// get the index for that random propotion we randomly generated
		int index=0;
		for(int i=0; i<p.length;i++) {
			if (p[i]>n) {
				index =i;
				break;				
			}
		}
		
		return index;
	}

	
	
	
	
	
	
	
	
	
	
	// get the pssm score for the whole seq with a index list
	
	public static double total_pssm(int[]index,int m ,String[]seq2) {
		
		
		//new motif for each seq
		
		List<String> motif = new ArrayList<String>();
		  for(int i=0; i< seq2.length;i++) {				
			  motif.add(get_motif2(index[i],m,seq2[i]));  			
		}
		
		  // get the new index list from update
		  
		  String [] motif2= new String [seq2.length-1];
		  double pssm_score[]=new double[seq2.length];
			for(int i=0;i<index.length;i++) {
				
			  List<String> motif_tem = new ArrayList<String>();				
				  for(int n=0; n< motif.size();n++) {
						if(n!=i) {
							motif_tem.add(motif.get(n));
						}
				  }				 
				// exclude the sequence you are looking at to generate the motif
				//convert it to string array
					motif2= motif_tem.toArray(motif2);
					// generate a pssm matrix for this motif
					//public static double[][] pssm  (int m, String[] motif, String seq)					
					// when you get the pssm you need to exclude the current seq motif
					double pssm[][]=pssm(m,motif2,seq2[i]);
			
					//get the pssm sum for each seq  (motif is motif2, seq=motif[i])
					
					String seq=motif.get(i);
					
					double temp_sum=0;
					for(int j=0;j<seq.length();j++) {
						if(seq.charAt(j)=='A') {
							temp_sum+=pssm[j][0];
						}
						if(seq.charAt(j)=='C') {
							temp_sum+=pssm[j][1];
						}
						if(seq.charAt(j)=='G') {
							temp_sum+=pssm[j][2];
						}
						if(seq.charAt(j)=='T') {
							temp_sum+=pssm[j][3];
						}		
					}
					
					
					pssm_score[i]=temp_sum;
					
				
			}
			// get the total of pssm score
			double total_pssm=0;
			for (int i=0;i<pssm_score.length;i++) {
				total_pssm+=pssm_score[i];
			}

		return total_pssm;
	}
		
	
	
	
	
	
	
	
	
	
	
	
}
