/* Author: Riku Katainen @ University of Helsinki
 * 
 * Tumor Genomics Group (http://research.med.helsinki.fi/gsb/aaltonen/) 	
 * Contact: riku.katainen@helsinki.fi / help@baseplayer.fi
 * 
 * LICENSE: 
 * 
 * GNU AFFERO GENERAL PUBLIC LICENSE
 * Version 3, 19 November 2007
 *  
 */
package base.BasePlayer;

import htsjdk.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JFrame;

public class Test {

	 	static HashMap<String, Integer> baseMap = new HashMap<String, Integer>();
	 	static String chrom = "";
		public static void main(String[] args) {
		String seq = "GGCGCCACCTTGTGGTT";
		String mutated ="GGTGCCACCTTGTGGTT";
		
		String pmf = "5423,2600,0,64,641,46,5823,87,0,37,10566,90,183,2,406,2660,4213;546,371,8052,0,9591,10107,5201,7091,8872,49,290,1616,5157,0,114,4936,3590;1733,11366,0,12480,8,0,207,0,1,23,1625,9828,0,8472,7463,143,2456;4208,632,505,86,46,0,241,3888,0,14151,1355,764,7952,121,2221,9766,1302";
			
		String[] split = pmf.split(";");
		String[] values = split[0].split(",");
		int[][] matrix = new int[4][values.length];
		
		for(int i = 1; i<4; i++) {
			for(int j = 0; j<values.length; j++) {
				matrix[i-1][j] = Integer.parseInt(values[j]);
			}
			values = split[i].split(",");
			
			
		}
		double sum;
		HashMap<Character, Integer> map = new HashMap<Character, Integer>();
		HashMap<Character, Double> background = new HashMap<Character, Double>();
		map.put('A', 0);
		map.put('C', 1);
		map.put('G', 2);
		map.put('T', 3);		
		background.put('A', 0.3);
		background.put('C', 0.2);
		background.put('G', 0.2);
		background.put('T', 0.3);		
		Double value, total = 0.0, mutTotal = 0.0;
		double pseudocount = 0.8;
		double mutatedvalue;
		
		for(int j = 0; j<matrix[0].length; j++) {
		try {
			sum = matrix[0][j] + matrix[1][j] +matrix[2][j] +matrix[3][j];
			
			value = matrix[map.get(seq.charAt(j))][j]/(double)sum;
			mutatedvalue = matrix[map.get(mutated.charAt(j))][j]/(double)sum;
			if(value == 0) {
				value = pseudocount;
				
			}
			if(mutatedvalue == 0) {
				mutatedvalue = pseudocount;
			}
			value = value*Math.log(value/background.get(seq.charAt(j)));
			mutatedvalue = mutatedvalue*Math.log(mutatedvalue/background.get(mutated.charAt(j)));
			total +=value;
			mutTotal += mutatedvalue;
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		}
		System.out.println(mutTotal-total);
			
			
		/*	if (args.length < 2) {
				System.out.println("Give file with break points e.g. 1:43243513 and bam file or folder including bam-files.");
				System.exit(0);
			}
			makeBasemap();
			try {
				BufferedReader file = new BufferedReader(new FileReader(args[0]));
				File bamfile = new File(args[1]);
				String line;
			
				while((line = file.readLine()) != null) {
					chrom = line.split(":")[0];
					if(bamfile.isDirectory()) {
						
						File[] bamfiles = bamfile.listFiles(new FilenameFilter() {
			  	    		 
		  	    	     public boolean accept(File dir, String name) {
		  	    	        return name.toLowerCase().matches(".*_" +chrom +"\\..*\\.bam");
		  	    	     }
		  	    	   	 });
						
						getConsSeq(line, bamfiles[0].getCanonicalPath());
					}
					else {
						getConsSeq(line, bamfile.getCanonicalPath());
					}
				}
				
				
			
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			*/
		}
		
		
	/*	public static void getConsSeq(String position, String bamfile ) {
			String[] pos = position.split(":");
			String chrom = pos[0];
			String breakpoint = pos[1];			
			File sample = new File(bamfile);
			
			ReadNode read;
		
			int minpos = Integer.MAX_VALUE, maxpos = 0;
			ArrayList<Object[]> readlist = new ArrayList<Object[]>();
			
			for(int i = 0; i<readsample.getreadHash().get(split).getReads().size();i++)  {
				
				read = readsample.getreadHash().get(split).getReads().get(i);
				do {
					if(read.getPosition() < centerpos && read.getPosition()+read.getWidth() > centerpos) {
						if(read.getMismatches() != null && read.getMismatches().size() > 10) {
							if(read.getPosition() < minpos) {
								minpos = read.getPosition();
							}
							if(read.getPosition()+read.getWidth() > maxpos) {
								maxpos = read.getPosition()+read.getWidth();
							}
							SAMRecord readsam = Main.fileReader.getRead(chromosomeDropdown.getSelectedItem().toString(),read.getPosition(), read.getPosition()+read.getWidth(),read.getName(), readsample.getreadHash().get(split));
							
							Object[] adder = {read.getPosition(), readsam.getReadString()};
							readlist.add(adder);
						}
					}
				
				}
				while ((read = read.getNext()) != null);				
			}
			
			int[][] matrix = new int[5][maxpos-minpos];
			
			for(int i = 0 ; i<5;i++) {
				for(int j = 0;j<maxpos-minpos;j++) {
					matrix[i][j] = 0;
				}
			}
			for(int i = 0 ; i<readlist.size();i++) {
				for(int j =0; j<readlist.get(i)[1].toString().length(); j++) {		
					
					matrix[baseMap.get(""+(readlist.get(i)[1].toString().charAt(j)))-1][(int)readlist.get(i)[0]-minpos+j]++;
				}
			}
			StringBuffer buffer = new StringBuffer("");//readsample.getName() +"\nPosition: " +chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(minpos) +"\nBreak point: " +chromosomeDropdown.getSelectedItem().toString() +":" +MethodLibrary.formatNumber(centerpos) +"\n" );
	
			read = null;
		
			buffer.append(">" +sample.getName() +"|BP="+position  +" (LeftPosition=" +chrom +":" +MethodLibrary.formatNumber(minpos) +")\n");
			int max=0,maxindex=0;
			String[] bases = {"A","C","G","T"};
			StringBuffer fasta = new StringBuffer("");
			for(int j =0; j<maxpos-minpos; j++) {
				max = 0;
				maxindex = 0;
				for(int i = 0; i<4; i++) {
					if(matrix[i][j] > max ) {
						max = matrix[i][j];
						maxindex = i;
					}
					else if(max > 0 && matrix[i][j] == max) {
						max = 0;
						maxindex = -1;
					}
				}
				if(maxindex > -1) {
					fasta.append(bases[maxindex]);
				}
				else {
					fasta.append("N");
				}
			}
			buffer.append(fasta.toString());
			System.out.println(buffer.toString());
	
		  
		}
		public static void makeBasemap() {
			baseMap.put("A", 1);
			baseMap.put("C", 2);
			baseMap.put("G", 3);
			baseMap.put("T", 4);
			baseMap.put("N", 5);
		}*/
}
