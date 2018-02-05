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
import java.io.BufferedReader;
import java.io.BufferedWriter;

import java.io.FileInputStream;

import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;


public class ConvertExons {
	public static void main(String[] args) {
		try {
			
			GZIPInputStream gzip = new GZIPInputStream(new FileInputStream("X:/cg8/Ensembl-update/fetch/exons37.temp.gz"));
			BufferedWriter writer = new BufferedWriter(new FileWriter("X:/cg8/Riku/exons.txt"));
			BufferedReader reader = new BufferedReader(new InputStreamReader(gzip));
			String line, transcript = "";
			String header = "#Chrom\tGeneStart\tGeneEnd\tName\tExonCount\tStrand\tENSG\tENST\tUniProt\tCanonical\tBiotype\tCodingStart\tCodingEnd\tExonStarts\tExonEnds\tStartPhases\tDescription\n";
			String[] split;
			int exoncount = 0;
			boolean first = true, started = false;
			line = reader.readLine();
			split = line.split("\t");
			String chrom = "", phase;
			String genestart, geneend, strand, uniprot = "-", canonical, biotype, codingstart, codingend, description, name, ensg;
			StringBuffer exonstarts, exonends, phases;
			writer.write(header);
			while(true) {
				
				if(!first) {
					line = reader.readLine();
					if(line == null) {
						break;
					}
				}
				
				
				if(first) {
					
					transcript = split[10];
					first = false;
					name = split[1];
					if(split[8].matches("23")) {
						chrom = "X";
					}
					else if(split[8].matches("24")) {
						chrom = "Y";
					}
					else if(split[8].matches("25")) {
						chrom = "MT";
					}
					else {
						chrom = split[8];
					}
					
				    exoncount = 0;
				    description = split[5];
				    ensg = split[0];
				    codingstart = "";
				    codingend = "";
				    
				    if(split[9].equals("1")) {
				    	 strand = "+";
				    }
				    else {
				    	strand = "-";
				    }
				   
				    biotype = split[2];
				    
				    if(split[4].equals(transcript)) {
				    	canonical ="1";
				    }
				    else {
				    	canonical = "-";
				    }
				    started = false;
				    exonstarts  = new StringBuffer();
				    exonends    = new StringBuffer();
				    phases  	= new StringBuffer();
				    
					while(transcript.equals(split[10]) ) {
						exonstarts.append("" +(Integer.parseInt(split[12])-1) +",");
						exonends.append(""+(Integer.parseInt(split[13])) +",");
						if(strand.equals("+")) {
							if(split.length < 17) {
								
								phases.append(split[14] +",");
								
							}
							else {
								if(!started) {
									phases.append("0,");
									codingstart = ""+(Integer.parseInt(split[12])+Integer.parseInt(split[16]) -2);
									codingend = ""+(Integer.parseInt(split[12])+Integer.parseInt(split[17].replaceAll(" ", "")) -1);
									started = true;
								}
								else if(started) {
									phases.append(split[14] +",");
								
									codingend = ""+(Integer.parseInt(split[12])+Integer.parseInt(split[17].replaceAll(" ", "")) -1);
								}
							}
						}
						else {
							phase = split[14];
							if(split.length < 17) {										
								phases.append(phase +",");							
							}
							else {
								if(!started ) {
									if( phase.equals("-1")){
										phases.append("0,");
									}
									else {
										phases.append(phase+",");
									}
									codingstart = ""+(Integer.parseInt(split[13])-Integer.parseInt(split[16]));
									started = true;
									codingend = ""+(Integer.parseInt(split[13])-Integer.parseInt(split[17].replaceAll(" ", ""))+1);
								}
								else if(started) {
									phases.append("0,");
									codingend = ""+(Integer.parseInt(split[13])-Integer.parseInt(split[17].replaceAll(" ", ""))+1);									
								}								
							}
						}
						exoncount++;						
						
						line = reader.readLine();
						if(line == null) {
							break;
						}
						split = line.split("\t");
					}
					first = true;
					genestart = exonstarts.toString().split(",")[0];
					geneend = exonends.toString().split(",")[exonends.toString().split(",").length-1];
					if(codingstart.equals("")) {
						codingstart = geneend;
						codingend = geneend;
					
					}
			//		System.out.print(chrom +"\t" +genestart +"\t" +geneend +"\t" +name +"\t" +exoncount +"\t" +strand +"\t" +ensg +"\t" +transcript +"\t" +uniprot +"\t" +canonical +"\t" +biotype +"\t" +codingstart +"\t" +codingend +"\t" +exonstarts +"\t" +exonends +"\t" +phases +"\t" +description +"\n");
					writer.write(chrom +"\t" +genestart +"\t" +geneend +"\t" +name +"\t" +exoncount +"\t" +strand +"\t" +ensg +"\t" +transcript +"\t" +uniprot +"\t" +canonical +"\t" +biotype +"\t" +codingstart +"\t" +codingend +"\t" +exonstarts +"\t" +exonends +"\t" +phases +"\t" +description +"\n");
					
				}
				if(line == null) {
					break;
				}				
				
			}
			writer.close();
			reader.close();
			gzip.close();
		}
		catch(Exception e) {
			
			e.printStackTrace();
		}
	}
}
