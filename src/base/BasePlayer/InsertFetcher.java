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
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
//import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

@SuppressWarnings("deprecation")
public class InsertFetcher {

	byte[] seqresult;
	private StringBuffer seqBuffer;
	private RandomAccessFile chromo;
	private String resString;
	int seqStart=0;
	String seq;
	private String readlen;
	static Hashtable<String, Long[]> chromIndex = new Hashtable<String, Long[]>();
	int[][] consensus5;
	char[] result5 = {};
	int[][] consensus3;
	char[] result3 = {};
	
	public static void main(String[] args) {
		if(args.length < 2) {
			System.out.println("Give bam-file and reference.fa");
			System.exit(0);
		}
		
		try {
			//File refFile = new File("C:/HY-Data/RKATAINE/Rikurator/Rikurator/annotation/reference/default/hs37d5.fa");
			File refFile = new File(args[1]);
			RandomAccessFile referenceFile = new RandomAccessFile(refFile, "r");
			
			File chromindex = null;
						
			chromindex = new File(refFile.getAbsolutePath() +".fai");		
			BufferedReader reader = new BufferedReader(new FileReader(chromindex));
			String line;
			String[] split;
			 
			while((line = reader.readLine()) != null) {
				split = line.split("\t");
				
				Long[] add = {Long.parseLong(split[2]), Long.parseLong(split[1]), Long.parseLong(split[3])};
				chromIndex.put(split[0], add);
			}
			reader.close();
			
			SamReader samFileReader  = SamReaderFactory.make().open(new File(args[0])); //= new SAMFileReader(new File(args[0]));
		//	SAMFileReader samFileReader = new SAMFileReader(new File("V:/cg8/projects/CRC/c9_5119_N_LP6005134-DNA_A01/wgspipeline/align/c9_5119_N_LP6005134-DNA_A01_1.final.bam"));
			InsertFetcher fetcher = new InsertFetcher();
			
			fetcher.searchInsSites2(samFileReader, referenceFile);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	void searchInsSites(SamReader samFileReader, RandomAccessFile reference) {
		
		Hashtable<Integer, Character> bases = new Hashtable<Integer,Character>();
		bases.put(0, 'A');
		bases.put(1, 'C');
		bases.put(2, 'G');
		bases.put(3, 'T');
		bases.put(4, 'N');
		
		Hashtable<Character, Integer> bases2 = new Hashtable<Character, Integer>();
		bases2.put('A', 0);
		bases2.put('C', 1);
		bases2.put('G', 2);
		bases2.put('T', 3);
		
	
			java.util.List<Integer> refReadQualities = Collections.synchronizedList(new ArrayList<Integer>());
						
			SAMRecord cur;
			int readLengths = 0, clipLengths = 0, consensusStart = 0;
			String sequence = "", repeat = "";
			int maxbase = 0, allbases=0, disbases=0, disbasetemp = 0;
			List<Integer> maxindex = new ArrayList<Integer>();
			int del = 0;
		
			String readString = "";
			Hashtable<Integer, Integer> startVote = new Hashtable<Integer, Integer>(), endVote = new Hashtable<Integer, Integer>();
			int readPos =0, miscounter=0, clusterStart = 0, startReadMisSum =0, endReadMisSum = 0, clusterEnd = 0, max = 0, sum = 0, refReads = 0, altReads5 = 0, altReads3 = 0, breakPointStart = 0, breakPointEnd =0, support5=0, support3=0;
			java.util.List<Integer> positionList = Collections.synchronizedList(new ArrayList<Integer>());
			java.util.List<Object[]> startReadList = Collections.synchronizedList(new ArrayList<Object[]>());
			java.util.List<Object[]> endReadList = Collections.synchronizedList(new ArrayList<Object[]>());
			String rate3 = "", rate5 = "";
			Cigar cigar;
			Iterator<SAMRecord> iter = samFileReader.iterator();
			List<CigarElement> elements;
			
			while (iter.hasNext()) {	
				cur = iter.next();
				if(cur.getReadUnmappedFlag() || cur.getMappingQuality() < 10) {
					continue;
				}
				
				try {
					
					
					
					miscounter=0;
					clusterStart = 0;
					clusterEnd = 0;
					startVote.clear();
					endVote.clear();
					refReads = 0;					
					refReadQualities.clear();
					startReadList.clear();
					endReadList.clear();
					cigar = cur.getCigar();
					elements = cigar.getCigarElements();
					readString = cur.getReadString();
					
					if((cigar.toString().contains("S"))) {
							
							breakPointStart = breakPointEnd = 0;
							clusterStart = cur.getUnclippedStart();
							clusterEnd = cur.getUnclippedStart() + 150;
							
						if(cigar.toString().endsWith("S")) {
							clipLengths = elements.get(elements.size()-1).getLength();
							readLengths = cur.getReadLength();
							
							Object[] read = {(cur.getUnclippedStart()+cur.getReadLength())-elements.get(elements.size()-1).getLength()-1, readString.substring(readString.length()-elements.get(elements.size()-1).getLength(), readString.length()) };
							endReadList.add(read);
							
						//	endVote.put((cur.getUnclippedStart()+cur.getReadLength())-cur.getCigar().getCigarElement(cur.getCigar().getCigarElements().size()-1).getLength(), 1);
						}
						else {
							clipLengths = elements.get(0).getLength();
							readLengths = cur.getReadLength();
							
							Object[] read = {cur.getUnclippedStart(), readString.substring(0, elements.get(0).getLength()) };
							startReadList.add(read);
					//		startVote.put(cur.getUnclippedStart()+cur.getCigar().getCigarElement(0).getLength(), 1);
						}
						if(iter.hasNext()) {
							cur = iter.next();							
							cigar = cur.getCigar();
							elements = cigar.getCigarElements();
							readString = cur.getReadString();
						}
						else {
							break;
						}
				//	while(cur.getUnclippedStart() < clusterEnd) {
						
					while(clipLengths/(double)readLengths > 0.10 && cur.getUnclippedStart() < clusterEnd) {	
						try {
						if(cur.getCigarLength() > 0 && cur.getReadNegativeStrandFlag() && elements.get(0).getOperator().compareTo(CigarOperator.SOFT_CLIP)== 0) {
						
							if(cur.getMappingQuality() >= 10) {
								
								clipLengths += elements.get(0).getLength();
								readLengths += cur.getReadLength();
								Object[] read = {cur.getUnclippedStart(), readString.substring(0, elements.get(0).getLength()) };
								startReadList.add(read);	
								if(startReadList.size() > 500) {
									startReadList.clear();
									endReadList.clear();
									break;
								}
								
							}
						//	System.out.println("Start: " +(MethodLibrary.formatNumber(cur.getUnclippedStart()+cur.getCigar().getCigarElement(0).getLength())));
						}
						else if(!cur.getReadNegativeStrandFlag() && cigar.toString().endsWith("S")) {
							
							if(cur.getMappingQuality() >= 10) {
								clipLengths += elements.get(elements.size()-1).getLength();
								readLengths += cur.getReadLength();
								Object[] read = {(cur.getUnclippedStart()+cur.getReadLength())-elements.get(elements.size()-1).getLength()-1, readString.substring(readString.length()-elements.get(elements.size()-1).getLength(), readString.length()) };
								endReadList.add(read);
								if(endReadList.size() > 500) {
									startReadList.clear();
									endReadList.clear();
									break;
								}
								
							}
						//	System.out.println("End: " +(MethodLibrary.formatNumber((cur.getUnclippedStart()+cur.getReadLength())-cur.getCigar().getCigarElement(0).getLength())));
						}
						else {
							
							readLengths += cur.getReadLength();
							refReadQualities.add(cur.getMappingQuality());
							refReads++;
						}
						if(iter.hasNext()) {
							cur = iter.next();
							cigar = cur.getCigar();
							elements = cigar.getCigarElements();
							readString = cur.getReadString();
						}
						else {
							break;
						}
				//		System.out.println(clipLengths/(double)readLengths);
						}
						catch(Exception ex) {
							ex.printStackTrace();
						}
					}
					
					if(startReadList.size() < 3 && endReadList.size() < 3) {
						
						continue;
					}
					
					clusterEnd = cur.getUnclippedStart() + cur.getReadLength();
									
					Collections.sort(refReadQualities);
					if(refReadQualities.size() > 0 && refReadQualities.get(refReadQualities.size()/2) < 10) {
						
						continue;
					}
					startReadMisSum = 0;
					this.seqStart = clusterStart-100;
					sequence = getSeq(cur.getReferenceName(),clusterStart-100, clusterEnd+100, reference);
					
					for(int r = 0; r< startReadList.size(); r++) {
						miscounter = 0;
						readPos = (Integer)startReadList.get(r)[0];
						readString = (String)startReadList.get(r)[1];
						
						for(int c = 0; c< readString.length(); c++) {												
							
							if(readPos+c > this.seqStart && (readPos+c)-this.seqStart <sequence.length()) {
								if(sequence.charAt((readPos+c)-this.seqStart) != readString.charAt(c)) {
									miscounter++;
								}
							}					
						}
					
						if(readString.length() < 10 || (miscounter/(double)readString.length() > 0.05)) {
							
							startReadMisSum++;
							if(startVote.containsKey(readPos+readString.length())) {
								startVote.put(readPos+readString.length(), startVote.get(readPos+readString.length())+1);								
							}
							else {
								startVote.put(readPos+readString.length(), 1);								
							}		
						}
						else {
							startReadList.remove(r);
							r--;
						}
					}
					
					endReadMisSum = 0;
					for(int r = 0; r < endReadList.size(); r++) {
						miscounter = 0;
						readPos = (Integer)endReadList.get(r)[0];
						readString = (String)endReadList.get(r)[1];
						
						
						for(int c = 0; c< readString.length(); c++) {												
							
							if(readPos+c > this.seqStart && (readPos+c)-this.seqStart < sequence.length()) {
								if(sequence.charAt((readPos+c)-this.seqStart+1) != readString.charAt(c)) {
									miscounter++;
								}
							}					
						}
						if(readString.length() < 10 || miscounter/(double)readString.length() > 0.05) {
							
							endReadMisSum++;
							if(endVote.containsKey(readPos+1)) {
								endVote.put(readPos+1, endVote.get(readPos+1)+1);
								
							}
							else {								
								endVote.put(readPos+1, 1);
							}
						}
						else {
							endReadList.remove(r);
							r--;
						}
						
					}
					
					if(endVote.size() < 3 && startVote.size() < 3) {
						continue;
					}
					
					Iterator<Map.Entry<Integer, Integer>> iterator = startVote.entrySet().iterator();
				    positionList.clear();			
					max = altReads5 = altReads3 = support5 = support3 = 0;
					
					while (iterator.hasNext()) {
						  Map.Entry<Integer, Integer> entry = iterator.next();						 
						  altReads5 += entry.getValue();
						  
						  if(entry.getValue() > max) {
							  max = entry.getValue();
							  
							  positionList.clear();
							  positionList.add(entry.getKey());
							 
						  }
						  else if(entry.getValue() == max) {
							  positionList.add(entry.getKey());
							 
						  }						  
					}
					
					support5 = max;
					sum = 0;
					for(int j = 0; j<positionList.size(); j++) {
						sum+=positionList.get(j);
						
					}
					breakPointStart = (int)(sum/(double)positionList.size()+0.5);
					
					positionList.clear();
					max = 0;
					iterator = endVote.entrySet().iterator();
					
					while (iterator.hasNext()) {
						  Map.Entry<Integer, Integer> entry = iterator.next();
						  altReads3 += entry.getValue();
						  
						  if(entry.getValue() > max) {
							  max = entry.getValue();
							  positionList.clear();
							  positionList.add(entry.getKey());
							
						  }
						  else if(entry.getValue() == max) {
							  positionList.add(entry.getKey());
						  }			
					}
					
					support3 = max;
					sum = 0;
					for(int j = 0; j<positionList.size(); j++) {
						sum+=positionList.get(j);
						
					}
					breakPointEnd = (int)(sum/(double)positionList.size()+0.5);					
					
					for(int r = 0; r < endReadList.size(); r++) {
						seq = (String)endReadList.get(r)[1];
						readPos = (Integer)endReadList.get(r)[0];
						if(readPos+seq.length()< breakPointEnd || readPos > breakPointEnd+50) {
							endReadList.remove(r);
							r--;
						}
					}
						
					if(support5 < 3 && support3 < 3) {
						continue;
					}
					
					if(startReadMisSum < 3 && endReadMisSum < 3) {
						continue;
					}
					
					if(support5 < 3 && startReadMisSum/(double)endReadMisSum < 0.3) {
						breakPointStart = -1;
					}
					if(support3 < 3 && endReadMisSum/(double)startReadMisSum < 0.3) {
						breakPointEnd = -1;
					}
					
					if(support5 == 1) {
						breakPointStart = -1;
					}
					if(support3 == 1) {
						breakPointEnd = -1;
					}
					if(breakPointStart <= 0) {
						breakPointStart = -1;
					}
					if(breakPointEnd <= 0) {
						breakPointEnd = -1;
					}
					
				//	System.out.println(breakPointStart +" " +breakPointEnd);
					
					if(breakPointStart > -1 && startReadList.size() > 0) {
						Collections.sort(startReadList, new MethodLibrary.ReadListSorter());
						consensusStart = (Integer)startReadList.get(0)[0];
						consensus5 = new int[4][breakPointStart-consensusStart+1];
						result5 = new char[consensus5[0].length];
						disbases = 0;
						
						
						for(int c=0; c< consensus5.length; c++) {
							for(int h=0;h<consensus5[c].length; h++) {
								consensus5[c][h] = 0;
							}
						}
						
						
						for(int r = 0; r< startReadList.size(); r++) {
							readPos = (Integer)startReadList.get(r)[0];
							seq = (String)startReadList.get(r)[1];
							del = 0;
							for(int s = 0; s<seq.length(); s++) {
								
								if(bases2.containsKey(seq.charAt(s))) {
									if(readPos+s-del-consensusStart < consensus5[0].length && readPos+s-del-consensusStart > -1) {									
										consensus5[bases2.get(seq.charAt(s))][readPos+s-del-consensusStart]++;										
									}
									else {
										break;
									}
								}
								else if(seq.charAt(s) == 'D') {
									
									del++;
								}
							
							}
						}
						for(int c=0; c < consensus5[0].length; c++) {
							disbasetemp = 0;
							maxbase = 0;
							maxindex.clear();
							for(int h=0;h<consensus5.length; h++) {
								allbases += consensus5[h][c];
								disbasetemp += consensus5[h][c];
								
								if(consensus5[h][c] > maxbase) {
									maxbase = consensus5[h][c];
									maxindex.clear();
									maxindex.add(h);
								}
								else if(consensus5[h][c] == maxbase) {
									maxindex.add(h);
								}
								
							
							}
							disbases += disbasetemp - maxbase;
							if(maxbase == 0) {
								result5[c] = 'N';
							}
							else {
								result5[c] = bases.get(maxindex.get((int)(Math.random()*maxindex.size())));
							//	System.out.print(bases.get(maxindex.get((int)(Math.random()*maxindex.size()))));
							}
							
						//	System.out.println();
						}
						//rate5 = ""+MethodLibrary.round(disbases/(double)allbases, 2);
						if(allbases > 0) {
							rate5 = ""+MethodLibrary.round(disbases/(double)allbases, 2);
						}
						else {
							rate5 = "0";
						}
					}
					
					if(breakPointEnd > -1 && endReadList.size() > 0) {
						Collections.sort(endReadList, new MethodLibrary.ReadListSorter());
						readlen = (String)endReadList.get(endReadList.size()-1)[1];
						consensus3 = new int[4][(Integer)endReadList.get(endReadList.size()-1)[0]+readlen.length()-breakPointEnd+1];
						result3 = new char[consensus3[0].length];
						
					//	System.out.println("Forward:");
						
					//	consensusStart = (Integer)endReadList.get(0)[0];
						
						maxbase = allbases = disbases = disbasetemp = 0;					
						
						for(int r = 0; r< endReadList.size(); r++) {
							
							readPos = (Integer)endReadList.get(r)[0];
							seq = (String)endReadList.get(r)[1];
							del = 0;
						
							for(int s = 0; s<seq.length(); s++) {
								
								if(readPos+s < breakPointEnd-1) {
									
									continue;
								}
								
								if(bases2.containsKey(seq.charAt(s))) {
									if(readPos+s-del-(breakPointEnd-1) < consensus3[0].length) {
										try {
											consensus3[bases2.get(seq.charAt(s))][readPos+s-del-(breakPointEnd-1)]++;
										}
										catch(Exception e) {
											
											e.printStackTrace();
										}
									}
									else {
										
										break;
									}
								}
								else if(seq.charAt(s) == 'D') {
									
									del++;
								}
							//	System.out.print(seq.charAt(s));
							}
							
						//	System.out.println();
						}
						for(int c=0; c < consensus3[0].length; c++) {
							disbasetemp = 0;
							maxbase = 0;
							maxindex.clear();
							for(int h=0;h<consensus3.length; h++) {
								allbases += consensus3[h][c];
								disbasetemp += consensus3[h][c];
								
								if(consensus3[h][c] > maxbase) {
									maxbase = consensus3[h][c];
									maxindex.clear();
									maxindex.add(h);
								}
								else if(consensus3[h][c] == maxbase) {
									maxindex.add(h);
								}
						
							
							
							}
							disbases += disbasetemp - maxbase;
							if(maxbase == 0) {
								result3[c] = 'N';
							}
							else {
								result3[c] = bases.get(maxindex.get((int)(Math.random()*maxindex.size())));
							//	System.out.print(bases.get(maxindex.get((int)(Math.random()*maxindex.size()))));
							}
							
						//	System.out.println();
						}						
						if(allbases > 0) {
							rate3 = ""+MethodLibrary.round(disbases/(double)allbases, 2);
						}
						else {
							rate3 = "0";
						}
						
					}
					
					if(breakPointStart > -1 && breakPointEnd > -1) {
						if(breakPointStart <= breakPointEnd) {
							repeat = sequence.substring(breakPointStart-clusterStart, breakPointEnd-clusterStart);
							System.out.println(cur.getReferenceName() +"\t" +breakPointStart +"\t" +breakPointEnd +"\t" +(breakPointEnd-breakPointStart) +"\t-\t+\t" +repeat +"\t" +(new String(result5)) +"\t" +rate5 +"\t" +(new String(result3)) +"\t" +rate3 +"\t" +(altReads5+altReads3) +"\t" +refReads +"\t" +support5 +"\t" +support3 );
						}
						else {
							repeat = sequence.substring(breakPointEnd-clusterStart, breakPointStart-clusterStart);
							System.out.println(cur.getReferenceName() +"\t" +breakPointEnd +"\t" +breakPointStart +"\t" +(breakPointStart-breakPointEnd) +"\t+\t-\t" +repeat +"\t" +(new String(result5)) +"\t" +rate5 +"\t" +(new String(result3)) +"\t" +rate3 +"\t" +(altReads5+altReads3) +"\t" +refReads +"\t" +support3 +"\t" +support5 );
							
						}
					}
					else {
						if(breakPointEnd == -1) {
							if(sequence.length() < breakPointStart-clusterStart + 20) {
								repeat = sequence.substring(breakPointStart-clusterStart);
							}
							else {
								repeat = sequence.substring(breakPointStart-clusterStart, breakPointStart-clusterStart + 20);
							}
							System.out.println(cur.getReferenceName() +"\t" +breakPointStart +"\t" +(breakPointStart+1) +"\t0\t-\t-\t" +repeat +"\t" +(new String(result5)) +"\t" +rate5 +"\t-\t-\t" +(altReads5+altReads3) +"\t" +refReads +"\t" +support5 +"\t" +support3 );
						}
						else {
							if(breakPointStart-clusterStart - 20 < 0) {
								repeat = sequence.substring(0, breakPointEnd-clusterStart);
							}
							else {
								repeat = sequence.substring(breakPointStart-clusterStart, breakPointStart-clusterStart + 20);
							}
						//	repeat = sequence.substring(breakPointEnd-clusterStart -20, breakPointEnd-clusterStart);
							System.out.println(cur.getReferenceName() +"\t" +breakPointEnd   +"\t" +(breakPointEnd+1)   +"\t0\t+\t+\t" +repeat +"\t-\t-\t" +(new String(result3)) +"\t" +rate3 +"\t" +(altReads5+altReads3) +"\t" +refReads +"\t" +support5 +"\t" +support3 );
							
						}
					}
					
				}
				}
				catch(Exception e) {
					e.printStackTrace();
				}
			}
			
		
	}
	
void searchInsSites2(SamReader samFileReader, RandomAccessFile reference) {		
									
			SAMRecord cur;			
			String sequence = "";			
			String readString = "";
			int interval = 1000;
			int seqStart = 0, seqEnd = interval;
			
			Iterator<SAMRecord> iter = samFileReader.iterator();
			List<CigarElement> elements;
			int count = 0, miscount = 0;
			int readpointer = 0, seqpointer = 0;
			while (iter.hasNext()) {						
				cur = iter.next();
				
				if(cur.getReadUnmappedFlag() || cur.getMappingQuality() < 10) {
					continue;
				}
				
				count++;
				if(count == 50) {
					break;
				}
				if (cur.getUnclippedEnd() > seqEnd) {			
					
					seqStart = cur.getUnclippedStart()-100;
					seqEnd = seqStart+interval;
					sequence = getSeq(cur.getReferenceName(), seqStart, seqEnd, reference); 
					
				}			
				
				try {
					miscount = 0;
					seqpointer = cur.getUnclippedStart()-seqStart;
					
					for(int i = 0; i<cur.getUnclippedStart()-seqStart; i++) {
						System.out.print(" ");
					}
					
					if(cur.getCigarLength() == 1) {
						readString = cur.getReadString();
						for(int i = 0; i<cur.getReadString().length(); i++) {
							if(readString.charAt(i) != sequence.charAt(seqpointer+i)) {
								miscount++;
							}
							System.out.print(readString.charAt(i));
						}
						
					}
					else  {
						
						readpointer = 0;
						elements = cur.getCigar().getCigarElements();
						readString = cur.getReadString();
						
						for(int i = 0; i<elements.size(); i++) {
							
							if(elements.get(i).getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
								for(int j = 0; j<elements.get(i).getLength(); j++) {
									if(readString.charAt(readpointer) != sequence.charAt(seqpointer+j)) {
										miscount++;										
									}
									
									System.out.print(readString.charAt(readpointer));	
									
									readpointer++;
								}
								seqpointer += elements.get(i).getLength();
								continue;
							}
							
							if(elements.get(i).getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH)== 0) {
								
								for(int j = 0; j<elements.get(i).getLength(); j++) {
									if(readString.charAt(readpointer) != sequence.charAt(seqpointer+j)) {
										miscount++;										
									}
									System.out.print(readString.charAt(readpointer));	
								
									readpointer++;
								}
								seqpointer += elements.get(i).getLength();
								continue;
							}
							if(elements.get(i).getOperator().compareTo(CigarOperator.INSERTION)== 0) {
								readpointer += elements.get(i).getLength();
							//	seqpointer++;
								continue;
														
								
							}
							if(elements.get(i).getOperator().compareTo(CigarOperator.DELETION)== 0) {
								
								for(int j = 0; j<elements.get(i).getLength(); j++) {
									System.out.print(" ");	
									
									continue;
								}				
								seqpointer+=elements.get(i).getLength();
								
							}
						}
					}
					
					System.out.println(" " +miscount +" " +cur.getCigarString());
					
					
				}
				catch(Exception e) {
					e.printStackTrace();
				}
			}
			
		
	}

	public String getSeq(String chrom, int start, int end, RandomAccessFile seqchrom) {
		
		try {
			seqresult = new byte[(end-start)+((end-start)/60)+1];
			
		}
		catch(Exception e) {
			e.printStackTrace();
			return "";
		}
		
		seqBuffer = new StringBuffer();
		chromo = seqchrom;
		
		try {
		
		
			try {
				chromo.seek((chromIndex.get(chrom)[0]+(start)+((start)/60))-1);
			}
			catch(Exception e) {
				
				e.printStackTrace();
			}
			chromo.readFully(seqresult);
			
			if(seqresult[0] == 10) {
				chromo.seek((chromIndex.get(chrom)[0]+(start-1)+((start)/60))-1);
				seqchrom.readFully(seqresult);
				
			}			
					
			for(int i= 0; i< seqresult.length; i++){
				
				if(seqresult[i] != 10) {
					seqBuffer.append((char)seqresult[i]);
			//		resString+=(char)seqresult[i];					
				}
				
				resString = seqBuffer.toString().toUpperCase();
				if(resString.length() >= end-start) {
					break;
				}			
			}
			
		}
		catch (Exception e) {
			e.printStackTrace();
		
		}	

		
		
		return resString;			
	}
}
