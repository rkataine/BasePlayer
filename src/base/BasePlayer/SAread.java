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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;

public class SAread {

	ReadNode read;
	int readlength = 0;
	ArrayList<Object[]> reads = new ArrayList<Object[]>();
	static int chrom = 0, pos = 1, forward = 2, relativepos = 3, relativelength = 4, SARead = 5; 
	
	public SAread(ReadNode read) {
		try {
		this.read = read;		
		String[] SAs = read.SA.split(";");
		String[] sa;
		String chr = read.split.chrom;
		htsjdk.samtools.Cigar cigar;
		this.readlength = read.getCigar().getReadLength();
		cigar = read.getCigar();
		double readpos= 0, relativelength;
		int pos = read.getPosition(), length = this.readlength;
		if(!read.getPrimary()) {
			length = read.getCigar().getReadLength();
			if(read.isForward()) {
				if(cigar.getCigarElement(0).getOperator().compareTo(CigarOperator.HARD_CLIP) == 0) {
					this.readlength+=cigar.getCigarElement(0).getLength();
					
				}
				if(cigar.getCigarElement(cigar.getCigarElements().size()-1).getOperator().compareTo(CigarOperator.HARD_CLIP) == 0) {
					this.readlength+=cigar.getCigarElement(cigar.getCigarElements().size()-1).getLength();
				}
				readpos = cigar.getCigarElement(0).getLength()/(double)this.readlength;
			}
			else {
				if(cigar.getCigarElement(0).getOperator().compareTo(CigarOperator.HARD_CLIP) == 0) {
					this.readlength+=cigar.getCigarElement(0).getLength();
					
				}
				if(cigar.getCigarElement(cigar.getCigarElements().size()-1).getOperator().compareTo(CigarOperator.HARD_CLIP) == 0) {
					this.readlength+=cigar.getCigarElement(cigar.getCigarElements().size()-1).getLength();
					readpos = cigar.getCigarElement(cigar.getCigarElements().size()-1).getLength()/(double)this.readlength;
				}
			}
		}
		else {		
			
			if(read.isForward()) {
				if(cigar.getCigarElement(0).getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
					length-=cigar.getCigarElement(0).getLength();
					readpos = cigar.getCigarElement(0).getLength()/(double)this.readlength;
				}
				if(cigar.getCigarElement(cigar.getCigarElements().size()-1).getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
					length-=cigar.getCigarElement(cigar.getCigarElements().size()-1).getLength();
				}
			}
			else {
				if(cigar.getCigarElement(0).getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
					length-=cigar.getCigarElement(0).getLength();
				}
				if(cigar.getCigarElement(cigar.getCigarElements().size()-1).getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
					length-=cigar.getCigarElement(cigar.getCigarElements().size()-1).getLength();
					readpos = cigar.getCigarElement(cigar.getCigarElements().size()-1).getLength()/(double)this.readlength;
				}	
			}
		}
		relativelength = length/(double)this.readlength;	
		
		reads.add(makeReadList(chr, pos, read.isForward(), readpos, relativelength, read));
		
		boolean forward = true;
		
		for(int i = 0; i<SAs.length; i++) {
			sa = SAs[i].split(",");
			chr = sa[0];			
			pos = Integer.parseInt(sa[1]);	
			cigar = TextCigarCodec.decode(sa[3]);
			length = cigar.getReadLength();
			readpos = 0;
			if(sa[2].equals("-")) {
				forward = false;				
				if(cigar.getCigarElement(0).getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
					length-=cigar.getCigarElement(0).getLength();
				}
				if(cigar.getCigarElement(cigar.getCigarElements().size()-1).getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
					length-=cigar.getCigarElement(cigar.getCigarElements().size()-1).getLength();
					readpos = cigar.getCigarElement(cigar.getCigarElements().size()-1).getLength()/(double)this.readlength;
				}				
			}
			else {
				forward = true;
				if(cigar.getCigarElement(0).getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
					length-=cigar.getCigarElement(0).getLength();
					readpos = cigar.getCigarElement(0).getLength()/(double)this.readlength;
				}
				if(cigar.getCigarElement(cigar.getCigarElements().size()-1).getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {
					length-=cigar.getCigarElement(cigar.getCigarElements().size()-1).getLength();
				}
			}
			relativelength = length/(double)this.readlength;			
			reads.add(makeReadList(chr, pos, forward, readpos, relativelength, null));
			
		}
		readSorter sorter = new readSorter();
		
		Collections.sort(reads, sorter);
		if(read.getMates() != null) {
			for(int i = 0 ; i<read.getMates().size();i++) {				
				
					for(int r = 0 ; r<reads.size(); r++) {
					
						if((int)reads.get(r)[SAread.pos] == read.getMates().get(i).getPosition()) {	
							
							reads.get(r)[SAread.SARead] = read.getMates().get(i);
							
						}
						
					}
			}
		}
						
		}
		catch(Exception e) {
			e.printStackTrace();
		}
}
public static class readSorter implements Comparator<Object[]> {
		
		public int compare(Object[] o1, Object[] o2) {  
		
	            
	        if ((double)o1[SAread.relativepos] <  (double)o2[SAread.relativepos]) {  
	                return -1;  
	        } 
	        else if((double)o1[SAread.relativepos] >  (double)o2[SAread.relativepos]) {  
	        		return 1;  
	        }
	        else {
	        	return 0;
	        } 	       
	}
}

Object[] makeReadList(String chr, int position, boolean forward, double startpos, double length, ReadNode read) {
	Object[] readArray = new Object[6];
	readArray[SAread.chrom] = chr;
	readArray[SAread.pos] = position;
	readArray[SAread.forward] = forward;
	readArray[SAread.relativepos] = startpos;
	readArray[SAread.relativelength] = length;
	readArray[SAread.SARead] = read;
	return readArray;
}
}
