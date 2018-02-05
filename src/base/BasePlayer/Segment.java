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

import java.util.ArrayList;
import java.util.Collections;

public class Segment implements Comparable<Segment> {
	int pos;
	String readname;
	String readString;
	//htsjdk.samtools.Cigar cigar;
	String cigarString;
	boolean inc;
	int freq;
	String mateRef;
	int matePos;
	boolean mateDir, check = false;
	int level = 0;
	int id;
	int chrom;
	boolean direction;
	int length;
	String mismatches;
	int insert;
	int end;
	int quality;
	Cigar cigar;
	double x = 0.0;
	static int counter = 0;
	int y=0;
	
	java.util.List<Segment> splits = Collections.synchronizedList(new ArrayList<Segment>());
	

	public Segment(int chrom, int pos, String readname, boolean direction, boolean inc, int length, String mismatches, int insert, Cigar cigar, String readstring, String mateRef, int matePos, boolean mateDir, int quality) {
		this.chrom = chrom;
		this.pos = pos;
		this.readname = readname;
		this.inc = inc;
		this.id = counter++;
		this.direction = direction;
		this.length = length;
		this.mismatches = mismatches;
		this.insert = insert;
		this.readString = readstring;
		this.cigar = cigar;
		this.mateRef = mateRef;
		this.matePos = matePos;
		this.mateDir = mateDir;
		this.quality = quality;
	}

	public Segment(int chrom, int pos, String readname, boolean direction, boolean inc, int length, int insert, int end, String readstring, String mateRef, int matePos, boolean mateDir, Cigar cigar, int quality) {
		this.chrom = chrom;
		this.pos = pos;
		this.readname = readname;
		this.inc = inc;
		this.id = counter++;
		this.direction = direction;
		this.length = length;
		this.cigar = cigar;
		this.insert = insert;
		this.readString = readstring;		
		this.mateRef = mateRef;
		this.matePos = matePos;
		this.mateDir = mateDir;
		this.end = end;
		this.quality = quality;
	}
	
	public int compareTo(Segment seg) {
		if (this.chrom < seg.chrom) return -1;
		if (this.chrom > seg.chrom) return 1;
		if (this.pos < seg.pos) return -1;
		if (this.pos > seg.pos) return 1;
	/*	if (this.read < seg.read) return -1;
		if (this.read > seg.read) return 1;
		*/
		if (this.id < seg.id) return -1;
		if (this.id > seg.id) return 1;
		return 0;
	}
	
}
