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

import htsjdk.samtools.util.BlockCompressedInputStream;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

public class Sample implements Serializable{
	
	
	private static final long serialVersionUID = 1L;
	private transient HashMap<SplitClass, Reads> readHash;
	Short maxCoverage = 0;
	Short prepixel = 0;
	boolean reading = false, multiVCF = false, multipart = false;
	private Short index;
	private String sampleName;
	private String tabixfile;
	private transient BlockCompressedInputStream inputStream;
	
	short infolocation = 0;
	int varcount = 0, homozygotes = 0, heterozygotes = 0, indels = 0, snvs=0, ins=0, coding = 0;
	File samFile;
	long vcfEndPos = 0;
	Integer longestRead = 0;
	String chr = "", vcfchr = "";
	private transient HashMap<String, ArrayList<ReadNode>> mates;
	double callrates = 0.0;
	double[] mutationTypes;
	short phred = 33;
	public boolean MD = false;
	public int sitioRate = 0, versioRate = 0;
	
	Sample(String sampleName, short index, String tabixfile) {		
		if(tabixfile != null) {
			this.tabixfile = tabixfile;		
			try {
				inputStream = new BlockCompressedInputStream(new File(this.tabixfile));
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		this.sampleName = sampleName;
		this.index = index;
	
	}
	
	void resetreadHash() {
		if(readHash == null) {
			readHash = new HashMap<SplitClass, Reads>();
		}
		for(int i = 0; i<Main.drawCanvas.splits.size(); i++) {
			Reads reads = new Reads();
			reads.sample = this;
			readHash.put(Main.drawCanvas.splits.get(i), reads);
		}
		

	}
/*	ArrayList<Reads> getreadHash() {
		return this.readHash;
	}*/
	BlockCompressedInputStream getVCFInput() {
		return this.inputStream;
	}
	void setInputStream() {
		try {
		 this.inputStream = new BlockCompressedInputStream(new File(this.tabixfile));
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	HashMap<SplitClass, Reads> getreadHash() {
		return this.readHash;
	}
	HashMap<String, ArrayList<ReadNode>> getMates() {
		return this.mates;
	}
	
	void setMates() {
		this.mates = new HashMap<String, ArrayList<ReadNode>>();
	}
	
	short getIndex() {
		return this.index;
	}
	void setIndex(short index) {
		this.index = index;
	}	
	String getName() {
		return this.sampleName;
	}
	String getTabixFile() {
		return this.tabixfile;
	}
		
}
