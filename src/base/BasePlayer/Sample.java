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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JOptionPane;

public class Sample implements Serializable{
	
	
	private static final long serialVersionUID = 1L;
	private transient HashMap<SplitClass, Reads> readHash;
	Short maxCoverage = 0;
	Short prepixel = 0;
	boolean reading = false, multiVCF = false, multipart = false, removed = false;
	private Short index;
	private String sampleName;
	private String tabixfile;
	private transient BlockCompressedInputStream inputStream;
	private transient BufferedReader vcfreader;
	short infolocation = 0;
	int varcount = 0, homozygotes = 0, heterozygotes = 0, indels = 0, snvs=0, ins=0, coding = 0;
	File samFile;
	long vcfEndPos = 0;
	Integer longestRead = 0;
	String chr = "", vcfchr = "";
	Short somaticColumn = null;
	private transient HashMap<String, ArrayList<ReadNode>> mates;
	double callrates = 0.0;
	double[] mutationTypes;
	short phred = 33;
	public boolean MD = false, CRAM = false, hasMates = false,complete = false;
	public int sitioRate = 0, versioRate = 0;
	
	Sample(String sampleName, short index, String tabixfile) {		
		if(tabixfile != null) {
			this.tabixfile = tabixfile;		
			try {
				setInputStream();
				//inputStream = new BlockCompressedInputStream(new File(this.tabixfile));
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		this.sampleName = sampleName;
		this.index = index;
	
	}
	
	void resetreadHash() {
		if(samFile != null) {
			if(samFile.getName().toLowerCase().endsWith(".cram")) {
				this.CRAM = true;
			}
		}
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
	BufferedReader getIDXInput() {
		return this.vcfreader;
	}
	void setInputStream() {
		try {
		if(this.tabixfile.endsWith(".gz")) {
			this.inputStream = new BlockCompressedInputStream(new File(this.tabixfile));
		}
		else {
			this.vcfreader= new BufferedReader(new FileReader(this.tabixfile));
		}
		}
		catch(Exception e) {
			
			this.inputStream =null;
			e.printStackTrace();
			JOptionPane.showMessageDialog(Main.chromDraw, e.getMessage() , "Error", JOptionPane.ERROR_MESSAGE);
			
		}
	}
	BufferedReader getVCFReader() {
		return this.vcfreader;
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
	short getMultiIndex() {
		if(this.multipart) {
			return (short)(this.index-1);
		}
		else {
			return this.index;
		}
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
