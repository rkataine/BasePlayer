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
import htsjdk.samtools.CRAMFileReader;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileReader;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
//import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.index.Block;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.util.BlockCompressedOutputStream;

import java.awt.Dimension;

import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;

import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.GZIPInputStream;
import javax.swing.JScrollPane;
import javax.swing.SwingWorker;

import org.apache.commons.io.FilenameUtils;

@SuppressWarnings("deprecation")
public class FileRead extends SwingWorker< String, Object >  {
	
	Boolean readVCF = false, changeChrom = false, readBAM = false,searchInsSites = false, varcalc = false, getreads = false;
	File[] files;
	String chrom = "0";
	int pos;
	static boolean searchingBams = false;
	int calls1, calls2;	
	static boolean nobeds = false;
	private boolean found, noref = false, multi = false;
	private int xpos;	
	private ReadNode mundane = null;
	
	static int searchwindow = 1000;
	private boolean left;
	boolean stop = false;
	private ReadNode lastAdded;	
	private ReadNode addNode;
	static VarNode lastVar = null, lastWriteVar = null, returnnode = null;
	static Gene currentGene = new Gene();
	static int currentGeneEnd = 0;
	static StringBuffer sampleString = new StringBuffer("");
	SamReader samFileReader;	
	CRAMFileReader CRAMReader = null;
	Iterator<SAMRecord> bamIterator = null;	
	SAMRecord samRecord, samRecordBuffer;
	private String[] info;
	private String[] coverages;
	private Sample sample;
	private Float quality, gq;
	static boolean novars = false;
	String basecontext = "";
	int start, end, viewLength;
	double pixel;
	public Reads readClass;	
	static final int headnode = 0, tailnode = 1;	
	public static VarNode head;
	VarNode current;
	static int firstReadPos;
	static int lastReadPos;
	private int startY;
	private boolean right;
	public SplitClass splitIndex;
	private Sample currentSample;
	public boolean statcalc;
	private boolean genotype;
	private int mutcount, linecounter = 0;
	private short firstallele;
	private short secondallele;
	private String altbase;
	private short refallele;
	private int refcalls;
	private short altallele;
	private int altcalls;
	private String altbase2;
	static String[] headersplit;
	private boolean first;
	private int samplecount;
	private int middle;
	private ReadBuffer currentread;
	private int searchPos;
	private boolean isDiscordant;
	private double[][] coverageArray;
	private int pointer;
	private int oldstart;
	private int addlength;
	private int readstart;
	private int readpos;
	private int mispos;	
	private CigarElement element;	
	boolean firstCov;	
	public boolean firstSample = false;
	private int gtindex;
	private int timecounter = 0;
	static boolean bigcalc=false;
	static boolean changing = false;
	static boolean readFiles;	
	public static boolean search;
	public static int searchStart;
	public static int searchEnd;
	public static boolean cancelvarcount= false;
	public static boolean cancelfileread = false;
	public static boolean cancelreadread= false;
	public static BufferedWriter output = null;
	public static BlockCompressedOutputStream outputgz =null;
	public static File outFile = null;
	public static TabixIndexCreator indexCreator;	
	public static int lastpos = 0;
	public static String outputName = "";
	public static long filepointer = 0;	
	HashMap<String, Integer[]> contexts;
	
	public FileRead(File[] files) {
		this.files = files;
	}
	public FileRead() {
		
	}
	
	protected String doInBackground() throws Exception {
		
		if(readVCF) {
			if(Main.drawCanvas.drawVariables.projectName.equals("Untitled")) {
				Main.frame.setTitle("BasePlayer - Untitled Project");
			}
			
			readVCF(files);
			readVCF = false;		
		
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
		}
		else if(readBAM) {
			if(Main.drawCanvas.drawVariables.projectName.equals("Untitled")) {
				Main.frame.setTitle("BasePlayer - Untitled Project");
			}
			Main.drawCanvas.loading("Loading samples");					
			readBAM(files);
			readBAM = false;
			Main.drawCanvas.ready("Loading samples");
		}
		
		else if(changeChrom) {
		/*	if(!Main.drawCanvas.loadingtext.equals("note")) {
				Main.drawCanvas.loading("Loading Chromosome " +chrom);					
			}*/
			changeChrom(chrom);	
			changeChrom = false;	
		/*	if(!Main.drawCanvas.loadingtext.equals("note")) {
				Main.drawCanvas.ready("Loading Chromosome " +chrom);
			}*/
			Draw.updatevars = true;
			Draw.updateReads = true;
			Main.drawCanvas.repaint();
		}		
		else if(varcalc) {
			Main.drawCanvas.loading("Processing variants...");		
			
			try {
			if(FileRead.head.getNext() == null) {
				
				if(VariantHandler.allChroms.isSelected() && !VariantHandler.allChromsfrom.isSelected() && Main.selectedChrom != 0) {
					varCalc();
				}
				else {
					
					varCalcBig();
				}
			}
			else {
				
				varCalc();
			}
			varcalc = false;
			}
			catch(Exception e) {
				e.printStackTrace();
				Main.drawCanvas.ready("Processing variants...");
			}
			Main.drawCanvas.ready("Processing variants...");
		}		
		else if(getreads) {
		
		
			Main.drawCanvas.loading("Loading reads");	
			
		
			readClass.loading = true;
			
			getReads(chrom, start, end, this.readClass);	
			readClass.loading = false;		
			splitIndex.updateReads = true;
	//		Main.drawCanvas.selectedSplit = readClass.sample.getReadClass().indexOf(readClass);
	
			getreads = false;			
			Main.drawCanvas.ready("Loading reads");
				
			Main.drawCanvas.repaint();
			
		}
		
		return null;
	}
	

	/*	File ref = new File("C:/HY-Data/RKATAINE/Rikurator/NewRator/genomes/hs37d5/hs37d5.fa");
		CRAMFileReader reader = new CRAMFileReader(new File("X:/cg8/Riku/HG00096.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram"), new File("X:/cg8/Riku/HG00096.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram.crai"),new ReferenceSource(ref));
		QueryInterval[] interval = { new QueryInterval(reader.getFileHeader().getSequence("chr2").getSequenceIndex(), 1000000, 2000000) };
		
		CloseableIterator<SAMRecord> iterator = reader.query(interval, true);
		
		while(iterator.hasNext()) {
			System.out.println(iterator.next().getSAMString());
		}*/
	//	fetchAnnotation();
//	fetchEnsembl();
	/*	try {
			
		      
			String tabixfile = "C:/HY-Data/RKATAINE/BasePlayer/BasePlayer/demo/Somatic-CRC/c32_5332_T_LP6005135-DNA_B01_somatic_GRCh37_20150513.vcf.gz";
		     
			
			
			TabixIndex index = new TabixIndex(new File(tabixfile+".tbi"));
			java.util.List<Block> blocks = index.getBlocks("1", 1000000, 4000000);
			System.out.println(blocks.get(0).getStartPosition());
			BlockCompressedInputStream input = new BlockCompressedInputStream(new File(tabixfile));
			input.seek(blocks.get(0).getStartPosition());
		//	gzip.skip(blocks.get(0).getStartPosition());
			 BufferedReader reader = new BufferedReader(new InputStreamReader(input));		
		//	 reader.skip(blocks.get(0).getStartPosition());
		      String line = reader.readLine();
		      System.out.println(line);
			
			reader.close();
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}*/
//	}
/*
TabixReader.Iterator getTabixIterator(String chrom, int start, int end, Sample sample) {
	try {
		
			tabixreader = new TabixReader(sample.getTabixFile());
			TabixReader.Iterator iterator = tabixreader.query(sample.vcfchr +chrom +":" +start+"-"+end);
				
			return iterator;
		
	}
	catch(Exception e) {
		
	//	JOptionPane.showMessageDialog(Main.chromDraw, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
		try {
			if(search) {
				if(chrom.equals("X")) {
					
					tabixreader = new TabixReader(sample.getTabixFile());
					iterator = tabixreader.query("23:" +searchStart+"-"+searchEnd);
				}
				else if(chrom.equals("Y")) {
					tabixreader = new TabixReader(sample.getTabixFile());
					iterator = tabixreader.query("24:" +searchStart+"-"+searchEnd);
				}
				else if(chrom.contains("M")) {
					tabixreader = new TabixReader(sample.getTabixFile());
					iterator = tabixreader.query("25:" +searchStart+"-"+searchEnd);
				}
			}
			else {
				if(chrom.equals("X")) {
					tabixreader = new TabixReader(sample.getTabixFile());
					iterator = tabixreader.query("23");
				}
				else if(chrom.equals("Y")) {
					tabixreader = new TabixReader(sample.getTabixFile());
					iterator = tabixreader.query("24");
				}
				else if(chrom.contains("M")) {
					tabixreader = new TabixReader(sample.getTabixFile());
					iterator = tabixreader.query("25");
				}
			}
		}
		catch(ArrayIndexOutOfBoundsException ex) {
			ErrorLog.addError(ex.getStackTrace());
			return null;
	//		System.out.println("Chromosome " +chrom +" not found in " +Main.drawCanvas.sampleList.get(i).getName());
		}
		catch(Exception ex) {
			FileRead.search = false;
			ErrorLog.addError(ex.getStackTrace());
			ex.printStackTrace();
			return null;
		}
		
	
		
	}
	return iterator;
}*/
static String setVCFFileStart(String chrom, int start, int end, Sample sample) {
	try {	Index index = null;	
			if(sample.getVCFReader() != null) {
				try {
					index = IndexFactory.loadIndex(sample.getTabixFile()+".idx");		
				}
				catch(Exception e) {
					index = IndexFactory.loadIndex(sample.getTabixFile()+".tbi");		
				}
			}
			else {
				index = new TabixIndex(new File(sample.getTabixFile()+".tbi"));
			}
			
			java.util.List<Block> blocks = null;
			
			if(index.containsChromosome(sample.vcfchr +chrom)) {
				chrom = sample.vcfchr +chrom;
			
				try {					
					blocks = index.getBlocks(chrom, start, end);	
					
				}
				catch(Exception e) {
				
					sample.vcfEndPos = 0;					
					return "";
				}
			
				if(blocks.size() > 0) {
					
					if(sample.getVCFReader() != null) {
						sample.setInputStream();
						sample.getVCFReader().skip(blocks.get(0).getStartPosition());
						
					}
					else {
						try {
							sample.getVCFInput().seek(0);
							sample.getVCFInput().seek(blocks.get(0).getStartPosition());
						}
						catch(Exception e) {
							ErrorLog.addError(e.getStackTrace());
							return "";
						}
					}
					
					
					sample.vcfEndPos = blocks.get(blocks.size()-1).getEndPosition();
				
					
				}
				else {
					
					if(sample.getVCFReader() != null) {
						sample.setInputStream();
					}
					else {
						sample.getVCFInput().seek(0);
					}
				}
				
			}
			else {
				
				if(index.containsChromosome(sample.vcfchr +(Main.chromosomeDropdown.getSelectedIndex()+1))) {
					try {
						blocks = index.getBlocks(sample.vcfchr+(Main.chromosomeDropdown.getSelectedIndex()+1), start, end);
					}
					catch(Exception e) {
						sample.vcfEndPos = 0;
						return "";
					}					
			
					if(blocks.size() > 0) {
						if(sample.getVCFReader() != null) {
							sample.setInputStream();
							sample.getVCFReader().skip(blocks.get(0).getStartPosition());
						}
						else {
							sample.getVCFInput().seek(0);
							sample.getVCFInput().seek(blocks.get(0).getStartPosition());
						}
					}
					else {
						if(sample.getVCFReader() != null) {
							sample.setInputStream();
						}
						else {
							sample.getVCFInput().seek(0);
						}
					}
					chrom = sample.vcfchr+(Main.chromosomeDropdown.getSelectedIndex()+1);
				}
				else {
					sample.vcfEndPos = 0;
				}
			}						
		
	}
	catch(Exception e) {
		e.printStackTrace();

		
	}
	return chrom;
}
void cancelFileRead() {
	changing = false;
	FileRead.search = false;
	bigcalc = false;
	Main.opensamples.setText("Add samples");
	head.putNext(null);
	current = null;
	Main.drawCanvas.current = FileRead.head;
	Main.drawCanvas.variantsStart = 0;
	Main.drawCanvas.variantsEnd = 1;
	Draw.updatevars = true;
	Main.drawCanvas.repaint();
	
	try {
		if(output != null) {
			output.close();
		}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	
}

void getVariantWindow(String chrom, int start, int end) {
	try {
	FileRead.lastpos = 0;
	removeNonListVariants();
	removeBedLinks();	
	
	for (int i = 0; i<Control.controlData.fileArray.size(); i++) {
	  Control.controlData.fileArray.get(i).controlled = false;
	}
	readFiles = true;	
	head.putNext(null);
	current = head;
	cancelvarcount = false;	
	cancelfileread = false;
	
	Main.drawCanvas.loadingtext = "Loading variants...";
	for(int i = 0; i<Main.samples; i++) {
		if(cancelvarcount || cancelfileread) {
			cancelFileRead();
			break;
		}
		if( Main.drawCanvas.sampleList.get(i).getTabixFile() == null ||  Main.drawCanvas.sampleList.get(i).multipart) {
			continue;
		}
		current = head;
		
		getVariants(chrom,start,end, Main.drawCanvas.sampleList.get(i));
		
	}
	readFiles =false;
	annotate();
	
	if(Control.controlData.controlsOn) {		    	
		Control.applyControl();
	}
	Main.bedCanvas.intersected = false;
	if(bigcalc) {
		Main.drawCanvas.calcClusters(FileRead.head);
	}
	else {
		Main.drawCanvas.calcClusters(FileRead.head,1);
	}	
	
	if(Main.bedCanvas.bedOn) {
		for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
			Main.bedCanvas.bedTrack.get(i).used = false;
			if(Main.bedCanvas.bedTrack.get(i).small && Main.bedCanvas.bedTrack.get(i).getZoomlevel() == null) {						
				if( Main.bedCanvas.bedTrack.get(i).intersect && !Main.bedCanvas.bedTrack.get(i).loading) {				
					Main.bedCanvas.annotate(Main.bedCanvas.bedTrack.get(i).getHead());					
				}							
				else if(Main.bedCanvas.bedTrack.get(i).intersect && Main.bedCanvas.bedTrack.get(i).loading) {
					Main.bedCanvas.bedTrack.get(i).waiting = true;							
				}
			}
			else if(Main.bedCanvas.bedTrack.get(i).intersect) {	
			
				BedCanvas.Annotator annotator = Main.bedCanvas.new Annotator(Main.bedCanvas.bedTrack.get(i));
				annotator.annotateVars();				
			}
		}
		Main.bedCanvas.intersected = true;
		if(FileRead.bigcalc) {
			Main.drawCanvas.calcClusters(FileRead.head);
		}
		else {
			Main.drawCanvas.calcClusters(FileRead.head,1);
			
		}
	}			
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

static String getVCFLine(String chrom, int start, int end, Sample sample) {
	if(sample.calledvariants) {
		StringBuffer altbases = new StringBuffer("");
		int alts=0, refs=0;
		
		for(int i = 0 ; i<Main.drawCanvas.varOverLap.vars.size()-1; i++) {
			altbases.append(Main.drawCanvas.varOverLap.vars.get(i).getKey());
			
			
		}
		alts = Main.drawCanvas.varOverLap.vars.get(0).getValue().get(0).getCalls();
		refs = (Main.drawCanvas.varOverLap.vars.get(0).getValue().get(0).getCoverage()-Main.drawCanvas.varOverLap.vars.get(0).getValue().get(0).getCalls());
		String genotype = "";
		if(alts/(double)refs > 0.95) {
			genotype = "1/1";
		}
		else {
			genotype = "0/1";
		}
		altbases.append(Main.drawCanvas.varOverLap.vars.get(Main.drawCanvas.varOverLap.vars.size()-1).getKey());
		return Main.drawCanvas.splits.get(0).chrom +"\t" +(Main.drawCanvas.varOverLap.getPosition()+1) +"\t.\t" +Main.getBase.get(Main.drawCanvas.varOverLap.getRefBase()) +"\t" +altbases +"\t99\tPASS\tINFO\tGT:AD:DP\t" +genotype +":" +refs +"," +alts +":" +(refs+alts);
	}
	String line = "";
	
	cancelfileread = false;
	cancelvarcount = false;
	
	try {	
		if(!(VariantHandler.hideIndels.isSelected() && VariantHandler.hideSNVs.isSelected())) { 		
			if(sample.multipart) {
				for(int i = sample.getIndex(); i>= 0; i--) {
					if(!Main.drawCanvas.sampleList.get(i).multipart) {
						sample = Main.drawCanvas.sampleList.get(i);
						break;
					}
				}
			}
			if(sample.getTabixFile() != null) {
			
				setVCFFileStart(chrom, start, end+3, sample);
				boolean vcf = sample.getVCFReader() != null;
				while(line != null) {
					if(vcf) {
						try {
							sample.getVCFReader().ready();
						}
						catch(IOException ex) {
							
						}
					}					
					try {
						if(vcf) {
							line = sample.getVCFReader().readLine();							
						}
						else {							
							line = sample.getVCFInput().readLine();
						}
						if(Integer.parseInt(line.split("\t")[1]) == start+1 ) {
							if(sample.getVCFReader() != null) {
								sample.getVCFReader().close();
							}
							return line;
						}
 						
					}
					catch(Exception ex) {							
						Main.showError(ex.getMessage(), "Error");
						
						ErrorLog.addError(ex.getStackTrace());
						ex.printStackTrace();
						Main.cancel();
						changing = false;
						
					}											
									
				}
				
				if(sample.getVCFReader() != null) {
					sample.getVCFReader().close();
				}
				
				return line;	
			}		
		}		 
	}
	catch(Exception exc) {
		Main.showError(exc.getMessage(), "Error");
		
		System.out.println(sample.getName());
		exc.printStackTrace();
		ErrorLog.addError(exc.getStackTrace());
		changing = false;		
	}
	return "";
}

void getVariants(String chrom, int start, int end, Sample sample) {
	if(sample.calledvariants) {
		Main.drawCanvas.variantsStart = start;
		Main.drawCanvas.variantsEnd = end;
		return;
	}
	String line;
	String[] split;
	cancelfileread = false;
	cancelvarcount = false;
	
	
	try {	
		if(!(VariantHandler.hideIndels.isSelected() && VariantHandler.hideSNVs.isSelected())) { 		
			readFiles = true;	
			Main.drawCanvas.splits.get(0).transStart = 0;
			Main.drawCanvas.variantsStart = start;
			Main.drawCanvas.variantsEnd = end;
			Main.drawCanvas.loadbarAll = (int)((sample.getIndex()/(double)(Main.samples))*100);
			
			linecounter = 0;
			if(cancelfileread) {		
				cancelFileRead();
				return;
			}
			if(sample.multipart){
				return;
			}
			
			if(sample.getTabixFile() != null) {
				
				String searchChrom = setVCFFileStart(chrom, start, end, sample);
				boolean vcf = sample.getVCFReader() != null;
				if(vcf) {
					try {
						sample.getVCFReader().ready();
					}
					catch(IOException ex) {
						return;
					}
				}
				
				sample.setMaxCoverage(0F);		
				current = head;
				line = "";
				first = true;
				stop = false;
			
				while(true) {	    
					if(stop) {
						
						break;
					}
					if(cancelfileread || cancelvarcount || !Main.drawCanvas.loading) {
						cancelFileRead();
						
						break;
			  		}
					
					try {
						
						if(vcf) {
							line = sample.getVCFReader().readLine();
						
							if(line != null && line.startsWith("#")) {
								continue;
							}
						}
						else {
							try {
								line = sample.getVCFInput().readLine();
								
								
							}
							catch(htsjdk.samtools.FileTruncatedException e) {
								e.printStackTrace();
								
							}
							
							
						}
						if(line == null) {
							break;
						}
					}
					catch(Exception ex) {							
						Main.showError(ex.getMessage(), "Error");
						
						ErrorLog.addError(ex.getStackTrace());
						ex.printStackTrace();
						Main.cancel();
						changing = false;
						break;
					}												
					
					split = line.split("\\t+");
					try {
					if(split[0].startsWith("#")) {
						continue;
					}
					if(split.length >2 && (Integer.parseInt(split[1]) > end || !split[0].equals(searchChrom))) {
						
						break;
					}
					}
					catch(Exception e) {
						String string = "";
						for(int i = 0 ; i< split.length; i++) {
							string += split[i]+"\t";
						}
						
						ErrorLog.addError(string);
					}
					
					if(sample.multiVCF) {
						
						readLineMulti(split, sample);							
					}
					else {
					
						readLine(split, sample);
						
					}
					if(first) {
						first = false;
					}				

					if(sample.getVCFInput() != null) {
						
							if(sample.getVCFInput().getFilePointer() > sample.vcfEndPos ) {		
							
								break;
							}
					}
					
				}
				if(sample.getVCFReader() != null) {
					sample.getVCFReader().close();
				}
				Draw.updatevars = true;
				if(sample.getIndex()*Main.drawCanvas.drawVariables.sampleHeight < Main.drawScroll.getViewport().getHeight()+Main.drawCanvas.drawVariables.sampleHeight) {
					if(Main.drawCanvas.splits.get(0).viewLength < 2000000) {						
						Main.chromDraw.updateExons = true;
						Main.chromDraw.repaint();		
					}
				}			
			}		
		}		 
	}
	catch(Exception exc) {
		Main.showError(exc.getMessage(), "Error");
		
		System.out.println(sample.getName());
		exc.printStackTrace();
		ErrorLog.addError(exc.getStackTrace());
		changing = false;		
	}
	
}
	
	public void changeChrom(String chrom) {
		
		try {
			
			nobeds = false;			
			cancelfileread = false;
			
		try {	
			Main.drawCanvas.loading("Loading annotation...");		
			
			Main.drawCanvas.splits.get(0).setGenes(getExons(chrom));	
			
			Main.chromDraw.updateExons = true;
			Main.chromDraw.repaint();
			
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());			
			e.printStackTrace();
		}
		
		Main.drawCanvas.ready("Loading annotation...");
		if(Main.bedCanvas.bedTrack.size() > 0) {
			
			Main.drawCanvas.loading("Loading BED-files...");
			if(search) {				
				
				for(int i = 0; i< Main.bedCanvas.bedTrack.size(); i++) {
					Main.bedCanvas.bedTrack.get(i).used = false;
					
					if(Main.bedCanvas.bedTrack.get(i).small && Main.bedCanvas.bedTrack.get(i).getBBfileReader() == null) {
						
							Main.bedCanvas.getBEDfeatures(Main.bedCanvas.bedTrack.get(i), 1, Main.drawCanvas.splits.get(0).chromEnd);
					}
					else {
						if(searchEnd- searchStart < Settings.windowSize) {
							Main.bedCanvas.getBEDfeatures(Main.bedCanvas.bedTrack.get(i), searchStart, searchEnd);
						}
					}
				}
				
			}
			else {
				
				for(int i = 0; i< Main.bedCanvas.bedTrack.size(); i++) {
					
					Main.bedCanvas.bedTrack.get(i).used = false;
					if(Main.bedCanvas.bedTrack.get(i).small  && Main.bedCanvas.bedTrack.get(i).getBBfileReader() == null) {			
						
						Main.bedCanvas.getBEDfeatures(Main.bedCanvas.bedTrack.get(i), 1, Main.drawCanvas.splits.get(0).chromEnd);
										
					}	
					if(nobeds) {
						Main.drawCanvas.ready("Loading BED-files...");
						
						return;
					}
				}
			}
			
		}
		
		Main.drawCanvas.ready("Loading BED-files...");
		
		if(novars) {
			Main.drawCanvas.variantsStart= 0;
			Main.drawCanvas.variantsEnd = 0;
			
		}
		else {
			changing = true;
		}
		 
		if(Main.varsamples > 0 && !novars && !bigcalc) {			
			removeNonListVariants();
			removeBedLinks();
			
			Main.drawCanvas.loading("Loading variants...");
			head.putNext(null);
			current = FileRead.head;		
			if(FileRead.head.getPosition() > 0) {
				FileRead.head =  new VarNode(0,  (byte)0,"N", (short)0, (short)0, false,(float)0,(float)0,null,null, null, null, null);     
			}
			Main.drawCanvas.current = head;		
			Main.chromDraw.varnode = null;
			Main.chromDraw.vardraw = null;
			
			for(int i = 0; i<Main.samples; i++) {
				if(nobeds) {
					return;
				}
				if(cancelfileread || !Main.drawCanvas.loading) {
					
					cancelFileRead();
					break;
				}
				if( Main.drawCanvas.sampleList.get(i).getTabixFile() == null ||  Main.drawCanvas.sampleList.get(i).multipart) {
					continue;
				}
				
				if(search) {		
					
					getVariants(chrom, FileRead.searchStart, FileRead.searchEnd, Main.drawCanvas.sampleList.get(i));
				}
				else {
					
					getVariants(chrom, 0, Main.drawCanvas.splits.get(0).chromEnd, Main.drawCanvas.sampleList.get(i));
				}	
				
			}
			
			annotate();
			

		  	Main.drawCanvas.loading("Applying controls...");
			if(Control.controlData.controlsOn) {		    	
			    Control.applyControl();
			}
			Main.drawCanvas.ready("Applying controls...");
			readFiles = false;
			Main.bedCanvas.intersected = false;
			if(bigcalc) {
				Main.drawCanvas.calcClusters(FileRead.head);
			}
			else {
				Main.drawCanvas.calcClusters(FileRead.head,1);
				
			}
			Main.bedCanvas.bedOn = true;
			boolean ison = false;
			for(int i = 0 ; i<Main.bedCanvas.bedTrack.size();i++) {
				if(Main.bedCanvas.bedTrack.get(i).intersect) {
				
					ison = true;
					break;
				}
			}
			if(!ison) {
				Main.bedCanvas.bedOn = false;
			}
		
			if(Main.bedCanvas.bedOn) {
			
				Main.drawCanvas.loadingtext = "Annotating variants";
				for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
					if(Main.bedCanvas.bedTrack.get(i).small && Main.bedCanvas.bedTrack.get(i).getBBfileReader() == null) {						
						if( Main.bedCanvas.bedTrack.get(i).intersect && !Main.bedCanvas.bedTrack.get(i).loading) {
							Main.bedCanvas.annotate(Main.bedCanvas.bedTrack.get(i).getHead());	
							Main.bedCanvas.intersected = true;
							
						}							
						else if(Main.bedCanvas.bedTrack.get(i).intersect && Main.bedCanvas.bedTrack.get(i).loading) {
							Main.bedCanvas.bedTrack.get(i).waiting = true;							
						}
					}
					else if(Main.bedCanvas.bedTrack.get(i).intersect) {		
						
						BedCanvas.Annotator annotator = Main.bedCanvas.new Annotator(Main.bedCanvas.bedTrack.get(i));
						annotator.annotateVars();
						
						Main.bedCanvas.intersected = true;
						
					}
				}	
				
			}
			
			if(FileRead.bigcalc) {
				Main.drawCanvas.calcClusters(FileRead.head);
			}
			else {
				Main.drawCanvas.calcClusters(FileRead.head,1);
				
			}
		}
		
		
		Draw.updatevars = true;		  	
	  	Main.drawCanvas.ready("Loading variants...");
	  	if(!Main.drawCanvas.loading) {
			Draw.calculateVars = true;
		}
		if(novars) {
			Main.drawCanvas.variantsStart= 0;
			Main.drawCanvas.variantsEnd = Main.drawCanvas.splits.get(0).chromEnd;
		}
		FileRead.novars = false;
		search = false;
		
		
		changing = false;
	    current = null;
		Main.chromDraw.updateExons = true;
		Main.chromDraw.repaint();
		
		}
		catch(Exception e) {
			e.printStackTrace();
			ErrorLog.addError(e.getStackTrace());
			changing = false;
		}
		
		Main.drawCanvas.loadbarAll = 0;
		Main.drawCanvas.loadBarSample= 0;
		Main.drawCanvas.repaint();
	}
	
	public ArrayList<Gene> getExons(String chrom) {	
		
		ArrayList<Gene> transcriptsTemp = new ArrayList<Gene>();
		
		try {			
			if(Main.genomehash.size() == 0 || Main.genomehash.get(Main.defaultGenome).size() == 0) {
				return new ArrayList<Gene>();
			}
			if(ChromDraw.exonReader != null) {
				ChromDraw.exonReader.close();
				
			}
			ChromDraw.exonReader = new TabixReader(Main.genomehash.get(Main.defaultGenome).get(Main.annotation).getCanonicalPath());
			
			if(chrom == null) {
				
				return null;
			}
			
			TabixReader.Iterator exonIterator = null;
				try {	
					
					if(!ChromDraw.exonReader.getChromosomes().contains(chrom)) {
						
						String[] gene = { Main.chromosomeDropdown.getSelectedItem().toString(), "1", ""+Main.drawCanvas.splits.get(0).chromEnd, Main.chromosomeDropdown.getSelectedItem().toString(), "1", "+", "-", "-", "-", "-","-","1","1","1",""+Main.drawCanvas.splits.get(0).chromEnd,"-1,","-"};
						Gene addGene = new Gene(gene);
						Transcript addtrans = null;
						try {
							addtrans = new Transcript(gene);
						}
						catch(Exception e) {
							e.printStackTrace();
						}
						addGene.addTranscript(addtrans);
						addGene.setLongest(addtrans);
						
						addtrans.setGene(addGene);
						transcriptsTemp.add(addGene);	
						
						return transcriptsTemp;
						//return new ArrayList<Gene>();
					}
					else {
						exonIterator = ChromDraw.exonReader.query(chrom);
					}				
				}
				catch(Exception e) {
					try {
						
						if(chrom.matches("\\w+")) {
							exonIterator = ChromDraw.exonReader.query("M");
							
						}
					}
					catch(Exception ex) {						
						System.out.println(chrom);		
						
						e.printStackTrace();
					}					
				}
				
			String s;
			String[] exonSplit;
			Transcript addtrans = null;
			Hashtable<String, Gene> genes = new Hashtable<String, Gene>();
			Gene setGene;
			while(exonIterator != null && (s = exonIterator.next()) != null) {	
				
				exonSplit = s.split("\t");
				if(exonSplit[0].equals("23")) {
					exonSplit[0] = "X";
				}
				else if(exonSplit[0].equals("24")) {
					exonSplit[0] = "Y";
				}
				else if (exonSplit[0].equals("25")) {					
					exonSplit[0] = "MT";
				}			
				addtrans = new Transcript(exonSplit);
				
				if(!genes.containsKey(exonSplit[6])) {	
					Gene addgene = new Gene(exonSplit);
					genes.put(exonSplit[6],addgene);				
					addgene.addTranscript(addtrans);
					addgene.setLongest(addtrans);
					addtrans.setGene(addgene);
					transcriptsTemp.add(addgene);														
				}
				else {
					setGene = genes.get(exonSplit[6]);
					setGene.addTranscript(addtrans);
					addtrans.setGene(setGene);		
					if(addtrans.getLength() > setGene.getLongest().getLength()) {
						setGene.setLongest(addtrans);
					}
				}					
			}
			genes.clear();
			ChromDraw.exonReader.close();
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
		}
		
			
		
		return transcriptsTemp;
	}
	
	
	
	private void readBAM(File[] files) {
		try {
		 
		 File addFile =null;
	  	 File[] addDir;
	  	 Sample currentSample = null;
	  	 Boolean added = false;
	
	  	 for(int i = 0; i<files.length; i++) {
	  		 if(!files[i].exists()) {
	  			 continue;
	  		 }
		  	if(files[i].isDirectory()) {
		  		addDir = files[i].listFiles();
		  		for(int f = 0; f<addDir.length; f++) { 	
		  			if(addDir[f].getName().endsWith(".bam") ||addDir[f].getName().endsWith(".cram")  ) {
		  				addFile = addDir[f];
		  				
		  				break;
		  			}
		  		}
		  		
		  	}
		  	else {
		  		addFile = files[i];
		  		
		  	}
		  
	  		 if(addFile != null) {
	  			  Main.drawCanvas.bam = true;
		      	  currentSample = new Sample(addFile.getName(), (Short)Main.samples, null);	      	 
		      	  Main.drawCanvas.sampleList.add(currentSample);	      	 
		      	  currentSample.samFile = addFile;
		      	  currentSample.resetreadHash();
		      	  if(addFile.getName().endsWith(".cram")) {
		      		  currentSample.CRAM = true;
		      	  }
		      	 
			      	if(currentSample.samFile.getName().endsWith(".cram")) {
			      		currentSample.readString = "CRAM";
				  	}
				  	else {
				  		currentSample.readString = "BAM";
				  	}
		      try {
		    	  if(currentSample.samFile.getName().endsWith(".cram")) {
		    		/*
		    		 CRAMReader = new CRAMFileReader(currentSample.samFile, new File(currentSample.samFile.getCanonicalFile() +".crai"),new ReferenceSource(Main.ref), ValidationStringency.LENIENT);
		    		 
		    		 if(CRAMReader.getFileHeader().getTextHeader().contains("SN:chr")) {
		    			 currentSample.chr = "chr";
		    		 }
		    		 if(CRAMReader.getFileHeader().getTextHeader().contains("cgatools") || CRAMReader.getFileHeader().getTextHeader().contains("CompleteGenomics")) {
		    			 currentSample.complete = true;
		    		 }
		    		
		    			CRAMReader.close();
		    */
		    	  }
		    	  else {
		    		 
		    	//	  samFileReader = SamReaderFactory.makeDefault().open(currentSample.samFile);
		    	/*	  samFileReader = new SAMFileReader(currentSample.samFile);
			      	  if(samFileReader.getFileHeader().getTextHeader().contains("SN:chr")) {
			      		  currentSample.chr = "chr";
			      	  }
			      	  if(samFileReader.getFileHeader().getTextHeader().contains("cgatools") || samFileReader.getFileHeader().getTextHeader().contains("CompleteGenomics")) {
			      		currentSample.complete = true;
		    		 }
			      	 samFileReader.close();
		    	  */
		    	  }
		      	 }
		      	 catch(Exception e) {
		      		ErrorLog.addError(e.getStackTrace());
		      		 e.printStackTrace();
		      	 }
		      		added = true;
		      		Main.readsamples++;
				  Main.samples++;
		         }   
	  		 
	  	 }
	  	if(!added) {
	  		return;
	  	}
	  	checkSamples();
	   	 Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.samples);
	  //	 Main.drawCanvas.checkSampleZoom();
	  	 
	  	 if(Main.drawScroll.getViewport().getHeight()/(Main.drawCanvas.sampleList.size()) > Draw.defaultSampleHeight) {
	 		  Main.drawCanvas.drawVariables.sampleHeight = Main.drawScroll.getViewport().getHeight()/Main.drawCanvas.sampleList.size();
	 	  }
	 	  else {
	 		  Main.drawCanvas.drawVariables.sampleHeight = Draw.defaultSampleHeight;  		  
	 	  }
	 	
	 	 if(Main.drawCanvas.getHeight() < (Main.drawCanvas.sampleList.size())*Main.drawCanvas.drawVariables.sampleHeight) {  		   

	 		 Main.drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), (int)((Main.drawCanvas.sampleList.size())*Main.drawCanvas.drawVariables.sampleHeight));
	 		 Main.drawCanvas.revalidate();
	 	 }
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
		}
		Draw.updateReads = true;
		 Draw.updatevars = true;
		 Main.drawCanvas.repaint();
}	  	

static boolean checkIndex(File file) {
	try {
		if(file.getName().endsWith(".vcf")) {
			return new File(file.getCanonicalPath() +".idx").exists() || new File(file.getCanonicalPath() +".tbi").exists();			
		}
		else if(file.getName().endsWith(".vcf.gz")) {
			return new File(file.getCanonicalPath() +".tbi").exists();
		}
		else if(file.getName().endsWith(".bam")) {
			if(!new File(file.getCanonicalPath() +".bai").exists()) {
				return new File(file.getCanonicalPath().replace(".bam", "") +".bai").exists();
			}
			else {
				return true;
			}			
		}
		else if(file.getName().endsWith(".bed.gz") || file.getName().endsWith(".gff.gz") || file.getName().endsWith(".bed") || file.getName().endsWith(".gff")) {
			
			return new File(file.getCanonicalPath() +".tbi").exists();
		}
		else {
			return true;
		}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	
	return false;
}

static void checkMulti(Sample sample) {
	try {
	  Sample addSample;
	  BufferedReader reader = null;
	  GZIPInputStream gzip = null;
	  String line;
	  
	  if(sample.getTabixFile().endsWith(".gz")) {
		  try {
			  gzip = new GZIPInputStream(new FileInputStream(sample.getTabixFile()));
			  reader = new BufferedReader(new InputStreamReader(gzip));		
		  }
		  catch(Exception e) {
			  Main.showError("Corrupted file: " +sample.getTabixFile() +"\nTry to bgzip and recreate index file.", "Error");
			  Main.drawCanvas.sampleList.remove(sample);
			  Main.varsamples--;
			  Main.samples--;
			
		  }
	  }
	  else {
		  reader = new BufferedReader(new FileReader(sample.getTabixFile()));		
	  }
	
      line = reader.readLine();
     
	if(!sample.multipart && line != null) {
		
		while(line != null && line.startsWith("#")) {
			try {
			if(line.startsWith("##INFO")) {
				if(line.contains("Type=Float") || line.contains("Type=Integer")) {					
					VariantHandler.addMenuComponents(line);					
				}
			}
			if(line.startsWith("##FILTER")) {
				if(line.contains("ID=") || line.contains("Description=")) {					
					VariantHandler.addMenuComponents(line);					
				}	
			}
			if(line.startsWith("#CHROM")) {
				headersplit = line.split("\t+");
				
				if(headersplit.length  > 11) {
					sample.multiVCF = true;
					Main.varsamples--;
					for(int h = 9; h<headersplit.length; h++) {
						 addSample = new Sample(headersplit[h], (short)(Main.samples), null);  
			        	 addSample.multipart = true;			        	 
			        	 Main.drawCanvas.sampleList.add(addSample);				          	
			          	 Main.samples++;
			          	 Main.varsamples++;
			          	 if(sampleString == null) {
			          		 sampleString = new StringBuffer("");
			          	 }
			          	 sampleString.append(addSample.getName() +";");		
			          	
					}
					VariantHandler.commonSlider.setMaximum(Main.varsamples);
		          	VariantHandler.commonSlider.setUpperValue(Main.varsamples);
		          	VariantHandler.geneSlider.setMaximum(Main.varsamples);
					Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.drawCanvas.sampleList.size());	
					Main.drawCanvas.checkSampleZoom();
				  	Main.drawCanvas.resizeCanvas(Main.drawScroll.getViewport().getWidth(),Main.drawScroll.getViewport().getHeight());
				}
				
			}					
			}
			catch(Exception ex) {
				ex.printStackTrace();
			}
			line = reader.readLine();
		}
	line = reader.readLine();
	if(line != null) {
		if(line.startsWith("chr")) {
			sample.vcfchr = "chr";
		}
		
		headersplit = line.split("\t");
		
		if(headersplit.length == 11) {
			
			if(headersplit[10].startsWith("0:") || (headersplit[10].charAt(0) == '0' && headersplit[10].charAt(2) == '0')) {
				sample.somaticColumn = 9;
			}
			else {
				sample.somaticColumn = 10;
			}
		}
	/*	if(headersplit.length > 7 && headersplit[7].contains("=")) {
			String[] infosplit = headersplit[7].split(";");
			for(int i = 0; i<infosplit.length; i++) {
				if(infosplit[i].contains("=")) {
					if(Main.drawCanvas.drawVariables.advQualities == null) {
						Main.drawCanvas.drawVariables.advQualities = new HashMap<String, Float>();
					}
					if(!Main.drawCanvas.drawVariables.advQualities.containsKey(infosplit[i].substring(0,infosplit[i].indexOf("=")))) {
						try {
							Float.parseFloat(infosplit[i].substring(infosplit[i].indexOf("=")+1));
							
						}
						catch(Exception e) {
							continue;
						}	
						VariantHandler.addMenuComponents(infosplit[i].substring(0,infosplit[i].indexOf("=")));						
					}
				}
			}			
		}*/
	}
	
	checkSamples();
	line = null;
	reader.close();
	if(gzip != null) {
		gzip.close();
	}
	}
	else {
		reader.close();
		if(gzip != null) {
			gzip.close();
		}
	}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	 
}

public static class SearchBamFiles extends SwingWorker<String, Object> {		
	boolean cram = false;
	SamReader samFileReader;	
	int sampletemp, fileindex;
	
	CRAMFileReader CRAMReader = null;
	File[] files;
	public SearchBamFiles(File[] files, int fileindex, int sampletemp) {
		this.files = files;
		this.fileindex = fileindex;
		this.sampletemp = sampletemp;
	}
	protected String doInBackground() {		
		try {
			File[] bamfiles = null;
			searchingBams = true;
			if(fileindex > files.length) {
				bamfiles = files[0].listFiles(new FilenameFilter() {
				 public boolean accept(File dir, String name) {
				        return name.toLowerCase().endsWith(".bam") || name.toLowerCase().endsWith(".cram") || name.toLowerCase().endsWith(".list");
				     }
		  		 });
				
			}
			else {
				 bamfiles = files[fileindex].getParentFile().listFiles(new FilenameFilter() {
			     public boolean accept(File dir, String name) {
			        return name.toLowerCase().endsWith(".bam") || name.toLowerCase().endsWith(".cram") || name.toLowerCase().endsWith(".list");
			     }
	  		 });
			}
	    /* File[] cramfiles = files[fileindex].getParentFile().listFiles(new FilenameFilter() {
	         public boolean accept(File dir, String name) {
	            return name.toLowerCase().endsWith(".cram");
	         }
	       	 });
	  		*/
	  	 if(bamfiles.length > 0) {
	  		
	  		int index = -1, sampleindex;
	  		 for(int i = 0; i<bamfiles.length; i++) {
	  			 
	  			 cram = false;
		  		 sampleindex = 0;
		  		 index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".")));
		  		
		  		
		  		
		  	 	 if (index < 0) continue; 
		  	 	if(!checkIndex(bamfiles[i])) {
		   			 Main.putMessage("Check Tools->View log");
		   			 ErrorLog.addError("No index file found for " +bamfiles[i].getName());
		   			 continue;
		  	    	}
		  	 	 for(char letter : sampleString.substring(0, index).toString().toCharArray()) {
		  	 		 if(letter == ';') sampleindex++;
		  	 	 }
		  	 	Main.drawCanvas.bam = true;
		  	 	Main.readsamples++;
		  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile = new File(bamfiles[i].getCanonicalPath());	 
		  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).resetreadHash();
		  	 	
		  	 	
		  	 	if(Main.readsamples==1) {
		  	 		checkSamples();
		  	 	}
		  	 	//samFileReader = SamReaderFactory.makeDefault().open(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
			  	  if(bamfiles[i].getName().endsWith(".cram")){
			  		  try {
			  			
			  		//	CRAMReader = new CRAMFileReader(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile, new File(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile.getCanonicalFile() +".crai"),new ReferenceSource(Main.ref), ValidationStringency.LENIENT);
			  			cram = true;
			    	/*	 if(CRAMReader.getFileHeader().getTextHeader().contains("SN:chr")) {
			    			 Main.drawCanvas.sampleList.get(sampleindex+sampletemp).chr = "chr";
			    		 }
			    		 if(CRAMReader.getFileHeader().getTextHeader().contains("cgatools") || CRAMReader.getFileHeader().getTextHeader().contains("CompleteGenomics")) {
			    			 Main.drawCanvas.sampleList.get(sampleindex+sampletemp).complete = true;
			    		 }
			    		
			    			CRAMReader.close(); */
			  	  
			  		  }
			  		  catch(Exception e) {
			  			  ErrorLog.addError(e.getStackTrace());
			  		  
			  		  }
			  	  }
			  	  else {
			  		cram = false;
			  		/*samFileReader = new SAMFileReader(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
			  	    if(samFileReader.getFileHeader().getTextHeader().contains("SN:chr")) {
		      	    	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).chr = "chr";
		      	    }
				  	  if(samFileReader.getFileHeader().getTextHeader().contains("cgatools") || samFileReader.getFileHeader().getTextHeader().contains("CompleteGenomics")) {
				  		Main.drawCanvas.sampleList.get(sampleindex+sampletemp).complete = true;
			  		 }
			  		 */
			  	  }
			  	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).CRAM = cram;
			  	if(cram) {
			  		Main.drawCanvas.sampleList.get(sampleindex+sampletemp).readString = "CRAM";
			  	}
			  	else {
			  		Main.drawCanvas.sampleList.get(sampleindex+sampletemp).readString = "BAM";
			  	}
			  	  if(samFileReader != null) {
			      	samFileReader.close();
			  	  }
			  	 }  	
	  	 	}
	  	 
		  	sampleString = null;
		 	 files = null;
		}
		catch(Exception e) {
			searchingBams = false;
			e.printStackTrace();
		}
		searchingBams = false;
		Main.drawCanvas.repaint();
		return "";
	}
	
}

private void readVCF(File[] files) {
		try {	 
		//  boolean cram = false, diffPaths = false;
	  	
	  	  Main.drawCanvas.loading("Loading samples...");
	  	  File[] addDir;
	  	  int sampletemp = Main.samples;
	  	  Boolean added = false;
		   	Sample addSample = null;
		  	sampleString = new StringBuffer(""); 
		  	
		  	int fileindex = -1;
		  	readFiles = true;
		  	cancelfileread = false;
	  	
	 if(Control.controlData.controlsOn) {	    	
	    Control.dismissControls(head);
	 }
	
  	 for(int fi = 0; fi<files.length; fi++) { 
  		
  		if(Main.cancel || !Main.drawCanvas.loading) {
  			current = null;
  			FileRead.head.putNext(null);
  			return;
  		}
  		if(!files[fi].exists()) {
  			continue;
  		}
  		
  	    addDir = null;  	   
  	   
  	    if(files[fi].isDirectory()) {  	    	
  	    	
  	    	addDir = files[fi].listFiles(new FilenameFilter() {
  	    	     public boolean accept(File dir, String name) {
  	    	        return name.toLowerCase().endsWith(".vcf.gz") || name.toLowerCase().endsWith(".vcf");
  	    	     }
  	    	});
  	    	
  	    	for(int f= 0; f<addDir.length; f++) {
  	    		if(cancelfileread || !Main.drawCanvas.loading) {
  	    			current = null;
  	    			FileRead.head.putNext(null);
  	    			return;
  	    		}
  	    		 if(!checkIndex(addDir[f])) {
  	    			 Main.putMessage("Check Tools->View log");
  	    			 ErrorLog.addError("No index file found for " +addDir[f].getName());	    			 
  	    		 }
  	    		 addSample = new Sample(addDir[f].getName(), (Short)Main.samples, addDir[f].getCanonicalPath());      
  	    		 Main.drawCanvas.sampleList.add(addSample);
  	        	 Main.varsamples++;
  	        	 Main.samples++;
  	        	 Main.drawCanvas.drawVariables.visiblesamples++;
  	        	 checkMulti(addSample);
  	        	 Main.drawCanvas.loadingtext = "Loading samples... " +f +"/" +addDir.length;
  	        	 
  	        	 added = true;
  	        	 VariantHandler.commonSlider.setMaximum(Main.varsamples);
  	        	 VariantHandler.commonSlider.setUpperValue(Main.varsamples);
  	        	 VariantHandler.geneSlider.setMaximum(Main.varsamples);
  	        	 sampleString.append(addSample.getName() +";");
  	        	 fileindex = f;
  	    		 }
  	    	/*	 File[] bamfiles = files[fi].listFiles(new FilenameFilter() {
  	    		 
  	    	     public boolean accept(File dir, String name) {
  	    	        return name.toLowerCase().endsWith(".bam") || name.toLowerCase().endsWith(".cram") || name.toLowerCase().endsWith(".list");
  	    	     }
  	    	   	 });
  	    		
  	    	  	 if(bamfiles.length > 0) {
  	    	  		
  	    	  		int index = -1, sampleindex;
  	    	  		 for(int i = 0; i<bamfiles.length; i++) {
  	    	  			cram = false;
  	    		  		 sampleindex = 0;
  	    		  		// if(bamfiles[i].getName().indexOf(".bam") > 0) {
  	    		  		 
  	    		  			index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".")));
  	    		  		if(bamfiles[i].getName().endsWith(".cram") ) {
  	    		  			cram = true;
  	    		  		}
  	    		  		
  	    		  	 	 if (index < 0) continue; 
  	    		  	 	 
  	    		  	 	 for(char letter : sampleString.substring(0, index).toString().toCharArray()) {
  	    		  	 		 if(letter == ';') sampleindex++;
  	    		  	 	 }
  	    		  	 	 Main.drawCanvas.bam = true;
  	    		  	 	 Main.drawCanvas.sampleList.get(sampleindex+sampletemp).CRAM = cram;
  	    		  	 	 Main.drawCanvas.sampleList.get(sampleindex+sampletemp).resetreadHash();
  	    		  	 	 Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile = new File(bamfiles[i].getAbsolutePath());	 
 	    			  	  //samFileReader = SamReaderFactory.makeDefault().open(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
 	    			  	  samFileReader = new SAMFileReader(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
  	    		      	  if(samFileReader.getFileHeader().getTextHeader().contains("SN:chr")) {
  	    		      		Main.drawCanvas.sampleList.get(sampleindex+sampletemp).chr = "chr";
  	    		      	  }  	    		      	 
  	    		  	 }  	
  	    	  	 }*/
  	    	  	/* else if(cramfiles.length > 0) {
  	    	  		int index = -1, sampleindex;
 	    	  		 for(int i = 0; i<bamfiles.length; i++) {
 	    		  		 sampleindex = 0;
 	    		  		 index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".cram")));
 	    		  	 	 if (index < 0) continue; 
 	    		  	 	 
 	    		  	 	 for(char letter : sampleString.substring(0, index).toString().toCharArray()) {
 	    		  	 		 if(letter == ';') sampleindex++;
 	    		  	 	 }
 	    		  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).CRAM = true;
 	    		  	 	 Main.drawCanvas.sampleList.get(sampleindex+sampletemp).resetreadHash();
 	    		  	 	 Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile = new File(bamfiles[i].getAbsolutePath());	 
	    			  	  //samFileReader = SamReaderFactory.makeDefault().open(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
	    			  	  samFileReader = new SAMFileReader(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
 	    		      	  if(samFileReader.getFileHeader().getTextHeader().contains("SN:chr")) {
 	    		      		Main.drawCanvas.sampleList.get(sampleindex+sampletemp).chr = "chr";
 	    		      	  }  	    		      	 
 	    		  	 }  	
  	    	  	 }*/
  	    	  /*	 else {
  	    		  	 File[] listfiles = files[fi].listFiles(new FilenameFilter() {
  	    		  	    public boolean accept(File dir, String name) {
  	    		  	        return name.toLowerCase().endsWith(".list");
  	    		  	    }
  	    		  	 });
  	    		  	if(listfiles.length > 0) { 
  	    		  	
  	    		  	 int index = -1, sampleindex = 0;
  	    		  	
  	    		  	 for(int i = 0; i<listfiles.length; i++) {
  	    		  		 sampleindex = 0;
  	    		  		 index = sampleString.indexOf(listfiles[i].getName().replace(".list", ""));
  	    		  	 	 if (index < 0) continue; 
  	    		  	 	 
  	    		  	 	 for(char letter : sampleString.substring(0, index).toString().toCharArray()) {
  	    		  	 		 if(letter == ';') sampleindex++;
  	    		  	 	 }
  	    		  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile = new File(listfiles[i].getAbsolutePath());	
  	    		  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).resetreadHash();
  	    		  	 }  	
  	    		  	}
  	    	  	 }*/
  	    	  	
  	    //	  	 sampleString = new StringBuffer();
  	    		/*
  	    		if(addDir[f].getName().endsWith(".bam")) {
  	    			if(vcffile == null) {
  	    				bamfile = addDir[f];
  	    			}
  	    			else {
  	    				tabix = true;
  	    	        //	tabixReader = new TabixReader(addDir[f].getCanonicalPath());   
  	    				addSample = new Sample(addDir[f].getName(), (Short)Main.samples, addDir[f].getCanonicalPath());
  	    				addSample.readClass.samFile = addDir[f];
  	    				samFileReader =  SamReaderFactory.makeDefault().open(addDir[f]);			
	     				addSample.readClass.sample = addSample;
	     				if(samFileReader.getFileHeader().getTextHeader().contains("CompleteGenomics")) {							
	     					addSample.readClass.complete = true;
	     				}
	     				break;
  	    			}
  	    		}
  	    		
  	    		if (addDir[f].getName().endsWith(".vcf.gz")) {
  	    			if(bamfile == null) {
  	    				vcffile = addDir[f];
  	    			}
  	    			else {
  	    				tabix = true;
  	    	        	//tabixReader = new TabixReader(addDir[f].getCanonicalPath());   
  	    	        	addSample = new Sample(addDir[f].getName(), (Short)Main.samples, addDir[f].getCanonicalPath());
  	    	        	addSample.readClass.samFile = bamfile;
	     				samFileReader = SamReaderFactory.makeDefault().open(addDir[f]);		
	     				addSample.readClass.sample = addSample;
	     				if(samFileReader.getFileHeader().getTextHeader().contains("CompleteGenomics")) {							
	     					addSample.readClass.complete = true;
	     				}
	     				break;
  	    			}
  	    		}*/
  	    	
  	    }
  	    else {
  	    	if(!files[fi].getName().endsWith(".vcf") && !files[fi].getName().endsWith(".vcf.gz")) {
  	  			continue;
  	  		}
  	    	if(!checkIndex(files[fi])) {
   			 Main.putMessage("Check Tools->View log");
   			 ErrorLog.addError("No index file found for " +files[fi].getName());
   			 Main.drawCanvas.loadingtext = "Creating index for " +files[fi].getName();
   			 if(files[fi].getName().endsWith(".vcf.gz")) {
   				 MethodLibrary.createVCFIndex(files[fi]);
   			 }
   			 else {
   				MethodLibrary.createVCFIndex2(files[fi]);
   			 }   			
  	    	}
  	    	if(fileindex > -1) {
  	    		if(!files[fi].getParent().equals(files[fileindex].getParent())) {
  	    		//	diffPaths = true;
  	    		}
  	    	}
  	    	 fileindex =fi;
  	    	 
        	 addSample = new Sample(files[fi].getName(), (Short)Main.samples, files[fi].getCanonicalPath());        	
        	 Main.drawCanvas.sampleList.add(addSample);        	 
          	 Main.varsamples++;
          	 Main.samples++;
          	 Main.drawCanvas.drawVariables.visiblesamples = Main.samples;
          	 sampleString.append(addSample.getName() +";");
          	 checkMulti(addSample);
          	 Main.drawCanvas.loadingtext = "Loading samples... " +fi +"/" +files.length;
          	 added = true;
          	 VariantHandler.commonSlider.setMaximum(Main.varsamples);
          	 VariantHandler.commonSlider.setUpperValue(Main.varsamples);
          	 VariantHandler.geneSlider.setMaximum(Main.varsamples);
          	
       	  }     	
  	 }
  	 if(!added) {
  		Main.drawCanvas.ready("Loading samples...");
  		 return;
  	 }
  	
  
  //	Main.drawCanvas.drawVariables.visibleend = (short)(Main.drawCanvas.sampleList.size()-1);
  	Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.drawCanvas.sampleList.size());
  	Main.drawCanvas.checkSampleZoom();
  	Main.drawCanvas.resizeCanvas(Main.drawScroll.getViewport().getWidth(),Main.drawScroll.getViewport().getHeight());
  	
  	 int loading = 0;
  	if(!(VariantHandler.hideIndels.isSelected() && VariantHandler.hideSNVs.isSelected())) { 
  	 Main.drawCanvas.loadingtext = "Loading variants...";
  	 for(int i = sampletemp; i < Main.drawCanvas.sampleList.size(); i++) {
  		 linecounter = 0;
  		if(cancelfileread || !Main.drawCanvas.loading) {
			cancelFileRead();
			break;
		}
  		
  		 if(( Main.drawCanvas.sampleList.get(i).getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight < Main.drawScroll.getViewport().getHeight()+Main.drawCanvas.drawVariables.sampleHeight) {
 		//	Main.drawCanvas.drawVariables.visibleend = (short)(Main.drawCanvas.sampleList.get(i).getIndex());		
 		//	Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.drawCanvas.sampleList.get(i).getIndex()+1);	
 			Main.drawCanvas.checkSampleZoom();
 		}
	      Main.drawCanvas.loadbarAll = (int)((loading/(double)(Main.drawCanvas.sampleList.size()-sampletemp))*100);    
	   //   if(!Main.drawCanvas.sampleList.get(i).multipart) {	    	  
	    	  Main.drawCanvas.sampleList.get(i).setMaxCoverage(0F);	
	      try {
	    	 
	   // 	  vcfreader = new VCFFileReader(new File(Main.drawCanvas.sampleList.get(i).getTabixFile()));
	    	/*  try {
	    		  tabixreader = new TabixReader(Main.drawCanvas.sampleList.get(i).getTabixFile());
	    	  }
	    	  catch(Exception ex) {
	    		  JOptionPane.showMessageDialog(Main.chromDraw, "Index file (tbi) not found for " +Main.drawCanvas.sampleList.get(i).getName(), "Error", JOptionPane.ERROR_MESSAGE);
	    		  ErrorLog.addError(ex.getStackTrace());
	    		  ex.printStackTrace();
	    	  }*/
		//      iterator=null;   
	    	  if(Main.drawCanvas.sampleList.get(i).getTabixFile() == null) {
	    			continue;
	    	 	}
	    
			
		  if(start > 10000 && end < Main.drawCanvas.splits.get(0).chromEnd-10000) {	
			 
			  if(Main.drawCanvas.variantsEnd > 0) {
				 // iterator = tabixreader.query(Main.drawCanvas.sampleList.get(i).vcfchr +Main.chromosomeDropdown.getSelectedItem().toString()+":" +Main.drawCanvas.variantsStart +"-" +Main.drawCanvas.variantsEnd);
				  getVariants(Main.chromosomeDropdown.getSelectedItem().toString(), Main.drawCanvas.variantsStart,Main.drawCanvas.variantsEnd, Main.drawCanvas.sampleList.get(i));
					 
			  }
			  else {
				  Main.drawCanvas.variantsStart = start;
				  Main.drawCanvas.variantsEnd = end;
				//  iterator = tabixreader.query(Main.drawCanvas.sampleList.get(i).vcfchr +Main.chromosomeDropdown.getSelectedItem().toString()+":" +Main.drawCanvas.variantsStart +"-" +Main.drawCanvas.variantsEnd);
				  getVariants(Main.chromosomeDropdown.getSelectedItem().toString(), Main.drawCanvas.variantsStart,Main.drawCanvas.variantsEnd, Main.drawCanvas.sampleList.get(i));
					 
			  }
		  }
		  else {
			//  iterator = tabixreader.query(Main.drawCanvas.sampleList.get(i).vcfchr +Main.chromosomeDropdown.getSelectedItem().toString());
			  Main.drawCanvas.variantsStart = 0;
			  Main.drawCanvas.variantsEnd = Main.drawCanvas.splits.get(0).chromEnd;
			  getVariants(Main.chromosomeDropdown.getSelectedItem().toString(), Main.drawCanvas.variantsStart,Main.drawCanvas.variantsEnd, Main.drawCanvas.sampleList.get(i));
				 
		  }		 
	  }
	  catch(Exception e) {
		  Main.showError( e.getMessage(), "Error");
		  ErrorLog.addError(e.getStackTrace());
		  e.printStackTrace();
		  
	  }    
    
     /*
      if(iterator == null) {
    	  tabixreader = null;
    	  continue;
      }*/
	//    }
	/*      else {
	    	  continue;
	      }
      line = null;
     
      current = head;      
      first = true;
    
		while(true) {	      				
			try {
				if(cancelfileread) {
					cancelFileRead();
		  		
		  		}
				try {
					
					line = iterator.next();
				//	vcfline = vcfIterator.next();
			
					if(line == null) {
						break;
					}
					
				//	line = vcfline.getSource();
					
				}
				catch(Exception ex) {
					JOptionPane.showMessageDialog(Main.chromDraw, ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
					Main.cancel();
		//			tabixreader.mFp.close();
					ex.printStackTrace();
					break;
				}
							
				
				split = line.split("\\t+");		
				if(Main.drawCanvas.sampleList.get(i).multiVCF) {
				
					readLineMulti(split, Main.drawCanvas.sampleList.get(i));
				}
				else {
					
					readLine(split, Main.drawCanvas.sampleList.get(i));	
				//	readLine(vcfline, Main.drawCanvas.sampleList.get(i));
				}
				if(first) {
					first = false;
				}
				}
				catch(Exception e) {
					e.printStackTrace();
//					tabixreader.mFp.close();
				}
		}		
//		tabixreader.mFp.close();
		Main.drawCanvas.current = FileRead.head.getNext();
		
		if(Main.drawCanvas.splits.get(0).viewLength < 2000000) {
			
			Main.chromDraw.updateExons = true;
			Main.chromDraw.repaint();
		}
		if(!Main.drawCanvas.loading) {
			Draw.calculateVars = true;
		}
		Draw.updatevars = true;
		 Draw.updateReads = true;
		Main.drawCanvas.repaint();
		loading++;*/
    }
  	
  	 
  	}
	Main.opensamples.setText("Add samples");
	checkSamples();
  	annotate();
  	 if(fileindex > -1) {
		 SearchBamFiles search = new SearchBamFiles(files, fileindex, sampletemp);
		 search.execute();
	 }
    readFiles = false;
 //   Main.drawCanvas.clusterCalc = true;
    if(Control.controlData.controlsOn) {		    	
		Control.applyControl();
	}
    
	
	Main.bedCanvas.intersected = false;
	if(bigcalc) {
		Main.drawCanvas.calcClusters(FileRead.head);
	}
	else {
		Main.drawCanvas.calcClusters(FileRead.head,1);
	}	
	
	
	
	if(Main.bedCanvas.bedOn) {
		Main.drawCanvas.loadingtext = "Annotating variants";
		for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
			if(Main.bedCanvas.bedTrack.get(i).intersect && !Main.bedCanvas.bedTrack.get(i).loading) {
				if(Main.bedCanvas.bedTrack.get(i).small) {
					Main.bedCanvas.annotate(Main.bedCanvas.bedTrack.get(i).getHead());
					
				}
				else {
					BedCanvas.Annotator annotator = Main.bedCanvas.new Annotator(Main.bedCanvas.bedTrack.get(i));
					annotator.annotateVars();
					
				}
			}		
			else if(Main.bedCanvas.bedTrack.get(i).intersect && Main.bedCanvas.bedTrack.get(i).loading) {
				Main.bedCanvas.bedTrack.get(i).waiting = true;
			}
		}
		Main.bedCanvas.intersected = true;
	}	
	
	if(bigcalc) {
		Main.drawCanvas.calcClusters(FileRead.head);
	}
	else {
		Main.drawCanvas.calcClusters(FileRead.head,1);
	}	
    if(Main.drawCanvas.splits.get(0).viewLength < 2000000) {
		
		Main.chromDraw.updateExons = true;
		Main.chromDraw.repaint();
	}
//	Main.drawCanvas.drawVariables.visibleend = Main.samples;
	//Main.drawCanvas.drawVariables.visiblesamples = Main.samples;
	Main.drawCanvas.checkSampleZoom();
	Main.drawCanvas.current = head;
  //	Draw.updatevars = true;
//	Main.drawCanvas.repaint();	
	current = null;
	Main.drawCanvas.ready("Loading samples...");
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
	
	
	 	
	}
	catch(Exception e)  {
		e.printStackTrace();
		ErrorLog.addError(e.getStackTrace());
	}	
		
	}
	
static void checkSamples() {
	/*if(Main.varsamples == 0) {
		Main.drawCanvas.varbuffer = MethodLibrary.toCompatibleImage(new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB));	
		Main.drawCanvas.g2v = (Graphics2D)Main.drawCanvas.varbuffer.getGraphics();
	}
	else {*/
//		Main.drawCanvas.varbuffer = MethodLibrary.toCompatibleImage(new BufferedImage(Main.screenSize.width, Main.screenSize.height, BufferedImage.TYPE_INT_ARGB));	
//		Main.drawCanvas.g2v = (Graphics2D)Main.drawCanvas.varbuffer.getGraphics();
//	}
	if(Main.varsamples < 2) {
		VariantHandler.commonSlider.setMaximum(1);
		VariantHandler.commonSlider.setValue(1);
		VariantHandler.commonSlider.setUpperValue(1);
		VariantHandler.geneSlider.setMaximum(1);
		VariantHandler.geneSlider.setValue(1);
		VariantHandler.filterPanes.setEnabledAt(VariantHandler.filterPanes.getTabCount()-1, false);
		VariantHandler.filterPanes.setToolTipTextAt(VariantHandler.filterPanes.getTabCount()-1, "Open more samples to compare variants.");
		
	}
	else {
		VariantHandler.commonSlider.setMaximum(Main.varsamples);
		VariantHandler.commonSlider.setValue(1);
		VariantHandler.commonSlider.setUpperValue(Main.varsamples);
		VariantHandler.geneSlider.setMaximum(Main.varsamples);
		VariantHandler.geneSlider.setValue(1);
		VariantHandler.filterPanes.setEnabledAt(VariantHandler.filterPanes.getTabCount()-1, true);
		VariantHandler.filterPanes.setToolTipTextAt(VariantHandler.filterPanes.getTabCount()-1, "Compare variants.");
	}
	
	if(Main.readsamples > 0) {
	/*	Main.drawCanvas.readbuffer = MethodLibrary.toCompatibleImage(new BufferedImage(Main.screenSize.width, Main.screenSize.height, BufferedImage.TYPE_INT_ARGB));	
		Main.drawCanvas.rbuf = (Graphics2D)Main.drawCanvas.readbuffer.getGraphics();
		Main.drawCanvas.backupr = Main.drawCanvas.rbuf.getComposite();
		Main.drawCanvas.rbuf.setRenderingHints(Draw.rh);
		Main.drawCanvas.coveragebuffer = MethodLibrary.toCompatibleImage(new BufferedImage(Main.screenSize.width, Main.screenSize.height, BufferedImage.TYPE_INT_ARGB));	
		Main.drawCanvas.cbuf = (Graphics2D)Main.drawCanvas.coveragebuffer.getGraphics();
		Main.drawCanvas.backupc = Main.drawCanvas.cbuf.getComposite();
		*/
		Main.average.setEnabled(true);
		
		Main.variantCaller.setEnabled(true);
		
		
	}
	else {
/*		Main.drawCanvas.readbuffer = MethodLibrary.toCompatibleImage(new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB));	
		Main.drawCanvas.rbuf = (Graphics2D)Main.drawCanvas.readbuffer.getGraphics();
		Main.drawCanvas.backupr = Main.drawCanvas.rbuf.getComposite();
		Main.drawCanvas.rbuf.setRenderingHints(Draw.rh);
		Main.drawCanvas.coveragebuffer = MethodLibrary.toCompatibleImage(new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB));	
		Main.drawCanvas.cbuf = (Graphics2D)Main.drawCanvas.coveragebuffer.getGraphics();
		Main.drawCanvas.backupc = Main.drawCanvas.cbuf.getComposite();
		*/
		Main.average.setEnabled(false);
		Main.average.setToolTipText("No bam/cram files opened");
		Main.variantCaller.setEnabled(false);
		Main.variantCaller.setToolTipText("No bam/cram files opened");
	}
	
}

static void annotate() {
		
		if(Main.drawCanvas.splits.get(0).getGenes().size() == 0) {
			return;
		}
		Transcript transcript;
		Gene gene, prevGene = Main.drawCanvas.splits.get(0).getGenes().get(0);
		VarNode current = FileRead.head.getNext();
		Transcript.Exon exon;
		int position = 0, baselength;
		
		boolean intronic = true;
		try {
			
		if(current != null) {
			
			for(int g=0; g<Main.drawCanvas.splits.get(0).getGenes().size(); g++ ) {
				gene = Main.drawCanvas.splits.get(0).getGenes().get(g);				
				
			/*	while(current != null && current.getPosition() < gene.getStart()) {
					if(!current.isInGene()) {
						if(current.getTranscripts() == null) {
							current.setTranscripts();
						}					
						current.getTranscripts().add(prevGene.getTranscripts().get(0));
						current.getTranscripts().add(gene.getTranscripts().get(0));
					}
					current = current.getNext();
				}*/
				
				if(current == null) {
					break;
				}
				for(int t = 0; t<gene.getTranscripts().size(); t++) {
					transcript = gene.getTranscripts().get(t);
				
					if(current != null && current.getPrev() != null) {
						while(current.getPrev().getPosition() >= transcript.getStart()) {
							if(current.getPrev() != null) {
								current = current.getPrev();
							}						
						}
					}
					
					position = current.getPosition();					
					
					if(current.indel) {
						position++;
						baselength = MethodLibrary.getBaseLength(current.vars, 1);		
					}
					while(position < transcript.getEnd()) {
						try {
						if(position >= transcript.getStart() && position <= transcript.getEnd()) {
							
							current.setInGene();					
							
							baselength = 0;						
							intronic = true;
							
							for(int e=0;e<transcript.getExons().length; e++) {
								
								exon = transcript.getExons()[e];							
								
								if(position+baselength >= exon.getStart()-2 && position < exon.getEnd()+1) {
									if(current.getExons() == null) {
										current.setExons();
									}
									intronic = false;	
									if(!current.getExons().contains(exon)) {
																	
										current.getExons().add(exon);
										if(exon.getStartPhase() > -1 && position+baselength >= exon.getTranscript().getCodingStart() && position < exon.getTranscript().getCodingEnd() ) {
											current.coding = true;
										}
										break;
									}								
									
								}
													
							}				
							if(intronic) {
								if(current.getTranscripts() == null) {
									current.setTranscripts();
								}		
							
									current.getTranscripts().add(transcript);
						
								}
						}					
						if(!current.isInGene()) {
							if(current.getTranscripts() == null) {
								current.setTranscripts();
								current.getTranscripts().add(prevGene.getTranscripts().get(0));
								current.getTranscripts().add(gene.getTranscripts().get(0));
								
							}							
						}
						
						if(current.getNext() != null) {
							current = current.getNext();	
							position = current.getPosition();						
						}
						else {
							
							break;
						}		
					}
					catch(Exception e) {
						System.out.println(position);
						e.printStackTrace();
						break;
					}
					}			
				}	
				
				if(gene.getEnd() > prevGene.getEnd()) {
					prevGene = gene;
				}
				
			}						
			while(current != null) {
				if(!current.isInGene()) {
					if(current.getTranscripts() == null) {
						current.setTranscripts();
						current.getTranscripts().add(prevGene.getTranscripts().get(0));
						current.getTranscripts().add(prevGene.getTranscripts().get(0));
					}					
					
				}
				current = current.getNext();
			}
		}
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
			
		}
		current = null;
		transcript = null;
		exon = null;
	}
	public void readLineMulti(String[] split, Sample sample) {
		
		samplecount = 0;
		try {
			pos = Integer.parseInt(split[1])-1;		
		}
		catch(Exception e) {
			return;
		}
		if(pos < Main.drawCanvas.variantsStart) {
			return;
		}
		else if(pos >=Main.drawCanvas.variantsEnd) {
			stop = true;
			return;
		}
		/*else if(sample.getVCFInput().getFilePointer() > sample.vcfEndPos) {
			
			return;
		}*/
		if(linecounter == 0 || pos -linecounter > (Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart)/100) {		
			Main.drawCanvas.loadBarSample = (int)(((pos-Main.drawCanvas.variantsStart)/(double)(Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart))*100);     
			Draw.updatevars = true;
			linecounter = pos; 				
		}
		
		first = true;
		
		while(current != null && current.getNext() != null && current.getPosition() < pos ){			
			current = current.getNext();		      				
		}
		if(current.getPosition() == pos) {
			first = false;
			
		}
		
		for(int s = 9; s < split.length;s++) {
			try {	
				currentSample = Main.drawCanvas.sampleList.get(sample.getIndex()+1+s-9 );		
				if(currentSample.removed) {
					
					continue;
				}
				
				info = split[s].split(":");								
				noref = false;
				
				if(split[8].contains("GT")) {	
					gtindex = split[8].indexOf("GT")/3;
					if(info[gtindex].contains(".") || info[gtindex].length() < 3) {
						continue;
					}
					firstallele = Short.parseShort(""+info[gtindex].split("/")[0]);
					secondallele = Short.parseShort(""+info[gtindex].split("/")[1]);			
					genotype = firstallele == secondallele;
					
					if(genotype && firstallele == 0) {
						continue;
					}
				
					if(split[8].contains("AD:") || split[8].contains("AD ")) {
						if(split[8].contains("RD")) {				
							refcalls = (short)(Short.parseShort(info[split[8].indexOf("RD")/3]));
							altcalls = Short.parseShort(info[split[8].indexOf("AD")/3]);
						}
						else {				
							coverages = info[split[8].indexOf("AD")/3].split(",");
							calls1 = (short)(Short.parseShort(coverages[firstallele]));
							calls2 = Short.parseShort(coverages[secondallele]);				
						}			
					}
					else {			
						calls1 = (short)20;
						calls2 = (short)20;
					}
					
					if(!genotype) {
						if(firstallele  == 0) {
							refallele = firstallele;
							refcalls = calls1;
							altallele = secondallele;
							altcalls = calls2;
						}
						else if(secondallele  == 0){
							refallele = secondallele;
							refcalls = calls2;
							altallele = firstallele;					
							altcalls = calls1;
						}
						else {
							noref = true;
							refallele = firstallele;
							refcalls = calls1;
							altallele = secondallele;
							altcalls = calls2;
						}
					}
					else {
						
						refcalls = calls1; //(short)(Short.parseShort(coverages[0]));
						altallele = secondallele;
						altcalls = calls2;
						
					}
				}
				
				if(!split[4].contains(",")) {
					altbase = getVariant(split[3], split[4]);
					
				}	
				else if(!noref){
					
					altbase = getVariant(split[3], split[4].split(",")[altallele-1]);				
					
				}
				else {
					altbase = getVariant(split[3], split[4].split(",")[altallele-1]);	
					altbase2 = getVariant(split[3], split[4].split(",")[refallele-1]);
					
				}
				if(altbase.contains("*") || (altbase2!=null && altbase2.contains("*"))) {
					continue;
				}
				
				if(split[8].contains("Q")) {
					if(split[8].contains("GQ") && !info[split[8].indexOf("GQ")/3].equals(".")) {						
						gq = (float)Double.parseDouble(info[split[8].indexOf("GQ")/3]);						
					}
					else if(split[8].contains(":BQ:") && !info[split[8].indexOf("BQ")/3].equals(".")) {						
						quality = (float)Double.parseDouble(info[split[8].indexOf("BQ")/3]);						
					}					
					else if(split[7].contains("SSC")) {
						quality = Float.parseFloat(split[7].substring(split[7].indexOf("SSC") +4).split(";")[0]);							
					}
					
				}
				if(quality == null) {
					if(split[5].matches("\\d+.?.*")) {
						quality = (float)Double.parseDouble(split[5]);
						
					}				
				}
				
				HashMap<String, Float> advancedQualities = null;				
				
				if(VariantHandler.freeze.isSelected()) {
					
					if(refcalls+altcalls < VariantHandler.coverageSlider.getValue()) {
						
						continue;
					}
					
					if(quality != null &&  quality < VariantHandler.qualitySlider.getValue()) {
						continue;
					}
					if(gq != null && gq < VariantHandler.gqSlider.getValue()) {
						continue;
					}
					if(altcalls/(double)(refcalls+altcalls) < VariantHandler.callSlider.getValue()/100.0) {
						
						continue;
					}
					
					if(VariantHandler.hideSNVs.isSelected() && altbase.length() == 1 ) {
						continue;
					}
					
					
					if(VariantHandler.hideIndels.isSelected() && altbase.length() > 1) {
						continue;
					}
					
					if(VariantHandler.rscode.isSelected() && !split[2].equals(".")) {
						continue;
					}	
					
				}				
				if(split.length > 6) {
					if(checkAdvFilters(split[6], altbase.length() > 1)) {
						return;
					}	
				}
				if(split.length > 7) {	
					if(checkAdvQuals(split[7], altbase.length() > 1)) {
						return;
					}				
				}
				if(refcalls+altcalls > VariantHandler.maxCoverageSlider.getMaximum()) {
					VariantHandler.maxCoverageSlider.setMaximum(refcalls+altcalls);
					VariantHandler.maxCoverageSlider.setValue(refcalls+altcalls);
				}
			/*	if(currentSample.getMaxCoverage() < calls1+calls2) {
					currentSample.setMaxCoverage((float)(calls1+calls2));	
				}
				*/
					if(first && current.getNext() == null) {							
							
							if(!split[2].equals(".")) {
								current.putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, gq, advancedQualities, split[2], currentSample, current, null));		
							}
							else {
								current.putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, gq,advancedQualities,  null, currentSample, current, null));		
							}
							if(noref ) {
								if(split.length > 6) {
									if(checkAdvFilters(split[6], altbase2.length() > 1)) {
										return;
									}	
								}
								if(split.length > 7) {	
									if(checkAdvQuals(split[7], altbase2.length() > 1)) {
										return;
									}				
								}
								
								current.getNext().addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, gq, advancedQualities, currentSample);
							}						
							
							current = current.getNext();	
							first = false;
						
					}						
					else if(pos == current.getPosition()){		
						
							current.addSample(altbase, refcalls+altcalls, altcalls, genotype, quality, gq, advancedQualities, currentSample);
							if(noref ) {
								if(split.length > 6) {
									if(checkAdvFilters(split[6], altbase2.length() > 1)) {
										return;
									}	
								}
								if(split.length > 7) {	
									if(checkAdvQuals(split[7], altbase2.length() > 1)) {
										return;
									}				
								}
								
								current.addSample(altbase2, refcalls+altcalls, refcalls, genotype, quality, gq, advancedQualities, currentSample);
							}
							if(current.isRscode() == null && !split[2].equals(".")) {							
								current.setRscode(split[2]);
							}					
							
							
					}
				
					else if(current.getPosition() > pos) {
						
							if(!split[2].equals(".")) {
								current.getPrev().putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, refcalls+altcalls, altcalls, genotype, quality, gq, advancedQualities, split[2], currentSample, current.getPrev(), current));
							}
							else {
								current.getPrev().putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, refcalls+altcalls, altcalls, genotype, quality, gq,advancedQualities,  null, currentSample, current.getPrev(), current));
								
							}
							if(noref ) {
								if(split.length > 7) {	
									if(split.length > 6) {
										if(checkAdvFilters(split[6], altbase2.length() > 1)) {
											return;
										}	
									}
									if(checkAdvQuals(split[7], altbase2.length() > 1)) {
										return;
									}				
								}
								
								current.getPrev().addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, gq,advancedQualities,  currentSample);
							}
							current.putPrev(current.getPrev().getNext());
							current = current.getPrev();
						
					}
					
				}
				catch(Exception ex) {
					ErrorLog.addError(ex.getStackTrace());
					ex.printStackTrace();
					for(int i=0; i<split.length; i++) {
						System.out.print(split[i] +" ");
					}
					System.out.println();
					break;
				}		
		
			if(first) {
				first = false;
			}
		}		
		
	}
	
	
	boolean checkAdvQuals(String split, boolean indel) {
		if(!VariantHandler.indelFilters.isSelected()) {	
			if(Main.drawCanvas.drawVariables.advQDraw != null && Main.drawCanvas.drawVariables.advQDraw.size() > 0) {
				
				for(int i = 0 ; i<Main.drawCanvas.drawVariables.advQDraw.size(); i++) {		
					
					if(split.contains(Main.drawCanvas.drawVariables.advQDraw.get(i).key)) {	
						if(Main.drawCanvas.drawVariables.advQDraw.get(i).format.equals("<")) {
							if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDraw.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDraw.get(i).key.length()+1).split(";")[0].trim()) < Main.drawCanvas.drawVariables.advQDraw.get(i).value) {
								return true;
							}
						}
						else if(Main.drawCanvas.drawVariables.advQDraw.get(i).format.equals("<=")) {
							if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDraw.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDraw.get(i).key.length()+1).split(";")[0].trim()) <= Main.drawCanvas.drawVariables.advQDraw.get(i).value) {
								return true;
							}
						}
						else if(Main.drawCanvas.drawVariables.advQDraw.get(i).format.equals(">")) {
							if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDraw.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDraw.get(i).key.length()+1).split(";")[0].trim()) > Main.drawCanvas.drawVariables.advQDraw.get(i).value) {
								return true;
							}
						}
						else {
							if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDraw.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDraw.get(i).key.length()+1).split(";")[0].trim()) >= Main.drawCanvas.drawVariables.advQDraw.get(i).value) {
								return true;
							}
						}
					}					
				}					
			}
		}
		else {
			if(indel) {
				if(Main.drawCanvas.drawVariables.advQDrawIndel != null && Main.drawCanvas.drawVariables.advQDrawIndel.size() > 0) {
					
					for(int i = 0 ; i<Main.drawCanvas.drawVariables.advQDrawIndel.size(); i++) {		
						
						if(split.contains(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key)) {	
							if(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).format.equals("<")) {
								if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key.length()+1).split(";")[0].trim()) < Main.drawCanvas.drawVariables.advQDrawIndel.get(i).value) {
									return true;
								}
							}
							else if(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).format.equals("<=")) {
								if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key.length()+1).split(";")[0].trim()) <= Main.drawCanvas.drawVariables.advQDrawIndel.get(i).value) {
									return true;
								}
							}
							else if(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).format.equals(">")) {
								if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key.length()+1).split(";")[0].trim()) > Main.drawCanvas.drawVariables.advQDrawIndel.get(i).value) {
									return true;
								}
							}
							else {
								if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key.length()+1).split(";")[0].trim()) >= Main.drawCanvas.drawVariables.advQDrawIndel.get(i).value) {
									return true;
								}
							}
						}					
					}					
				}
			}
			else {
				if(Main.drawCanvas.drawVariables.advQDraw != null && Main.drawCanvas.drawVariables.advQDraw.size() > 0) {
					
					for(int i = 0 ; i<Main.drawCanvas.drawVariables.advQDraw.size(); i++) {		
						
						if(split.contains(Main.drawCanvas.drawVariables.advQDraw.get(i).key)) {	
							if(Main.drawCanvas.drawVariables.advQDraw.get(i).format.equals("<")) {
								if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDraw.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDraw.get(i).key.length()+1).split(";")[0].trim()) < Main.drawCanvas.drawVariables.advQDraw.get(i).value) {
									return true;
								}
							}
							else if(Main.drawCanvas.drawVariables.advQDraw.get(i).format.equals("<=")) {
								if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDraw.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDraw.get(i).key.length()+1).split(";")[0].trim()) <= Main.drawCanvas.drawVariables.advQDraw.get(i).value) {
									return true;
								}
							}
							else if(Main.drawCanvas.drawVariables.advQDraw.get(i).format.equals(">")) {
								if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDraw.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDraw.get(i).key.length()+1).split(";")[0].trim()) > Main.drawCanvas.drawVariables.advQDraw.get(i).value) {
									return true;
								}
							}
							else {
								if(Float.parseFloat(split.substring(split.indexOf(Main.drawCanvas.drawVariables.advQDraw.get(i).key+"=")+Main.drawCanvas.drawVariables.advQDraw.get(i).key.length()+1).split(";")[0].trim()) >= Main.drawCanvas.drawVariables.advQDraw.get(i).value) {
									return true;
								}
							}
						}					
					}					
				}
			}
		}
		
		return false;
	}
	boolean checkAdvFilters(String split, boolean indel) {
		if(!VariantHandler.indelFilters.isSelected()) {	
			if(Main.drawCanvas.drawVariables.advQDraw != null && Main.drawCanvas.drawVariables.advQDraw.size() > 0) {
				
				for(int i = 0 ; i<Main.drawCanvas.drawVariables.advQDraw.size(); i++) {		
					
					if(split.equals(Main.drawCanvas.drawVariables.advQDraw.get(i).key)) {	
						return true;
					}					
				}					
			}
		}
		else {
			if(indel) {
				if(Main.drawCanvas.drawVariables.advQDrawIndel != null && Main.drawCanvas.drawVariables.advQDrawIndel.size() > 0) {
					
					for(int i = 0 ; i<Main.drawCanvas.drawVariables.advQDrawIndel.size(); i++) {		
						
						if(split.equals(Main.drawCanvas.drawVariables.advQDrawIndel.get(i).key)) {	
							
							return true;
						}					
					}					
				}
			}
			else {
				if(Main.drawCanvas.drawVariables.advQDraw != null && Main.drawCanvas.drawVariables.advQDraw.size() > 0) {
					
					for(int i = 0 ; i<Main.drawCanvas.drawVariables.advQDraw.size(); i++) {			
						
						if(split.equals(Main.drawCanvas.drawVariables.advQDraw.get(i).key)) {	
							
							return true;
						}					
					}					
				}
			}
		}
		
		return false;
	}
	public void readLine(String[] split, Sample sample) {		
		try {
			
		if(split.length < 3) {			
			return;
		}	
		
		if(split[0].startsWith("#") || split[4].equals("*")) {
			return;
		}
		pos = Integer.parseInt(split[1])-1;	
		
		if(pos < Main.drawCanvas.variantsStart) {
			return;
		}
		else if(pos >= Main.drawCanvas.variantsEnd) {
			stop = true;
			return;
		}
		if(linecounter  == 0 || pos -linecounter > (Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart)/100) {
			
			 Main.drawCanvas.loadBarSample = (int)(((pos-Main.drawCanvas.variantsStart)/(double)(Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart))*100);     
			// Draw.updatevars = true;
			 linecounter = pos; 
			 if(search) {
				if(!Main.drawCanvas.loading) {
					Draw.calculateVars = true;
				}
			 }
			
	}

			
	/*	if(first) {
			
			info = split[split.length-1].split(":");
			
			if(info[0].length() == 1 || info[0].charAt(0) == '0' && info[0].charAt(2) =='0') {
				info = split[split.length-2].split(":");
				sample.infolocation = (short)(split.length-2);
			}
			else {
				sample.infolocation = (short)(split.length-1);
			}
		}
		else {
		
		*/	
		//}
		if(sample.somaticColumn != null) {
			if(sample.somaticColumn > split.length-1) {
				sample.somaticColumn = null;
				info = split[split.length-1].split(":");
			}
			else {
				info = split[sample.somaticColumn].split(":");
			}		
		}
		else {
			info = split[split.length-1].split(":");
		}
		noref = false;
		multi = false;
		if(split.length > 8 && split[8].contains("GT")) {	
			gtindex = split[8].indexOf("GT")/3;
			if(info[gtindex].contains(".")) {
				return;
			}
			try {
				if(info[gtindex].contains("|")) {
					
					firstallele = Short.parseShort(""+info[gtindex].split("|")[0]);
					secondallele = Short.parseShort(""+info[gtindex].split("|")[2]);	
				}
				else {
					firstallele = Short.parseShort(""+info[gtindex].split("/")[0]);
					secondallele = Short.parseShort(""+info[gtindex].split("/")[1]);	
				}
			}
			catch(Exception ex) {
				System.out.println(info[gtindex] +" " +info[gtindex].split("|").length);
				
			}
			genotype = firstallele == secondallele;
			
			if(genotype && firstallele == 0) {
				return;
			}
			
			
			if(split[8].contains("AD:") || split[8].contains("AD ")) {
				
				if(split[8].contains("RD")) {		
					//TODO jaettuna kolmella ei välttämättä pidä paikkaansa, jos format kentässä on esim: GT:GQ:FREQ:AD:RD....
					calls1 = Short.parseShort(info[split[8].indexOf("RD")/3]);
					if(info[split[8].indexOf("AD")/3].contains(",")) {
						calls2 = Short.parseShort(info[split[8].indexOf("AD")/3].split(",")[1]);
					}
					else {
						calls2 = Short.parseShort(info[split[8].indexOf("AD")/3]);
					}			
				}
				else {				
					coverages = info[split[8].indexOf("AD")/3].split(",");
					if(genotype) {
						calls1 = (Short.parseShort(coverages[0]));
					}
					else {
						calls1 = (Short.parseShort(coverages[firstallele]));
					}
						calls2 = Short.parseShort(coverages[secondallele]);					
				}	
				if(!genotype) {
					if(firstallele  == 0) {
						refallele = firstallele;
						refcalls = calls1;
						altallele = secondallele;
						altcalls = calls2;
					}
					else if(secondallele  == 0){
						refallele = secondallele;
						refcalls = calls2;
						altallele = firstallele;					
						altcalls = calls1;
					}
					else {					
						noref = true;
						refallele = firstallele;
						refcalls = calls1;
						altallele = secondallele;
						altcalls = calls2;
					}
				}
				else {
					try {
						refcalls =calls1; // Short.parseShort(info[split[8].indexOf("AD")/3].split(",")[0]);
						altallele = secondallele;
						altcalls = calls2;
					}
					catch(Exception e) {
						for(int i = 0 ; i<split.length; i++) {
							System.out.print(split[i] +"\t");
						}
						System.out.println();
						e.printStackTrace();
					}
				}
			}
			else if (split[7].contains("DP4")) {
				coverages = split[7].substring(split[7].indexOf("DP4")+4).split(";")[0].split(",");
				refallele = firstallele;					
				altallele = secondallele;	
				refcalls = (short)(Short.parseShort(coverages[0])+Short.parseShort(coverages[1]));
				altcalls = (short)(Short.parseShort(coverages[2])+Short.parseShort(coverages[3]));
			}
			else {		
				if(firstallele  == 0) {
					refallele = firstallele;					
					altallele = secondallele;					
				}
				else if(secondallele  == 0){
					refallele = secondallele;					
					altallele = firstallele;					
					
				}
				else {					
					
					refallele = firstallele;					
					altallele = secondallele;					
				}
				refcalls = (short)20;
				altcalls = (short)20;
			}			
		}
		else {
			refcalls = (short)20;
			altcalls = (short)20;
		}
		
		if(!split[4].contains(",")) {
			altbase = getVariant(split[3], split[4]);
			
		}	
		else if(!noref){			
			if (altallele > 0) {
				altbase = getVariant(split[3], split[4].split(",")[altallele-1]);				
			}
			else {
				altbase = getVariant(split[3], split[4].split(",")[0]);		
				multi = true;
			}
		}
		else {
			altbase = getVariant(split[3], split[4].split(",")[altallele-1]);		
			altbase2 = getVariant(split[3], split[4].split(",")[refallele-1]);			
			
		}
		quality = null;
		
		if(split.length > 8 && split[8].contains("Q")) {
			
			if(split[8].contains("GQ") && !info[split[8].indexOf("GQ")/3].equals(".")) {
				
				gq = (float)Double.parseDouble(info[split[8].indexOf("GQ")/3]);
				
			}
			else if(split[8].contains(":BQ:") && !info[split[8].indexOf("BQ")/3].equals(".")) {
				
				quality = (float)Double.parseDouble(info[split[8].indexOf("BQ")/3]);
				
			}
			
			else if(split[7].contains("SSC")) {
				quality = Float.parseFloat(split[7].substring(split[7].indexOf("SSC") +4).split(";")[0]);				
					
			}
			
		}
		if(quality == null) {
			if(split.length > 5 && split[5].matches("\\d+.?.*")) {
				quality = (float)Double.parseDouble(split[5]);			
				
			}			
		}
		HashMap<String, Float> advancedQualities = null;
		
	
		if(VariantHandler.freeze.isSelected()) {
			if(refcalls+altcalls < VariantHandler.coverageSlider.getValue()) {				
				return;
			}
			if(quality < VariantHandler.qualitySlider.getValue()) {				
				return;
			}
			if(gq != null && gq < VariantHandler.gqSlider.getValue()) {
				return;
			}
			if(altcalls/(double)(refcalls+altcalls) < VariantHandler.callSlider.getValue()/100.0) {
				return;
			}
			
			if(VariantHandler.hideSNVs.isSelected() && altbase.length() == 1) {
				return;
			}
			
			if(VariantHandler.hideIndels.isSelected() && altbase.length() > 1) {
				return;
			}
			
			if(VariantHandler.rscode.isSelected() && !split[2].equals(".")) {
				return;
			}			
		}		
		if(split.length > 6) {
			if(checkAdvFilters(split[6], altbase.length() > 1)) {
				return;
			}	
		}
		if(split.length > 7) {	
			if(checkAdvQuals(split[7], altbase.length() > 1)) {
				return;
			}				
		}
		if(refcalls+altcalls > VariantHandler.maxCoverageSlider.getMaximum()) {
			VariantHandler.maxCoverageSlider.setMaximum(refcalls+altcalls);
			VariantHandler.maxCoverageSlider.setValue(refcalls+altcalls);
		}
			
		
		while(current.getNext() != null && current != null && current.getPosition() < pos ){	
		
			current = current.getNext();		      				
		}		
		
		if(current.getNext() == null && current.getPosition() < pos) {	
			
			try {			
				
				if(!split[2].equals(".")) {
					current.putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, gq, advancedQualities, split[2], sample, current, null));		
				}
				else {
					
					current.putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, gq, advancedQualities, null, sample, current, null));		
					
				}
				
				if(noref ) {	
					if(split.length > 6) {
						if(checkAdvFilters(split[6], altbase.length() > 1)) {
							return;
						}	
					}
					if(split.length > 7) {	
						if(checkAdvQuals(split[7], altbase2.length() > 1)) {
							return;
						}				
					}
					
					current.getNext().addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, gq, advancedQualities, sample);				
				}
				
				if(multi) {
					String[] altsplit = split[4].split(",");
					
					for(int i = 1; i<altsplit.length; i++) {
						altbase = getVariant(split[3],altsplit[i]);
						if(split.length > 6) {
							if(checkAdvFilters(split[6], altbase.length() > 1)) {
								continue;
							}	
						}
						if(split.length > 7) {	
							if(checkAdvQuals(split[7], altbase.length() > 1)) {
								continue;
							}				
						}
						
						current.getNext().addSample(altbase, (short)(refcalls+altcalls), refcalls, genotype, quality, gq,advancedQualities,  sample);		
						
					}
				}
			}
			catch(Exception e) {
				ErrorLog.addError(e.getStackTrace());
				e.printStackTrace();
				Main.drawCanvas.ready("all");
			}			
		}		
		else if(current.getPosition() == pos) {
				
			current.addSample(altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, gq,advancedQualities,  sample);
			
			if(noref ) {
				if(split.length > 6) {
					if(checkAdvFilters(split[6], altbase2.length() > 1)) {
						return;
					}	
				}
				if(split.length > 7) {	
					if(checkAdvQuals(split[7], altbase2.length() > 1)) {
						return;
					}				
				}
				
				current.addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, gq, advancedQualities, sample);
			}
			if(current.isRscode() == null && !split[2].equals(".")) {							
				current.setRscode(split[2]);
			}
			if(multi) {
				
				String[] altsplit = split[4].split(",");
				for(int i = 1; i<altsplit.length; i++) {
					altbase = getVariant(split[3],altsplit[i]);
					if(split.length > 6) {
						if(checkAdvFilters(split[6], altbase.length() > 1)) {
							continue;
						}	
					}
					if(split.length > 7) {	
						if(checkAdvQuals(split[7], altbase.length() > 1)) {
							continue;
						}				
					}
					
					current.getNext().addSample(altbase, (short)(refcalls+altcalls), refcalls, genotype, quality, gq, advancedQualities, sample);		
					
				}
			}
		}
		else if(current.getPosition() > pos) {
			
				if(!split[2].equals(".")) {
					current.getPrev().putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, gq, advancedQualities, split[2], sample, current.getPrev(), current));
				}
				else {
					if(current.getPrev() == null) {
						System.out.println(current.getPosition());
						Main.cancel();;
					}
					current.getPrev().putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, gq,advancedQualities,  null, sample, current.getPrev(), current));
					
				}
				
				current.putPrev(current.getPrev().getNext());
				current = current.getPrev();
				if(noref ) {
					if(split.length > 6) {
						if(checkAdvFilters(split[6], altbase2.length() > 1)) {
							return;
						}	
					}
					if(split.length > 7) {	
						if(checkAdvQuals(split[7], altbase2.length() > 1)) {
							return;
						}				
					}
					
					current.addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, gq, advancedQualities, sample);
				}
				if(multi) {
					
					String[] altsplit = split[4].split(",");
					for(int i = 1; i<altsplit.length; i++) {
						altbase = getVariant(split[3],altsplit[i]);
						if(split.length > 6) {
							if(checkAdvFilters(split[6], altbase.length() > 1)) {
								continue;
							}	
						}
						if(split.length > 7) {	
							if(checkAdvQuals(split[7], altbase.length() > 1)) {
								continue;
							}				
						}
						
						current.getNext().addSample(altbase, (short)(refcalls+altcalls), refcalls, genotype, quality, gq, advancedQualities, sample);		
						
					}
				}
		}
		
		}
		catch(Exception ex) {
			//System.out.println(split[8] +" " +split[10] +" " +split[3] +" " +split[4]);
			ErrorLog.addError(ex.getStackTrace());
			ex.printStackTrace();
			
		}
		
	}
	
	File getListFile(File listfile, String chrom) {
		
		try {
			if(!listfile.exists()) {
				return new File(listfile.getCanonicalPath().replace(".list", ".cram"));
			}
		 BufferedReader reader = new BufferedReader(new FileReader(listfile));
		 String mnt = "/mnt";
	  	 if(!listfile.getAbsolutePath().startsWith("/mnt")) {
	  		if(listfile.getAbsolutePath().indexOf("\\cg") < 0) {	  	
	  			if(listfile.getAbsolutePath().indexOf("/cg") > -1) {
	  				mnt = listfile.getAbsolutePath().substring(0,listfile.getAbsolutePath().indexOf("/cg"));
	  			}
	  			else {
	  				mnt = "X:";
	  			}
	  		}
	  		else {
	  			mnt = listfile.getAbsolutePath().substring(0,listfile.getAbsolutePath().indexOf("\\cg"));
	  		}
	  	 }
		 String listline;			
		
  	 	 while((listline = reader.readLine()) != null) {
  	 		if(Main.selectedChrom > 24) {
  	 			 if(listline.contains("_" +"GL" +".")) {
  	 			
  	  	 			reader.close();  	 		
  	  	 			return new File(listline.replace("/mnt", mnt));  	 
  	  	 		 }
  	 		 }
  	 		 else if(listline.contains("_" +chrom +".")) {
  	 			reader.close();  	 		
  	 			return new File(listline.replace("/mnt", mnt));  	 			
  	 		 }  	 		
  	 	 }	 
  	 	reader.close();	
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();			
		}		
  	 	return null;
	}
	
	static String getVariant(String ref, String alt) {
		
		
		if(ref.length() == 1 && alt.length() == 1) {
			return alt;
		}
		else if (ref.length() == 1 || alt.length() == 1) {
			if(ref.length() == 2) {				
				return "del" +ref.substring(1);				
			}
			else if(alt.length() >= 2) {
					return "ins" +alt.substring(1);					
			}
			else {				
					return "del" +(ref.length()-1);						
			}
			
		}
		else if(alt.length() == ref.length()) {
			
			if(ref.contains("-")) {
				
				return getVariant(ref, ""+alt.charAt(0));
			}
			else if(alt.contains("-")) {
				return getVariant(ref, ""+alt.replaceAll("-", ""));
			}
			else {
				return "" +alt.charAt(0);
			}
		}
		else {		
		
			if(ref.length() < alt.length()) {
				
					return "ins" +alt.substring(1, alt.length()-(ref.length()-1));
				
				
			}
			else {
				
				if(ref.length() -alt.length() == 1) {
					return "del" +ref.substring(1, ref.length()-(alt.length()-1));
				}
				else {
					return "del" +(ref.length() -alt.length());
				}
				
			}
			
		}
				
	
	}
	void addToSequence(SplitClass split, StringBuffer buffer, int length) {
		
	/*	if(length < 0) {
			
			buffer.insert(0, Main.chromDraw.getSeq(split.chrom, (int)(split.readSeqStart-length)-1, split.readSeqStart, Main.referenceFile));
			split.readSeqStart = (int)(split.readSeqStart-length)-1; 
			
			if(split.readSeqStart < 0) {
				split.readSeqStart = 0;
			}				
			readSeqStart = split.readSeqStart;
		}
		else {
			split.readSequence.append(Main.chromDraw.getSeq(split.chrom, split.readSeqStart+split.readSequence.length(), (int)(split.readSeqStart+split.readSequence.length()+length+200), Main.referenceFile));
			
		}	*/		
	}

	Iterator<SAMRecord> getBamIterator(Reads readClass, String chrom, int startpos, int endpos) {
		
		try {
			if(samFileReader != null) {
				samFileReader.close();
			}
			if(CRAMReader != null) {				
				CRAMReader.close();
			}
			
			if(readClass.sample.samFile.getName().endsWith(".list")) {
				File file = getListFile(readClass.sample.samFile,chrom);
				if(file != null) {
					if(file.getName().endsWith(".cram")) {
						 readClass.sample.samFile = file;
						 readClass.sample.CRAM = true;
						 readClass.sample.readString = "CRAM";
						 CRAMReader = new CRAMFileReader(readClass.sample.samFile, new File(readClass.sample.samFile.getCanonicalFile() +".crai"),new ReferenceSource(Main.ref),Main.referenceFile, ValidationStringency.SILENT);
							
					}
					else {
					//samFileReader = SamReaderFactory.make().open(file);
						samFileReader = new SAMFileReader(file);
					}
				}
			}
			else {
				if(readClass.sample.CRAM) {					
					
					 CRAMReader = new CRAMFileReader(readClass.sample.samFile, new File(readClass.sample.samFile.getCanonicalFile() +".crai"),new ReferenceSource(Main.ref),Main.referenceFile, ValidationStringency.SILENT);
					 if(CRAMReader != null && !CRAMReader.hasIndex()) {			
							Main.showError("Index file is missing (.crai) for " +readClass.sample.samFile.getName(), "Note");
							return null;
						}
				}
				else {
					try {
					samFileReader = new SAMFileReader(readClass.sample.samFile);
					if(samFileReader != null && !samFileReader.hasIndex()) {	
						Main.showError("Index file is missing (.bai) for " +readClass.sample.samFile.getName(), "Note");
						return null;
					}
					}
					catch(Exception e) {
						Main.showError(e.getMessage(), "Note");
						e.printStackTrace();
					}
				}				
			}			
		}
		catch(Exception e) {
			e.printStackTrace();
		}		
		
		try {			
			if(readClass.sample.CRAM) {
				
				QueryInterval[] interval = { new QueryInterval(CRAMReader.getFileHeader().getSequence(chrom).getSequenceIndex(), startpos, endpos) };				
				Iterator<SAMRecord> value = CRAMReader.query(interval, false);
				if(!readClass.sample.chrSet) {
					if(!value.hasNext()) {
						value = CRAMReader.query(interval, false);
						if(value.hasNext()) {
							readClass.sample.chr = "chr";
							readClass.sample.chrSet = true;
						}
					}
					else {
						readClass.sample.chrSet = true;
					}
				}
				return value;
			}
			else {
				Iterator<SAMRecord> value = null;
			//	SAMReadGroupRecord record = new SAMReadGroupRecord("fixed_My6090N_sorted");
			//	samFileReader.getFileHeader().addReadGroup(record);
				try {
					SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
					value = samFileReader.queryOverlapping(chrom, startpos, endpos);	
				}
				catch(htsjdk.samtools.SAMFormatException e) {
					e.printStackTrace();
				}
				if(!readClass.sample.chrSet) {
					if(!value.hasNext()) {
						value = samFileReader.queryOverlapping("chr" +chrom, startpos, endpos);
						if(value.hasNext()) {
							readClass.sample.chr = "chr";
							readClass.sample.chrSet = true;
						}
					}
					else {
						readClass.sample.chrSet = true;
					}
				}
				
				return value;
			}
		}
		catch(Exception e) {		
		
			try {
				if(readClass.sample.CRAM) {
					QueryInterval[] interval = { new QueryInterval(CRAMReader.getFileHeader().getSequence(readClass.sample.chr + "M").getSequenceIndex(), startpos, endpos) };
					
					Iterator<SAMRecord> value = CRAMReader.query(interval, false);
					
					if(!readClass.sample.chrSet) {
						if(!value.hasNext()) {
							QueryInterval[] interval2 = { new QueryInterval(CRAMReader.getFileHeader().getSequence("chrM").getSequenceIndex(), startpos, endpos) };
							
							value = CRAMReader.query(interval2, false);
							if(value.hasNext()) {
								readClass.sample.chr = "chr";
								readClass.sample.chrSet = true;
							}
						}
						else {
							readClass.sample.chrSet = true;
						}
					}
					return value;
				}
				else {
					Iterator<SAMRecord> value = samFileReader.queryOverlapping(readClass.sample.chr+"M", startpos, endpos);
					if(!readClass.sample.chrSet) {
						if(!value.hasNext()) {
							value = samFileReader.queryOverlapping("chrM", startpos, endpos);
							if(value.hasNext()) {
								readClass.sample.chr = "chr";
								readClass.sample.chrSet = true;
							}
						}
						else {
							readClass.sample.chrSet = true;
						}
					}
					return value;
				}
				
			}
			catch(Exception ex) {
				
				ex.printStackTrace();
				ErrorLog.addError(e.getStackTrace());
				return null;
			}			
		}	
	}
	public SAMRecord getRead(String chrom, int startpos, int endpos, String name, Reads readClass) {		
		
		Iterator<SAMRecord> bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos-1,  endpos);
		SAMRecord samRecord = null;
		
		while(bamIterator != null && bamIterator.hasNext()) {	
			
				try {
					samRecord = bamIterator.next(); 
					if(samRecord.getUnclippedStart() < startpos) {
						continue;
					}
					if(samRecord.getUnclippedStart() == startpos && samRecord.getReadName().equals(name)) {
						
						return samRecord;
					}
				}
				catch(Exception e) {
					e.printStackTrace();
				}
			
		}
		return null;
	}
	
	public ArrayList<ReadNode> getReads(String chrom, int startpos, int endpos, Reads readClass) {
	//	ReadAdder readrunner = null;
		try {
		
		//Sample sample = readClass.sample;
		
		if(Main.drawCanvas.drawVariables.sampleHeight > 100) {
			if(viewLength > Settings.readDrawDistance) {				
				
				if(firstCov || readClass.getCoverageStart() == 0 || (readClass.getCoverageStart() < startpos && readClass.getCoverageEnd()  > endpos) || (readClass.getCoverageStart() > startpos && readClass.getCoverageEnd() < endpos) ) {
					
					double[][] coverages = new double[(int)Main.frame.getWidth()][8];
					readClass.setCoverages(coverages);
					readClass.setMaxcoverage(1);
					readClass.setCoverageStart(startpos);
					oldstart = endpos;
					readClass.setCoverageEnd(endpos);
					bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos,  endpos);
					firstCov = false;
				}
				else if(readClass.getCoverageStart()  > startpos) {
				/*	
					int[][] coverages = new int[(int)(readClass.getCoverages().length +((readClass.getCoverageStart()-startpos)*pixel))][8];
					
					pointer = coverages.length-1;
					for(int i = readClass.getCoverages().length-1; i>=0; i--) {				
						coverages[pointer][0] = readClass.getCoverages()[i][0];						
						pointer--;
					}
					left = true;
					oldstart = readClass.getCoverageStart();				
					bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos, readClass.getCoverageStart());
					readClass.setCoverageStart(startpos);
					readClass.setCoverages(coverages);
					*/
					double[][] coverages = new double[(int)Main.frame.getWidth()][8];
					
					readClass.setCoverages(coverages);
					readClass.setMaxcoverage(1);
					bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos, endpos);
			//		readClass.setReadEnd(endpos);
				
					readClass.setCoverageStart(startpos);
					readClass.setCoverageEnd(endpos);
					oldstart = endpos;
				}
				else if(readClass.getCoverageEnd()  < endpos) {
					double[][] coverages = new double[(int)Main.frame.getWidth()][8];
					
					readClass.setCoverages(coverages);
					readClass.setMaxcoverage(1);
					bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos, endpos);
			//		readClass.setReadEnd(endpos);
				
					readClass.setCoverageStart(startpos);
					readClass.setCoverageEnd(endpos);
					oldstart = endpos;
				}
				else {
				
					return null;
				}
				
			}
			else {
				firstCov = true;
			
			if(readClass.getReads().isEmpty()) {
					
					right = false;
					left = false;		
					
					if(viewLength <= Settings.readDrawDistance && !chrom.contains("M")) {
						searchwindow = 10000;						
					}
					else {						
						searchwindow = 0;
					}
					
					if(firstSample) {	
						
						splitIndex.setReference(new ReferenceSeq(splitIndex.chrom,startpos-searchwindow-1,  endpos+searchwindow +200, Main.referenceFile));
						
					}
					try {
						
						bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos-searchwindow,  endpos+searchwindow  );
						
					}
					catch(Exception e) {
						e.printStackTrace();
					/*	try {
						
						}
						catch(java.lang.IllegalArgumentException ex) {
							readClass.loading = false;
							System.out.println("Chromosome " +chrom +" not found in " +readClass.sample.samFile.getName());
							ErrorLog.addError(ex.getStackTrace());
							return null;
						}	*/					
					}
					
					startpos = startpos - searchwindow;
					endpos = endpos + searchwindow;
					
					if(startpos < 1) {
						startpos = 1;
					}
					
					readClass.setCoverageStart(startpos);
					readClass.startpos = startpos;
					readClass.endpos = endpos;
					readClass.setReadStart(startpos);
					readClass.setReadEnd(endpos);
					
					if(splitIndex.getMinReadStart() > startpos) {
						splitIndex.setMinReadStart(startpos);
					}
					if(splitIndex.getMaxReadEnd() < endpos) {
						splitIndex.setMaxReadEnd(endpos);
					}
					double[][] coverages = new double[endpos-startpos+searchwindow*2][8];
					readClass.setCoverageEnd(startpos + coverages.length);
					readClass.setCoverages(coverages);
					readClass.setMaxcoverage(1);
					
			}
			else if(!readClass.getReads().isEmpty() && readClass.getReadStart() /*(readClass.getFirstRead().getPosition()*/ < startpos && readClass.getReadEnd() /*readClass.getLastRead().getPosition()*/ > endpos){
				left = false;
				right = false;
				
				if(viewLength <= Settings.readDrawDistance && !chrom.contains("M")) {
					searchwindow = viewLength;
				}
				else {
					searchwindow = 0;
				}
			
				return readClass.getReads();
			}
			else if(readClass.getLastRead() != null && readClass.getReadEnd() < endpos) {
				
				if(viewLength <= Settings.readDrawDistance && !chrom.contains("M")) {
					searchwindow = viewLength;
				}
				else {
					searchwindow = 0;			
				}				
				
				right = true;
				left = false;								
				startpos = readClass.getReadEnd();				
				endpos = endpos+searchwindow;
				addlength = 1000;
			//	readClass.reference.append( endpos+200);
				if(firstSample) {	
					if(splitIndex.getReference() == null) {
						splitIndex.setReference(new ReferenceSeq(splitIndex.chrom,startpos-searchwindow-1,  endpos+searchwindow +200, Main.referenceFile));
						
					}
					else {
						splitIndex.getReference().append( endpos+200);
					}
				}
				bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos,  endpos );
				
				if(readClass.getCoverageEnd() < endpos) {
					if(searchwindow < 1000) {
						addlength = 2000;
					}
					else if(searchwindow < 10000) {
						addlength = 20000;
					}
					else {
						addlength = 200000;
					}
				}
				double[][] coverages = new double[(int)(readClass.getCoverages().length +(addlength))][8];				
				
				for(int i =0; i<readClass.getCoverages().length; i++) {
			
					coverages[i][0] = readClass.getCoverages()[i][0];
					coverages[i][1] = readClass.getCoverages()[i][1];
					coverages[i][2] = readClass.getCoverages()[i][2];
					coverages[i][3] = readClass.getCoverages()[i][3];
					coverages[i][4] = readClass.getCoverages()[i][4];
					coverages[i][5] = readClass.getCoverages()[i][5];
					coverages[i][6] = readClass.getCoverages()[i][6];
					coverages[i][7] = readClass.getCoverages()[i][7];
			
					
				}
				readClass.endpos = endpos;
				if(splitIndex.getMaxReadEnd() < endpos) {
					splitIndex.setMaxReadEnd(endpos);
				}
				readClass.setReadEnd(endpos);				
				readClass.setCoverageEnd(readClass.getCoverageStart()+coverages.length);
				readClass.setCoverages(coverages);
			
			
			}
			else if(readClass.getFirstRead() != null && readClass.getFirstRead().getPosition() > startpos) {
				left = true;
				right = false;
				
				
				if(viewLength <= Settings.readDrawDistance && !chrom.contains("M")) {
					searchwindow = viewLength;
				}
				else {
					searchwindow = 0;
				}					
				
				endpos = readClass.getReadStart();				
				startpos = startpos-searchwindow;
				
				if(firstSample) {			
					splitIndex.getReference().appendToStart( startpos-1);
				}
				bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos,  endpos+100 );
				addlength = 0;
				if(searchwindow < 1000) {
					addlength = 2000;
				}
				else if(searchwindow < 10000) {
					addlength = 20000;
				}
				else {
					addlength = 200000;
				}
				double[][] coverages = new double[(int)(readClass.getCoverages().length +addlength)][8];
				readClass.setCoverageStart(readClass.getCoverageStart()-addlength);
				pointer = coverages.length-1;
				
				for(int i = readClass.getCoverages().length-1; i>=0; i--) {			
					coverages[pointer][0] = readClass.getCoverages()[i][0];
					coverages[pointer][1] = readClass.getCoverages()[i][1];
					coverages[pointer][2] = readClass.getCoverages()[i][2];
					coverages[pointer][3] = readClass.getCoverages()[i][3];
					coverages[pointer][4] = readClass.getCoverages()[i][4];
					coverages[pointer][5] = readClass.getCoverages()[i][5];
					coverages[pointer][6] = readClass.getCoverages()[i][6];
					coverages[pointer][7] = readClass.getCoverages()[i][7];
					pointer--;					
				}
				left = true;				
				readClass.startpos = startpos;			
				readClass.setReadStart(startpos);
				if(splitIndex.getMinReadStart() > startpos) {
					splitIndex.setMinReadStart(startpos);
				}
				
				readClass.setCoverages(coverages);				
			}
			}
			
			coverageArray = readClass.getCoverages();			
			if(readClass.getReadStart() == null) {
				readClass.setReadStart(startpos);
			}
			
			while(bamIterator != null && bamIterator.hasNext()) {	
				
				try {	
					
				if(cancelreadread || !Main.drawCanvas.loading ||readClass.sample.getIndex() < Main.drawCanvas.drawVariables.visiblestart || readClass.sample.getIndex() > Main.drawCanvas.drawVariables.visiblestart + Main.drawCanvas.drawVariables.visiblesamples) {
					
		  			return null;
		  		}
				try {
					samRecord = bamIterator.next(); 	
					
				}
				catch(htsjdk.samtools.SAMFormatException ex) {
					ex.printStackTrace();					
				}
				
				if(samRecord.getReadUnmappedFlag()) {						
					
					continue;
				}
				if(Draw.variantcalculator) {
					if(readClass.getHeadAndTail().size() > Settings.readDepthLimit) {
						
						break;
					}
				}
				//TODO jos on pienempi ku vika unclipped start		
				
				
				if(samRecord.getUnclippedEnd() < startpos) { //this.readSeqStart+1) {
					
					continue;
				}
				
				if(samRecord.getUnclippedStart() >= endpos){
					Main.drawCanvas.loadBarSample = 0;  
					Main.drawCanvas.loadbarAll  =0;
					
					break;
				}		
				
				if(readClass.sample.longestRead <samRecord.getCigar().getReferenceLength()) {
					readClass.sample.longestRead = samRecord.getCigar().getReferenceLength();
					
				}
				if(readClass.sample.getComplete() == null) {
					if(samRecord.getReadName().startsWith("GS")) {
						readClass.sample.setcomplete(true);
					}
					else {
						readClass.sample.setcomplete(false);
					}
				}
				
			if(viewLength < Settings.readDrawDistance) {
				/*
					if(splitIndex.reference == null) {
						splitIndex.reference = new ReferenceSeq(chrom, samRecord.getUnclippedStart()-2000, end+1000, Main.referenceFile);
						
					}
					else if(samRecord.getUnclippedStart()-1 < splitIndex.reference.getStartPos()) {
						splitIndex.reference.appendToStart(samRecord.getUnclippedStart()-1);
					}
					else if(samRecord.getUnclippedEnd() > splitIndex.reference.getEndPos()) {
						splitIndex.reference.append(samRecord.getUnclippedEnd());
					}
				*/
				if(samRecord.getUnclippedEnd() > readClass.getCoverageEnd()-1) {
					
					double[][] coverages = new double[(int)(readClass.getCoverages().length +(samRecord.getUnclippedEnd()-samRecord.getUnclippedStart() + (int)splitIndex.viewLength))][8];
					
					//pointer = coverages.length-1;
					for(int i =0; i<readClass.getCoverages().length; i++) {
						
						coverages[i][0] = readClass.getCoverages()[i][0];
						coverages[i][1] = readClass.getCoverages()[i][1];
						coverages[i][2] = readClass.getCoverages()[i][2];
						coverages[i][3] = readClass.getCoverages()[i][3];
						coverages[i][4] = readClass.getCoverages()[i][4];
						coverages[i][5] = readClass.getCoverages()[i][5];
						coverages[i][6] = readClass.getCoverages()[i][6];
						coverages[i][7] = readClass.getCoverages()[i][7];
					//	pointer--;
						
					}						
					readClass.setCoverageEnd(readClass.getCoverageStart()+coverages.length);
					readClass.setCoverages(coverages);
					coverageArray = readClass.getCoverages();
					
				}
				else if(samRecord.getUnclippedStart() > 0 && samRecord.getUnclippedStart() < readClass.getCoverageStart() ) {
					
					/*
					double[][] coverages = null;
					if(readClass.getCoverageStart()-(readClass.getCoverageStart()-samRecord.getUnclippedStart())-1000 < 1) {
						coverages = new double[(int)(readClass.getCoverages().length + readClass.getCoverageStart())][8];
						readClass.setCoverageStart(1);
					}
					else {
						coverages = new double[(int)(readClass.getCoverages().length + (readClass.getCoverageStart()-samRecord.getUnclippedStart())+1000)][8];
						readClass.setCoverageStart(readClass.getCoverageStart()-(readClass.getCoverageStart()-samRecord.getUnclippedStart())-1000);
						
					}
					pointer = coverages.length-1;
					for(int i = readClass.getCoverages().length-1; i>=0; i--) {			
						coverages[pointer][0] = readClass.getCoverages()[i][0];
						coverages[pointer][1] = readClass.getCoverages()[i][1];
						coverages[pointer][2] = readClass.getCoverages()[i][2];
						coverages[pointer][3] = readClass.getCoverages()[i][3];
						coverages[pointer][4] = readClass.getCoverages()[i][4];
						coverages[pointer][5] = readClass.getCoverages()[i][5];
						coverages[pointer][6] = readClass.getCoverages()[i][6];
						coverages[pointer][7] = readClass.getCoverages()[i][7];
						pointer--;					
					}
					
					readClass.setCoverages(coverages);			
					*/
				}
			}
			
			
			 if(((viewLength < Settings.readDrawDistance && samRecord.getUnclippedStart() >= startpos) || (viewLength >= Settings.readDrawDistance && samRecord.getUnclippedEnd() >= startpos)) && samRecord.getUnclippedStart() <= endpos) {					
					java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> mismatches = null;
					
					if(samRecord.getMappingQuality() >= Settings.mappingQ) {
						
						if(viewLength > Settings.readDrawDistance) {
							
							if( Main.drawCanvas.loadBarSample != (int)(((samRecord.getUnclippedStart()-startpos)/(double)(oldstart-startpos))*100) ) {
								
								 Main.drawCanvas.loadBarSample = (int)(((samRecord.getUnclippedStart()-startpos)/(double)(oldstart-startpos))*100);   
								 Main.drawCanvas.loadbarAll = (int)(((samRecord.getUnclippedStart()-startpos)/(double)(oldstart-startpos))*100);
							}
					
							for(int i = 0; i<(int)(samRecord.getReadLength()*splitIndex.pixel +1); i++) {
								if((samRecord.getUnclippedStart()-start)*splitIndex.pixel +i < 0) {
									continue;
								}
								if((samRecord.getUnclippedStart()-start)*splitIndex.pixel + i > coverageArray.length-1) {
									break;
								}
							
								if(samRecord.getReadLength()*splitIndex.pixel >= 1) {
								
									coverageArray[(int)((samRecord.getUnclippedStart()-start)*splitIndex.pixel+i)][0]+=(1/splitIndex.pixel); 
								}
								else {
									
									coverageArray[(int)((samRecord.getUnclippedStart()-start)*splitIndex.pixel+i)][0]+=samRecord.getReadLength(); 
									
									
								}
								if(coverageArray[(int)((samRecord.getUnclippedStart()-start)*splitIndex.pixel+i)][0]/*pixel*/ > readClass.getMaxcoverage()) {
									readClass.setMaxcoverage((coverageArray[(int)((samRecord.getUnclippedStart()-start)*splitIndex.pixel+i)][0]/*pixel*/));
								
								}
								
							}
							
						}
						else {
							if( Main.drawCanvas.loadBarSample != (int)(((samRecord.getUnclippedStart()-startpos)/(double)(endpos-startpos))*100) ) {
								
								 Main.drawCanvas.loadBarSample = (int)(((samRecord.getUnclippedStart()-startpos)/(double)(endpos-startpos))*100);   
								 Main.drawCanvas.loadbarAll = (int)(((samRecord.getUnclippedStart()-startpos)/(double)(endpos-startpos))*100);
							}
							
							
							
						}
					}
					
					if(viewLength <= Settings.readDrawDistance) {	
					
						if(right && readClass.getLastRead() != null && (samRecord.getUnclippedStart() <= startpos)) {
							
							continue;
						}
						
						if(right && samRecord.getUnclippedStart() > endpos) {
							right = false;
							
						}
						mismatches = getMismatches(samRecord, readClass, coverageArray);
						if(left) {
							
							if(readClass.getFirstRead().getPosition() > samRecord.getUnclippedStart()) {							
								
								ReadNode addNode = new ReadNode(samRecord, readClass.sample.getComplete(), chrom, readClass.sample, splitIndex,readClass, mismatches);
								readClass.setFirstRead(addNode);
								startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-(readClass.readHeight+2)));					
						//		addNode.setRectBounds((int)((addNode.getPosition()-start)*pixel), startY, (int)(samRecord.getReadLength()*pixel), Main.drawCanvas.readHeight);
														
								addNode.setPrev(null);
								addNode.setNext(readClass.getHeadAndTail().get(0)[headnode]);
								readClass.getHeadAndTail().get(0)[headnode].setPrev(addNode);
								lastAdded = addNode;
						
							}
							else {
								
								ReadNode addNode = new ReadNode(samRecord, readClass.sample.getComplete(), chrom, readClass.sample, splitIndex,readClass, mismatches);
								startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-(readClass.readHeight+2)));					
								addNode.setPrev(lastAdded);
								addNode.setNext(lastAdded.getNext());
								lastAdded.setNext(addNode);
								lastAdded = addNode;
								readClass.getHeadAndTail().get(0)[headnode].setPrev(addNode);
								addNode.setNext(readClass.getHeadAndTail().get(0)[headnode]);
								
							}
						}
						else {	
							
							try {
								
								readSam(chrom, readClass, samRecord, mismatches);
								
							}
							catch(Exception e) {
								System.out.println(samRecord +" " +currentread);
								e.printStackTrace();
								
								break;
							}
						}
					}
				}
				
				else {
				
					continue;
				}
				}
				catch(Exception e) {
					if(readClass.getReadStart() > startpos) {
						readClass.setReadStart(startpos);
					}
					
					if(readClass.getReadEnd() < endpos) {
						readClass.setReadEnd(endpos);
					}
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
					break;
				}				
			}	
			if(left) {
				
				Main.drawCanvas.loadBarSample = 0;  
				Main.drawCanvas.loadbarAll  =0;
				setLevels(lastAdded, readClass.sample, readClass);					
				lastAdded = null;
				
			}
		}
		
		}
		catch(Exception ex) {
			if(readClass == null) {
				
				sample.resetreadHash();
				
			}
			if(readClass.getReadStart() == null) {
				readClass.setReadStart(startpos);
			}
			else if(readClass.getReadStart() > startpos) {
				readClass.setReadStart(startpos);
			}
			if(readClass.getReadEnd() < endpos) {
				readClass.setReadEnd(endpos);
			}
			
			
			ErrorLog.addError(ex.getStackTrace());
			ex.printStackTrace();		
			
			return null;
	//		readrunner.end = true;
		}
		
		return null;
	}
	double[][] coverageArrayAdd(int addlength, double[][] coverageArray) {
	
		double[][] coverages = new double[(int)(coverageArray.length +(addlength)+1000)][8];
	
	//pointer = coverages.length-1;
	for(int i =0; i<coverageArray.length; i++) {

		coverages[i][0] = coverageArray[i][0];
		coverages[i][1] = coverageArray[i][1];
		coverages[i][2] = coverageArray[i][2];
		coverages[i][3] = coverageArray[i][3];
		coverages[i][4] = coverageArray[i][4];
		coverages[i][5] = coverageArray[i][5];
		coverages[i][6] = coverageArray[i][6];
		coverages[i][7] = coverageArray[i][7];
		
	}
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
	readClass.setCoverageEnd(readClass.getCoverageStart()+coverages.length);
	readClass.setCoverages(coverages);
	return coverages;
	
}
	
	
java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> getMismatches(SAMRecord samRecord, Reads readClass, double[][] coverageArray) {	
	
	if(readClass == null) {
		sample.resetreadHash();
	}
	java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> mismatches = null;
	String MD = samRecord.getStringAttribute("MD");	
	String SA = samRecord.getStringAttribute("SA");
	
	try {
	
	if(MD == null) {
		
		if(readClass.sample.MD) {
			readClass.sample.MD = false;
		}
	}
	
	if((readClass.sample.MD || MD!=null) && (samRecord.getCigarString().contains("S") && Settings.softClips == 0)/* && ((!samRecord.getCigarString().contains("S") && SA == null) || SA !=null)*/) {		
		/*if(samRecord.getReadName().contains("1104:8947:11927")) {
			System.out.println(readClass.sample.MD +" " +MD +" " +samRecord.getCigarString().contains("S") +" " +Settings.softClips);
		}*/
		if(!readClass.sample.MD) {
			readClass.sample.MD = true;
		}
	
		java.util.ArrayList<java.util.Map.Entry<Integer,Integer>> insertions = null;
		int softstart = 0;
		if(samRecord.getCigarLength() > 1) {
			
			readstart = samRecord.getUnclippedStart();	
			if(readClass.getCoverageStart()>readstart) {
				return null;
			}
			readpos= 0;
			for(int i = 0; i<samRecord.getCigarLength(); i++) {
				element = samRecord.getCigar().getCigarElement(i);
				if(element.getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH)== 0) {
					for(int r = readpos; r<readpos+element.getLength(); r++) {					
						
						coverageArray[((readstart+r)-readClass.getCoverageStart())][0]++; 
						
						if(coverageArray[((readstart+r)-readClass.getCoverageStart())][0] > readClass.getMaxcoverage()) {
							readClass.setMaxcoverage(coverageArray[((readstart+r)-readClass.getCoverageStart())][0]);
							
						}
						mispos++;
					}
					
					readpos+=element.getLength();
				}
				else if(element.getOperator().compareTo(CigarOperator.DELETION)== 0) {
					for(int d = 0; d<element.getLength(); d++) {
						readClass.getCoverages()[(readstart+readpos+d)- readClass.getCoverageStart()][Main.baseMap.get((byte)'D')]++;
						
					}
					readpos+=element.getLength();
				}
				else if(element.getOperator().compareTo(CigarOperator.INSERTION)== 0) {	
					if(insertions == null) {
						insertions = new java.util.ArrayList<java.util.Map.Entry<Integer,Integer>>();						
					}
					java.util.Map.Entry<Integer, Integer> ins = new java.util.AbstractMap.SimpleEntry<>(readpos, element.getLength());
					insertions.add(ins);
					readClass.getCoverages()[(readstart+readpos)- readClass.getCoverageStart()][Main.baseMap.get((byte)'I')]++;
		
					mispos+=element.getLength();
					
				
				}
				else if(element.getOperator().compareTo(CigarOperator.SOFT_CLIP)== 0 ) {
					
						if(i==0 && SA == null) {						
							softstart = element.getLength();
						}				
					
						if(Settings.softClips == 1) {
							for(int r = readpos; r<readpos+element.getLength(); r++) {
								if(SA == null) {
									coverageArray[((readstart+r)-readClass.getCoverageStart())][0]++; 
								}
								mispos++;
							}
							readpos+=element.getLength();	
							
						}
						else {
							if(i == 0) {
							
								readstart = samRecord.getAlignmentStart();
							
								mispos+=element.getLength();
							}
						}													
					
				}
				else if (element.getOperator().compareTo(CigarOperator.HARD_CLIP)== 0) {
					if(i == 0) {
						readstart = samRecord.getAlignmentStart();
					}
				}
				else if(element.getOperator().compareTo(CigarOperator.SKIPPED_REGION)== 0) {										
					readpos+=element.getLength();
				}
			}
		}
		else {
			readstart = samRecord.getUnclippedStart();
			
			for(int i = 0; i<samRecord.getReadLength(); i++) {
				try {			
					if(readClass.getCoverageStart()>readstart) {
						break;
					}
						coverageArray[(readstart+i)-readClass.getCoverageStart()][0]++; 
					
					if(coverageArray[(samRecord.getUnclippedStart()-readClass.getCoverageStart())+i][0] > readClass.getMaxcoverage()) {
						readClass.setMaxcoverage(coverageArray[(samRecord.getUnclippedStart()-readClass.getCoverageStart())+i][0]);
					}
				}
				catch(Exception e) {
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
					FileRead.cancelreadread = true;
					break;
				}
			}		
		}
		if(MD != null) {
			try {
				Integer.parseInt(MD);
				return null;
			}
			catch(Exception ex) {
				
				readstart = samRecord.getAlignmentStart()-softstart;		
				
				char[] chars = MD.toCharArray();
				String bases = samRecord.getReadString();
				String qualities = samRecord.getBaseQualityString();
				
			/*	if(SA == null) {
					softstart = 0;
				}*/
				int readpos = softstart;
				int mispos = softstart;
						
				int add =0;				
				int digitpointer = 0, insertionpointer = 0;
				boolean first = true;
				
				for(int i = 0; i<chars.length; i++) {
					try {
					if(insertions != null) {
						if(insertionpointer < insertions.size() && insertions.get(insertionpointer).getKey() <= readpos) {
							while(insertionpointer < insertions.size() && insertions.get(insertionpointer).getKey() <= readpos) {
								mispos+=insertions.get(insertionpointer).getValue();
								insertionpointer++;
							}
						}						
					}
					if(Character.isDigit(chars[i])) {
						digitpointer = i+1;
						while(digitpointer < chars.length && Character.isDigit(chars[digitpointer])) {
							digitpointer++;
						}
					
						if(digitpointer== chars.length) {
							//System.out.println(MD.substring(i));
							if((mispos +Integer.parseInt(MD.substring(i,digitpointer))) > samRecord.getReadLength()) {
								if(!first) {
									break;
								}
									first = false;								
									readpos = 0;
									mispos = 0;
									i = -1;
									digitpointer = 0;
									insertionpointer = 0;
									mismatches.clear();
									continue;
							}
						
							break;
						}
						else {
							add = Integer.parseInt(MD.substring(i,digitpointer));
							readpos += add;		
							mispos += add;
						}						
						i = digitpointer-1;
					}
					else if(chars[i] != '^') {
						
						if(mismatches == null) {
							mismatches = new java.util.ArrayList<java.util.Map.Entry<Integer,Byte>>();
						}
						
						readClass.getCoverages()[(readstart+readpos)- readClass.getCoverageStart()][Main.baseMap.get((byte)bases.charAt(mispos))]++;
						
							if(samRecord.getBaseQualityString().length() != 1 && (int)qualities.charAt(mispos)-readClass.sample.phred < Settings.baseQ) {
								
								java.util.Map.Entry<Integer, Byte> mismatch = new java.util.AbstractMap.SimpleEntry<>(readpos, (byte)Character.toLowerCase(bases.charAt(mispos)));
								mismatches.add(mismatch);
							}
							else {
								
								java.util.Map.Entry<Integer, Byte> mismatch = new java.util.AbstractMap.SimpleEntry<>(readpos, (byte)bases.charAt(mispos));
								mismatches.add(mismatch);
							}
						mispos++;
						readpos++;
					}
					else {
					
						digitpointer = i+1;
					
						while(!Character.isDigit(chars[digitpointer])) {
							digitpointer++;
							readpos++;
						}
						
						i=digitpointer-1;
					}
					}
					catch(Exception exc) {
						if(!first) {
							break;
						}
						first = false;			
						readpos = 0;
						mispos = 0;
						i = -1;
						digitpointer = 0;
						insertionpointer = 0;
						mismatches.clear();
					}
				}				
			}			
		}
	}
	else {			
		timecounter = 0;
		while(splitIndex.getReference() == null) {
			Thread.sleep(100);
			timecounter++;
			if(timecounter > 20) {
				break;
			}
		}
		Thread.sleep(0);
		if(samRecord.getCigarLength() > 1) {			
			readstart = samRecord.getUnclippedStart();			
			readpos= 0;
			mispos= 0;
			if(readClass.getCoverageStart()>readstart) {
				return null;
			}
			if(mismatches == null) {
				mismatches = new java.util.ArrayList<java.util.Map.Entry<Integer,Byte>>();
			}
			for(int i = 0; i<samRecord.getCigarLength(); i++) {
				element = samRecord.getCigar().getCigarElement(i);
				if(element.getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH)== 0) {
					for(int r = readpos; r<readpos+element.getLength(); r++) {																	
						if(((readstart+r)-splitIndex.getReference().getStartPos()-1) > splitIndex.getReference().getSeq().length-1) {
							splitIndex.getReference().append(samRecord.getUnclippedEnd()+1000);
						}
						if(((readstart+r)-splitIndex.getReference().getStartPos()-1) < 0) {
							splitIndex.getReference().appendToStart(samRecord.getUnclippedStart()-1000);
						}
						
						if(samRecord.getReadBases()[mispos] !=  splitIndex.getReference().getSeq()[((readstart+r)-splitIndex.getReference().getStartPos()-1)]) {	
							readClass.getCoverages()[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get((byte)samRecord.getReadBases()[mispos])]++;
							
							if(samRecord.getBaseQualityString().length() != 1 && (int)samRecord.getBaseQualityString().charAt(mispos)-readClass.sample.phred < Settings.baseQ ) {
								
								java.util.Map.Entry<Integer, Byte> mismatch = new java.util.AbstractMap.SimpleEntry<>(r, (byte)Character.toLowerCase((char)samRecord.getReadBases()[mispos]));
								mismatches.add(mismatch);
							}
							else {
								
								java.util.Map.Entry<Integer, Byte> mismatch = new java.util.AbstractMap.SimpleEntry<>(r, samRecord.getReadBases()[mispos]);
								mismatches.add(mismatch);
							}
						}
					
						coverageArray[((readstart+r)-readClass.getCoverageStart())][0]++; 
						
						if(coverageArray[((readstart+r)-readClass.getCoverageStart())][0] > readClass.getMaxcoverage()) {
							readClass.setMaxcoverage(coverageArray[((readstart+r)-readClass.getCoverageStart())][0]);
						}
						mispos++;
					}
					
					readpos+=element.getLength();
				}
				else if(element.getOperator().compareTo(CigarOperator.DELETION)== 0) {
			//		readWidth+=element.getLength();						
					for(int d = 0; d<element.getLength(); d++) {
						readClass.getCoverages()[(readstart+readpos+d)- readClass.getCoverageStart()][Main.baseMap.get((byte)'D')]++;
						
					}
					readpos+=element.getLength();
				}
				else if(element.getOperator().compareTo(CigarOperator.INSERTION)== 0) {										
					readClass.getCoverages()[(readstart+readpos)- readClass.getCoverageStart()][Main.baseMap.get((byte)'I')]++;
					mispos+=element.getLength();					
				
				}
				else if(element.getOperator().compareTo(CigarOperator.SOFT_CLIP)== 0 ) {
						if(Settings.softClips == 1) {
							for(int r = readpos; r<readpos+element.getLength(); r++) {
								if(((readstart+r)-splitIndex.getReference().getStartPos()-1) > splitIndex.getReference().getSeq().length-1) {
									splitIndex.getReference().append(samRecord.getUnclippedEnd()+1000);
								}
								if(((readstart+r)-splitIndex.getReference().getStartPos()-1) < 0) {
									splitIndex.getReference().appendToStart(samRecord.getUnclippedStart()-1000);
								}
								/*	if(((readstart+r)-readClass.reference.getStartPos()-1) < 0) {
									continue;
								}*/
								if(samRecord.getReadBases()[mispos] != splitIndex.getReference().getSeq()[((readstart+r)-splitIndex.getReference().getStartPos()-1)]) {													
									if(SA == null) {
									
										readClass.getCoverages()[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get(samRecord.getReadBases()[mispos])]++;
									}
									if(mismatches == null) {
										mismatches = new java.util.ArrayList<java.util.Map.Entry<Integer,Byte>>();
									}
									java.util.Map.Entry<Integer, Byte> mismatch = new java.util.AbstractMap.SimpleEntry<>(r, samRecord.getReadBases()[mispos]);													
									mismatches.add(mismatch);												
								}
								if(SA == null) {
									coverageArray[((readstart+r)-readClass.getCoverageStart())][0]++; 
								}
								mispos++;
							}
							readpos+=element.getLength();	
							
						}
						else {
							if(i == 0) {							
								readstart = samRecord.getAlignmentStart();							
								mispos+=element.getLength();
							}
						}				
														
					
				}
				else if (element.getOperator().compareTo(CigarOperator.HARD_CLIP)== 0) {
					if(i == 0) {
						readstart = samRecord.getAlignmentStart();
					}
				}
				else if(element.getOperator().compareTo(CigarOperator.SKIPPED_REGION)== 0) {										
					readpos+=element.getLength();										
		//			this.readWidth+=element.getLength();
				}
			}
		}
		else {
			readstart = samRecord.getUnclippedStart();
			if(readClass.getCoverageStart()>readstart) {
				return null;
			}
			for(int r = 0; r<samRecord.getReadLength(); r++) {
				try {									
					if(((readstart+r)-splitIndex.getReference().getStartPos()-1) > splitIndex.getReference().getSeq().length-1) {
						splitIndex.getReference().append(samRecord.getUnclippedEnd()+1000);
					}
					if(((readstart+r)-splitIndex.getReference().getStartPos()-1) < 0) {
						splitIndex.getReference().appendToStart(samRecord.getUnclippedStart()-1000);
					}
				/*	if(((readstart+r)-readClass.reference.getStartPos()-1) > readClass.reference.getSeq().length-1) {
						readClass.reference.append(samRecord.getUnclippedEnd()+1000);
					}
					if(((readstart+r)-readClass.reference.getStartPos()-1) < 0) {
						continue;
					}*/
					if(samRecord.getReadBases()[r] != splitIndex.getReference().getSeq()[((readstart+r)-splitIndex.getReference().getStartPos()-1)]) {									
						
						if((readstart+r)-readClass.getCoverageStart() < readClass.getCoverages().length-1) {
							readClass.getCoverages()[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get(samRecord.getReadBases()[r])]++;
							if(mismatches == null) {
								mismatches = new java.util.ArrayList<java.util.Map.Entry<Integer,Byte>>();
							}
							java.util.Map.Entry<Integer, Byte> mismatch = new java.util.AbstractMap.SimpleEntry<>(r, samRecord.getReadBases()[r]);												
							mismatches.add(mismatch);
						}									
					}										
					try {
						coverageArray[(readstart+r)-readClass.getCoverageStart()][0]++; 
					}
					catch(Exception e) {
					//	System.out.println(readClass.getCoverageStart());
						e.printStackTrace();
					}
					if(coverageArray[(samRecord.getUnclippedStart()-readClass.getCoverageStart())+r][0] > readClass.getMaxcoverage()) {
						readClass.setMaxcoverage(coverageArray[(samRecord.getUnclippedStart()-readClass.getCoverageStart())+r][0]);
					}
				}
				catch(Exception e) {
				//	System.out.println(splitIndex.readSeqStart +" " +splitIndex.readSequence.length());
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
					
					break;
				}
			}
		
		}
	}
	}
	catch(Exception e) {
		
		e.printStackTrace();
		return null;
	}
	return mismatches;
}
	/*
public class ReadAdder extends SwingWorker<String, Object> {
		
		boolean end = false;
		String chrom;		
		Reads readClass;
		ReadBuffer buffer;
		public ReadAdder(String chrom, Reads readClass, ReadBuffer buffer) {
			this.chrom = chrom;			
			this.readClass = readClass;
			this.buffer = buffer;
		}
		
		protected String doInBackground() {
			
			while(true) {
				try {	
					
					if(buffer.prev != null ) {						
						readSam(chrom, readClass,buffer.record);	
					//	buffer.prev.next = null;
						buffer.prev = null;
					}	
					if(buffer.next != null) {
						buffer = buffer.next;
					}
					
					if(end && buffer.next == null) {
						
						buffer.prev = null;
						break;
					}
				}
				catch(Exception e) {
					
					e.printStackTrace();
					end = true;
					break;
				}
				
				}
				Draw.updateReads = true;
				buffer.prev = null;
				buffer = null;
				
				return "";
			}

}
	*/
void readBED(File[] files) {

	try {	  

	 File bedfile;  
	 
	 for(int i = 0; i<files.length; i++) {			
		
		  bedfile = files[i];		 
		 
		  if(!checkIndex(files[i])) {
			
   			 Main.drawCanvas.loadingtext = "Indexing " +files[i].getName();
   			 if(files[i].getName().endsWith(".gz")) {
   				
   				 MethodLibrary.createBEDIndex(files[i]);
   			 }
   			 else {
   				 MethodLibrary.createBEDIndex2(files[i]);
   			 }	   			
	  	  }
		  if(!checkIndex(files[i])) {			  
			
		      Main.showError("Could not open file. Is the file sorted?", "Note");
			  continue;
		  }
		 	
	      BedTrack addTrack = new BedTrack(bedfile,Main.bedCanvas.bedTrack.size());		
	      if(bedfile.getName().endsWith("tsv.gz")) {
				 addTrack.getZerobased().setSelected(false);
				 addTrack.iszerobased = 1;
				 addTrack.getSelector().frame.setVisible(true);				
			 }	
	      Main.bedCanvas.bedTrack.add(addTrack);
	      Main.bedCanvas.trackDivider.add(0.0);
	      addTrack.minvalue = 0;
	      
	      if(bedfile.length() / 1048576 < Settings.settings.get("bigFile")) {
	    	  addTrack.small = true;		    	
	    	  Main.bedCanvas.getBEDfeatures(addTrack, 1, Main.drawCanvas.splits.get(0).chromEnd);		    	 
	      }	
	      else if(bedfile.getName().toLowerCase().endsWith(".bedgraph") || bedfile.getName().toLowerCase().endsWith(".bedgraph.gz")) {	    	  
	      
	    	  Main.bedCanvas.pressGraph(addTrack);
	    	  addTrack.getSelectorButton().setVisible(false);
	      }
			 else if(bedfile.getName().toLowerCase().endsWith("bigwig") || bedfile.getName().toLowerCase().endsWith("bw")) {
				 addTrack.small = true;
				 Main.bedCanvas.pressGraph(addTrack);
		    	  addTrack.getSelectorButton().setVisible(false);
		    	  Main.bedCanvas.getBEDfeatures(addTrack, (int)Main.drawCanvas.splits.get(0).start, (int)Main.drawCanvas.splits.get(0).end);
			 }
	      else {	    	 
	    	  setBedTrack(addTrack);
	      }   
	     
		 setTable(addTrack);
		
	 }		
	 if(Main.bedCanvas.bedTrack.size() == 0) {
		 return;
	 }
	 if(!Main.trackPane.isVisible()) {
		  Main.trackPane.setVisible(true);			 
		  Main.varpane.setDividerSize(3);	  
		  if(files.length < 5) {
			  Main.varpane.setResizeWeight(files.length *0.1);
			  Main.trackPane.setDividerLocation(files.length *0.1);
			  Main.varpane.setDividerLocation(files.length *0.1);			  
		  }
		  else {
			  Main.varpane.setResizeWeight(files.length *0.5);
			  Main.trackPane.setDividerLocation(files.length *0.5);
			  Main.varpane.setDividerLocation(files.length *0.5);	
		  }
		  Main.trackPane.revalidate();
		  Main.varpane.revalidate();
		  
	  }
	  else {
		  Main.trackPane.setDividerLocation(Main.trackPane.getDividerLocation()+70);
		  Main.varpane.setDividerLocation(Main.varpane.getDividerLocation()+70);
		  Main.varpane.revalidate();
		  Main.trackPane.revalidate();
		  
		  if(Main.controlScroll.isVisible()) {
			  Main.trackPane.setDividerSize(3);
		  }
	  }
	  Main.bedScroll.setVisible(true);
	  Main.bedCanvas.setVisible(true);
	 for(int i = 0 ; i<Main.bedCanvas.trackDivider.size(); i++) {
		 Main.bedCanvas.trackDivider.set(i, ((i+1)*(Main.varpane.getDividerLocation()/(double)Main.bedCanvas.trackDivider.size())));
	 }		
	
	}
	catch(Exception e) {			
		e.printStackTrace();
		ErrorLog.addError(e.getStackTrace());
	}		
		
}

static void setBedTrack(BedTrack addTrack) {
	try {
	 if(!addTrack.file.getName().endsWith(".gz")) {
		 return;
	 }
	  String[] split;
	  GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(addTrack.file));
	  BufferedReader reader = new BufferedReader(new InputStreamReader(gzip));
	  int count = 0;
	  if(addTrack.file.getName().endsWith(".gff.gz")) {
		  while(count < 10) {
			  if(reader.readLine().startsWith("#")) {
				  continue;
			  }
			  split = reader.readLine().split("\\t");
			  
			  
			  if(split.length > 5) {
				  if(!Double.isNaN(Double.parseDouble(split[5]))) {		    					 
					  addTrack.hasvalues = true;
				  }
				 
			  }
			  if(Main.SELEXhash.containsKey(split[2])) {								
				  addTrack.selex = true;
				  addTrack.getAffinityBox().setVisible(true);
				  addTrack.iszerobased = 1;
				  addTrack.getZerobased().setSelected(false);
			  }
			  count++;
		  }
	  }
	  else  if(addTrack.file.getName().endsWith(".bed.gz")) {
		  while(count < 10) {
			  if(reader.readLine().startsWith("#")) {
				  continue;
			  }
			  split = reader.readLine().split("\\t");			  
			  
			  if(split.length > 4) {
				  if(!Double.isNaN(Double.parseDouble(split[4]))) {		    					 
					  addTrack.hasvalues = true;
				  }
				 
			  }
			  if(Main.SELEXhash.containsKey(split[3])) {								
				  addTrack.selex = true;
				  addTrack.getAffinityBox().setVisible(true);
				  addTrack.iszerobased = 1;
				  addTrack.getZerobased().setSelected(false);
			  }
			  count++;
		  }
		  
	  }
	  else if(addTrack.file.getName().endsWith(".tsv.gz")) {
		  if(addTrack.valuecolumn != null) {
			 
		  
			  while(count < 10) {
				  if(reader.readLine().startsWith("#")) {
					  continue;
				  }
				  split = reader.readLine().split("\\t");				  
				 			
				  if(!Double.isNaN(Double.parseDouble(split[addTrack.valuecolumn]))) {	
					 
					  addTrack.hasvalues = true;
					  break;
				  }
				  	  
				 
				  count++;
			  	}
			  }
		 }
	  if(!addTrack.hasvalues) {
		  
		  addTrack.getLimitField().setVisible(false);
	  }
	  else {
		  addTrack.getLimitField().setVisible(true); 
	  }
	  gzip.close();
	  reader.close();
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

static String getInfoValue(HashMap<String, String> map, String key ) {
	
	if(map.containsKey(key)) {
		return map.get(key); 
	}
	else {
		return "-";
	}
	
}

static HashMap<String, String> makeHash(String[] line) {
	if(line == null || line.length < 9) {
		return null;
	}
	HashMap<String, String> hash = new HashMap<String, String>();
	
	String[] infosplit = line[8].split(";");
	hash.put("seqid", line[0].replace("chr", ""));
	hash.put("source", line[1]);
	hash.put("type", line[2]);
	hash.put("start", ""+(Integer.parseInt(line[3])-1));
	hash.put("end", line[4]);
	hash.put("score", line[5]);
	hash.put("strand", line[6]);
	if(line[7].equals("2")) {
		line[7] = "1";
	}
	else if(line[7].equals("1")) {
		line[7] = "2";
	}
	
	hash.put("phase", line[7]);
	
	
	for(int i = 0 ; i<infosplit.length; i++) {
		
		hash.put(infosplit[i].split("=")[0].toLowerCase(), infosplit[i].split("=")[1]);
	}
	return hash;
}
static void readGFF(File infile, String outfile,  SAMSequenceDictionary dict) {	
	BufferedReader reader = null;	
	GZIPInputStream gzip = null;
	String line= "", chrom = "-1";
	HashMap<String, String> lineHash;
	HashMap<String,Gene> genes = new HashMap<String, Gene>();
	HashMap<String,Transcript> transcripts = new HashMap<String, Transcript>();
	Gene addgene;
	Transcript addtranscript;
	
	try {
		
		if(infile.getName().endsWith(".gz")) {			
			gzip = new GZIPInputStream(new FileInputStream(infile));
			reader = new BufferedReader(new InputStreamReader(gzip));	
		}
		else {
			reader = new BufferedReader(new FileReader(infile));	
		}
	//	line = reader.readLine();
		  while((line = reader.readLine()) != null) { 
		    	
		    			    	
		    	if(line.startsWith("#")) {
		    		
		    		continue;
		    	}
		    	
		    	lineHash = makeHash(line.split("\t"));
		    	if(lineHash.get("type").startsWith("region")) {		
		    		
			    		if(line.contains("unlocalized")) {
			    			chrom = "unlocalized";
			    		}
			    		else if(lineHash.get("chromosome") != null) {
		    				chrom = lineHash.get("chromosome").replace("chr", "");	    
		    			}
		    			else if(lineHash.get("name") != null) {
		    				chrom = lineHash.get("name").replace("chr", "");	    				
		    			}	    			
		    			    		
		    		continue;
		    	}
		    	
		    	if(!lineHash.containsKey("parent")) {
		    		if(!lineHash.get("type").contains("gene")) {
		    			
		    			continue;
		    		}
		    		
		    		Gene gene = new Gene(chrom, lineHash);
		    		
		    		genes.put(getInfoValue(lineHash, "id"),gene);
		    		/*try {
		    			System.out.println(genes.get(genes.size()-1).getName());
		    		}
		    		catch(Exception e) {
		    			System.out.println(line);
		    		}*/
		    		continue;
		    	}
		    	if(genes.containsKey(lineHash.get("parent"))) {
		    		
		    		addgene = genes.get(lineHash.get("parent"));
		    		transcripts.put(getInfoValue(lineHash, "id"), new Transcript(lineHash,addgene));
		    		if(lineHash.get("type").equals("exon")) {
		    			addtranscript = transcripts.get(getInfoValue(lineHash, "id"));
		    			addtranscript.addExon(lineHash, addtranscript);
		    		}
		    		if(addgene.getDescription().equals("-")) {
		    			if(lineHash.containsKey("product")) {
		    				addgene.setDescription(lineHash.get("product"));
		    			}
		    		}
		    		continue;
		    	}
		    	if(transcripts.containsKey(lineHash.get("parent"))) {
		    	
		    		addtranscript = transcripts.get(lineHash.get("parent"));
		    		addtranscript.addExon(lineHash, addtranscript);
		    		continue;
		    	}
		    	
		  }
		
	}
	catch(Exception e) {
		System.out.println(line);
		e.printStackTrace();		
		System.exit(0);
	}
	try {
		
		Transcript transcript;
		Gene gene;
		StringBuffer exStarts, exEnds, exPhases;
	    Iterator<Map.Entry<String,Gene>> it = genes.entrySet().iterator();
	    ArrayList<String[]> geneArray = new ArrayList<String[]>();
	    
	    while (it.hasNext()) {
	        Map.Entry<String,Gene> pair = (Map.Entry<String,Gene>)it.next();
	        gene = pair.getValue();
	       
	      
	        for(int i = 0 ; i<gene.getTranscripts().size(); i++) {
	        	transcript = gene.getTranscripts().get(i);
	        	exStarts = new StringBuffer("");
	 	        exEnds = new StringBuffer("");
	 	        exPhases = new StringBuffer("");
	        	for(int e =0;e<transcript.exonArray.size(); e++) {
	        		exStarts.append(transcript.exonArray.get(e).getStart() +",");
	        		exEnds.append(transcript.exonArray.get(e).getEnd() +",");
	        		exPhases.append(transcript.exonArray.get(e).getStartPhase()+",");
	        	}
	        	
	        	String[] row = {gene.getChrom(), ""+transcript.getStart(), ""+transcript.getEnd(), gene.getName(), ""+transcript.exonArray.size(),
						MethodLibrary.getStrand(gene.getStrand()), gene.getID(), transcript.getENST(), transcript.getUniprot(),
						"-", transcript.getBiotype(), ""+transcript.getCodingStart(), ""+transcript.getCodingEnd(),
						exStarts.toString(), exEnds.toString(), exPhases.toString(), transcript.getDescription() };        	
	        	geneArray.add(row);
	        }
	       
	        it.remove();
	    }
	
	    gffSorter gffsorter = new gffSorter();
	    Collections.sort(geneArray, gffsorter);	    
	    
	    if(outfile != null) {
	    	MethodLibrary.blockCompressAndIndex(geneArray, outfile, false, dict);	  
	    }
	 
	    geneArray.clear();	   
		reader.close();
		
		if(gzip != null) {
			gzip.close();
		}
	   
	}
	catch(Exception e) {
		e.printStackTrace();
	}

}

public static class gffSorter implements Comparator<String[]> {

	public int compare(String[] o1, String[] o2) {  
		
		if(o1[0].compareTo(o2[0]) == 0) {
			if(Integer.parseInt(o1[1]) < Integer.parseInt(o2[1])) {
				return -1;
				
			}
			else if(Integer.parseInt(o1[1]) > Integer.parseInt(o2[1])) {
				return 1;
			}
			else {
				return 0;
			}			
			
		}		
		else {
			return o1[0].compareTo(o2[0]);
		}		
		  
	}
}
public static class geneSorter implements Comparator<Gene> {

	public int compare(Gene o1, Gene o2) {  
		
		if(o1.getChrom().compareTo(o2.getChrom()) == 0) {
			if(o1.getStart() < o2.getStart()) {
				return -1;				
			}
			else if(o1.getStart() > o2.getStart()) {
				return 1;
			}
			else {
				return 0;
			}			
		}		
		else {
			return -o1.getChrom().compareTo(o2.getChrom());
		}		  
	}
}
void readBED(String url, String index) {

	try {	  

	 
	  if(!Main.trackPane.isVisible()) {
		  Main.trackPane.setVisible(true);
		  Main.varpane.setDividerSize(3);	  
		  Main.varpane.setResizeWeight(0.1);
		  Main.varpane.setDividerLocation(0.1);			  
		 
	  }
	  else {
		  Main.varpane.setDividerLocation(Main.varpane.getDividerLocation()+100);
		  Main.varpane.revalidate();
		  Main.trackPane.setDividerLocation(Main.varpane.getDividerLocation()/2);			  
		  
		  if(Main.controlScroll.isVisible()) {
			  Main.trackPane.setDividerSize(3);
		  }
	  }
	  Main.bedScroll.setVisible(true);
	  Main.bedCanvas.setVisible(true);
	  BedTrack addTrack = null;
	  if(index.equals("nan")) {
		  addTrack = new BedTrack(new URL(url), null,Main.bedCanvas.bedTrack.size());	
			
	  }
	  else {
		  addTrack = new BedTrack(new URL(url), new URL(index),Main.bedCanvas.bedTrack.size());		
	  }
	 
	  //Main.bedCanvas.getBEDfeatures(addTrack, 1, Main.drawCanvas.splits.get(0).chromEnd);	
	 if(addTrack != null) {
		 Main.bedCanvas.trackDivider.add(0.0);
		 Main.bedCanvas.bedTrack.add(addTrack);
		 setTable(addTrack);
		 if(url.toLowerCase().endsWith(".bedgraph")) {
	    	  Main.bedCanvas.pressGraph(addTrack);
	      }
		 if(url.toLowerCase().endsWith("bigwig")  || url.toLowerCase().endsWith("bw")) {
			 addTrack.small = true;
			 Main.bedCanvas.getBEDfeatures(addTrack, (int)Main.drawCanvas.splits.get(0).start, (int)Main.drawCanvas.splits.get(0).end);		
			
		 }
		
		
	//	 Main.bedCanvas.getBEDfeatures(addTrack, 1, Main.drawCanvas.splits.get(0).chromEnd);	
	 }
	 
	}
	catch(Exception e) {			
		e.printStackTrace();
		ErrorLog.addError(e.getStackTrace());
	}  	 
}
static void setTable(BedTrack track) {
	
	  JScrollPane addpane = new JScrollPane();
	  addpane.setPreferredSize(new Dimension(500,400));
	  BedTable add = new BedTable((int)Main.screenSize.width, (int)Main.screenSize.height, addpane, track);
		
	  addpane.getViewport().add(add);
	  addpane.getVerticalScrollBar().addAdjustmentListener(				
			new AdjustmentListener() {					
				public void adjustmentValueChanged(AdjustmentEvent event) {
					if(VariantHandler.tabs.getSelectedIndex() > 2) {							
						VariantHandler.tables.get(VariantHandler.tabs.getSelectedIndex()-3).repaint();
					}	
					
				}										
			}
	  );
	 
	  VariantHandler.tables.add(add);	  
	  VariantHandler.tablescrolls.add(addpane);
	  track.setTable(add);
	  add.resizeTable(VariantHandler.tableScroll.getViewport().getWidth());
	  if(track.file != null) {
		  if(track.file.getName().length() > 10) {
			  VariantHandler.tabs.add(track.file.getName().substring(0, 10) +"...", addpane);
			  add.setName(track.file.getName().substring(0, 10) +"...");
		  }
		  else {
			  VariantHandler.tabs.add(track.file.getName(), addpane);
			  add.setName(track.file.getName());
		  }
	  }
	  else {
		  if(FilenameUtils.getName(track.url.getFile()).length() > 10) {
			  VariantHandler.tabs.add(FilenameUtils.getName(track.url.getFile()).substring(0, 10) +"...", addpane);
			  add.setName(FilenameUtils.getName(track.url.getFile()).substring(0, 10) +"...");
		  }
		  else {
			  VariantHandler.tabs.add(FilenameUtils.getName(track.url.getFile()), addpane);
			  add.setName(FilenameUtils.getName(track.url.getFile()));
		  }
		
	  }
	  //VariantHandler.clusterTable.addHeaderColumn(track);
	  MethodLibrary.addHeaderColumns(track);
	  VariantHandler.tabs.revalidate();
}
static void removeTable(BedTrack track) {
	try {
		
	if(VariantHandler.tables.indexOf(track.getTable()) > -1) {
		
		 VariantHandler.tabs.remove(track.getTable().tablescroll);
		 VariantHandler.tablescrolls.remove(track.getTable().tablescroll);
		 VariantHandler.tables.remove(track.getTable());
		 VariantHandler.tabs.revalidate();
		 VariantHandler.tabs.repaint();
		 
	}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

/*int binarySearch(ArrayList<int[]> list, int value, int middle) {
	
	if(middle == 0) {
		
		return 0;
	}
	if(middle == list.size()-1) {
		return list.size()-1;
	}
	if(list.get(middle-1)[1] <= value) {
		
		middle = middle/2;	
		return binarySearch(list, value, middle);
	}
	if(list.get(middle-1)[1] <= value && list.get(middle)[1] >= value) {
		return middle-1;
	}
	if(list.get(middle-1)[1] > value ) {
		
		return binarySearch(list, value, middle);
	}
	
	return -1;	
}*/
int searchIndex(ArrayList<int[]> list, int value, int low, int high) {
	
		middle = (high+low)/2;
		if(middle == list.size()-1) {
			
			return list.size();
		}
		
		if(list.get(middle)[1] <= value && list.get(middle+1)[1] >= value) {
			return middle+1;
		}
		if(list.get(middle)[1] > value) {
			
	
			return searchIndex(list, value, low, middle);
		}
		if(list.get(middle)[1] < value) {
			
			
			return searchIndex(list, value, middle, high);
		}
		return -1;
		
}

void readSam(String chrom, Reads readClass, SAMRecord samRecordBuffer, java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> mismatches) {
	try {
		
		if(readClass.getReads().isEmpty()) {
			
			addNode = new ReadNode(samRecordBuffer, readClass.sample.getComplete(), chrom, readClass.sample, splitIndex,readClass, mismatches);			
			readClass.setFirstRead(addNode);			
			startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-(readClass.readHeight+2)));					
			addNode.setRectBounds((int)((addNode.getPosition()-start)*pixel), startY, (int)(samRecordBuffer.getReadLength()*pixel), readClass.readHeight);		
			readClass.getReads().add(addNode);			
			ReadNode[] addList = new ReadNode[2];
			addList[headnode] = addNode;
			addList[tailnode] = addNode;
			readClass.getHeadAndTail().add(addList);
			readClass.setLastRead(addNode);
			
	
		}
		else {		
						
			mundane = null;
			found = false;			
			
			isDiscordant = MethodLibrary.isDiscordant(samRecordBuffer, readClass.sample.getComplete());		
			
				if(samRecordBuffer.getCigarLength() > 1) {
					if(samRecordBuffer.getCigar().getCigarElement(0).getOperator().compareTo(CigarOperator.HARD_CLIP) == 0) {
						searchPos =  samRecordBuffer.getAlignmentStart();
					}
					else {
						searchPos =  samRecordBuffer.getUnclippedStart();
					}
				}
				else {
					searchPos =  samRecordBuffer.getUnclippedStart();
				}
				
					for(int i = 0; i<readClass.getHeadAndTail().size(); i++) {						
						if(i>Settings.readDepthLimit) {
							found = true;
							mundane = null;
							addNode = null;							
							break;
						}
						if(!isDiscordant) {		
								if(readClass.getHeadAndTail().get(i)[tailnode].getPosition() +readClass.getHeadAndTail().get(i)[tailnode].getWidth() +2 < searchPos) { //.getPosition()) {
									try {										
										if(mundane == null) {
											xpos = (int)((searchPos-start)*pixel);	
											addNode = new ReadNode(samRecordBuffer, readClass.sample.getComplete(), chrom, readClass.sample, splitIndex,readClass, mismatches);	
											readClass.setLastRead(addNode);	
										}									
																				
										startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((i+1)*(readClass.readHeight+2))));							
										addNode.setRectBounds(xpos, startY, (int)(addNode.getWidth()*pixel), readClass.readHeight);
										if(i > readClass.getHeadAndTail().size()-1) {
											return;
										}
										addNode.setPrev(readClass.getHeadAndTail().get(i)[tailnode]);						
										readClass.getHeadAndTail().get(i)[tailnode].setNext(addNode);
										readClass.getHeadAndTail().get(i)[tailnode] = addNode;						
										found = true;
										break;
									}
									catch(Exception e) {
										ErrorLog.addError(e.getStackTrace());
										e.printStackTrace();
									}
									
								}
								continue;
						}						
						else  {
							
						
							
							if(readClass.getHeadAndTail().get(i)[tailnode].getPosition()+readClass.getHeadAndTail().get(i)[tailnode].getWidth() +2 < searchPos) {
								/*if(readClass.readNames.containsKey(addNode.getName())) {
									readClass.readNames.get(addNode.getName())[1] = addNode;
									
								}
								else {
									ReadNode[] nameAdd = new ReadNode[2];
									nameAdd[0] = addNode;
									nameAdd[1] = null;
									readClass.readNames.put(addNode.getName(), nameAdd);
								}*/
								
								xpos = (int)((searchPos-start)*pixel);	
								addNode = new ReadNode(samRecordBuffer, readClass.sample.getComplete(), chrom, readClass.sample, splitIndex,readClass, mismatches);	
								
								startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((i+1)*(readClass.readHeight+2))));					
								addNode.setRectBounds(xpos, startY, (int)(addNode.getWidth()*pixel), readClass.readHeight);
								addNode.setPrev(readClass.getHeadAndTail().get(i)[tailnode]);
								readClass.getHeadAndTail().get(i)[tailnode].setNext(addNode);
								readClass.getHeadAndTail().get(i)[tailnode] = addNode;		
								readClass.setLastRead(addNode);	
								found = true;
								
								break;
							}
							else if(!readClass.getHeadAndTail().get(i)[tailnode].isDiscordant()) {
								xpos = (int)((searchPos-start)*pixel);	
								addNode = new ReadNode(samRecordBuffer, readClass.sample.getComplete(), chrom, readClass.sample, splitIndex,readClass, mismatches);	
								
								mundane = readClass.getHeadAndTail().get(i)[tailnode];
								/*
								if(readClass.readNames.containsKey(addNode.getName())) {
									readClass.readNames.get(addNode.getName())[1] = addNode;
									
								}
								else {
									ReadNode[] nameAdd = new ReadNode[2];
									nameAdd[0] = addNode;
									nameAdd[1] = null;
									readClass.readNames.put(addNode.getName(), nameAdd);
								}*/
								
								startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((i+1)*(readClass.readHeight+2))));					
								addNode.setRectBounds(xpos, startY, (int)(addNode.getWidth()*pixel), readClass.readHeight);
								
								if(mundane.getPrev() != null) {
									addNode.setPrev(mundane.getPrev());
									mundane.getPrev().setNext(addNode);									
								}
								if(readClass.getHeadAndTail().get(i)[headnode].equals(mundane)) {
									readClass.getHeadAndTail().get(i)[headnode] = addNode;
								}
								readClass.getHeadAndTail().get(i)[tailnode] = addNode;
								found = false;
								readClass.setLastRead(addNode);	
								
								try {
									readClass.getReads().set(i,addNode);
								}
								catch(Exception e) {
									System.out.println(readClass.getHeadAndTail().size() +" " +readClass.getReads().size());
								}
								addNode = mundane;
								searchPos = addNode.getPosition();
								isDiscordant = false;
								addNode.setNext(null);
								addNode.setPrev(null);
								
								continue;
							}					
				
						//	continue;
						}				
					
				}
				
				if(!found) {	
					
				/*	if(addNode.isDiscordant()) {
						if(readClass.readNames.containsKey(addNode.getName())) {
							readClass.readNames.get(addNode.getName())[1] = addNode;
							}
						else {
							ReadNode[] nameAdd = new ReadNode[2];
							nameAdd[0] = addNode;
							nameAdd[1] = null;
							readClass.readNames.put(addNode.getName(), nameAdd);
						}
					}*/
					if(mundane == null) {
						xpos = (int)((searchPos-start)*pixel);	
						addNode = new ReadNode(samRecordBuffer, readClass.sample.getComplete(), chrom, readClass.sample, splitIndex,readClass, mismatches);	
						
					}
					ReadNode[] addList = new ReadNode[2];
					addList[headnode] = addNode;
					addList[tailnode] = addNode;
					addNode.setPrev(null);
					addNode.setNext(null);
					readClass.getHeadAndTail().add(addList);
					readClass.getReads().add(addNode);
					readClass.setLastRead(addNode);	
					startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((readClass.getReads().size())*(readClass.readHeight+2))));					
					addNode.setRectBounds(xpos, startY, (int)(addNode.getWidth()*pixel), readClass.readHeight);
					
				}
				
			
		}
		//preRecord = readClass.samRecordBuffer;
		
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
		
	}
	
	}

void setLevels(ReadNode node, Sample sample, Reads readClass) {
	try {
		
	while(node != null) {
		mundane = null;
		addNode = node.getPrev();
		found = false;		
		xpos = (int)((node.getPosition()-start)*pixel);
		
		for(int i = 0; i<readClass.getReads().size(); i++) {
			if(i > Settings.readDepthLimit) {
				if(node.getNext() != null) {
					node.getNext().setPrev(node.getPrev());;
				}
				if(node.getPrev() != null) {
					node.getPrev().setNext(node.getNext());
				}
				mundane = null;
				node = null;
			//	addNode = null;
				found = true;
				break;
			}
		
			if(!node.isDiscordant()) {
				
				if(node.getPosition()+node.getWidth()+2 < readClass.getHeadAndTail().get(i)[headnode].getPosition()) {
					if(i > 0) {
						
						startY = (int)(((sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((i+1)*(readClass.readHeight+2))));					
						node.setRectBounds(xpos, startY, (int)(node.getWidth()*pixel), readClass.readHeight);				
												
						if(mundane == null) {
							if(node.getPrev()!=null) {
								node.getPrev().setNext(node.getNext());
							}
							node.getNext().setPrev(node.getPrev());
						}
						
						node.setNext(readClass.getHeadAndTail().get(i)[headnode]);						
						readClass.getHeadAndTail().get(i)[headnode].setPrev(node);						
						//node.setPrev(sample.headAndTail.get(i)[headnode].getPrev());
						
						readClass.getHeadAndTail().get(i)[headnode] = node;		
						readClass.getHeadAndTail().get(i)[headnode].setPrev(null);
						
					}
					else {
					
						readClass.getHeadAndTail().get(i)[headnode] = node;
						
					}
					
					found = true;
					break;
				
				}				
				continue;					
				
			}
			else  {
				
				if(node.getPosition()+node.getWidth()+2 < readClass.getHeadAndTail().get(i)[headnode].getPosition()) {
					
					if(i > 0) {
						
						startY = (int)(((sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((i+1)*(readClass.readHeight+2))));					
						node.setRectBounds(xpos, startY, (int)(node.getWidth()*pixel), readClass.readHeight);				
					/*	if(node.getPrev()!=null) {
							node.getPrev().setNext(node.getNext());
						}
						
						node.getNext().setPrev(node.getPrev());
						*/
						node.setNext(readClass.getHeadAndTail().get(i)[headnode]);
						
						readClass.getHeadAndTail().get(i)[headnode].setPrev(node);
						node.setPrev(null);
						readClass.getHeadAndTail().get(i)[headnode] = node;					
					}
					else {
						
						readClass.getHeadAndTail().get(i)[headnode] = node;
					}
					
					found = true;
					break;
					
				}
				else if(readClass.getHeadAndTail().get(i)[headnode].isDiscordant()) {
					
					if(i==0) {
						if(node.getPrev() != null) {
							readClass.getHeadAndTail().get(i)[headnode].setPrev(node.getPrev());
							node.getPrev().setNext(readClass.getHeadAndTail().get(i)[headnode]);
							
						}						
					}					
					
					continue;
				}
				else {	
					
						mundane = readClass.getHeadAndTail().get(i)[headnode];
						
						startY = (int)(((sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((i+1)*(readClass.readHeight+2))));					
						node.setRectBounds(xpos, startY, (int)(node.getWidth()*pixel), readClass.readHeight);						
						node.setNext(mundane.getNext());	
						if(mundane.getNext() != null) {
							mundane.getNext().setPrev(node);	
						}
						mundane.setNext(null);
						mundane.setPrev(null);
						node.setPrev(null);
						readClass.getHeadAndTail().get(i)[headnode] = node;	
						readClass.getReads().set(i, node);
						node = mundane;								
						found = false;				
						continue;					
					
				}					
			}			
		}
		if(!found) {			
			
			ReadNode[] addList = new ReadNode[2];
			addList[headnode] = node;
			addList[tailnode] = node;
			readClass.getHeadAndTail().add(addList);	
			
			if(node.getPrev() != null) {
				node.getPrev().setNext(node.getNext());
				
			}			
			if(node.getNext() != null) {
				node.getNext().setPrev(node.getPrev());		
			}
			node.setNext(null);
			node.setPrev(null);
			readClass.getReads().add(node);
			startY = (int)(((sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((readClass.getReads().size())*(readClass.readHeight+2))));					
			node.setRectBounds(xpos, startY, (int)(node.getWidth()*pixel), readClass.readHeight);		
			
		}	
		node = addNode;
		
	}
	
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
	mundane = null;
	
}
void removeNonListBeds(BedNode node, int topos) {
	
	BedNode tempNode = null;
	if(node.equals(node.getTrack().getHead())) {
		node = node.getNext();
	}
	while(node != null ) {
		if(node.getPosition() >= topos) {
			
			node = null;
			break;
		}
		if(!node.inVarlist) {
			if(node.getPosition()+node.getLength()+1 < topos) {
				tempNode = node.getNext();
				node.removeNode();						
				node = tempNode;			
			}
			else {
				
				node = node.getNext();
			}
			continue;
		}
		else {
			
			node = node.getNext();
		}
		
	}
	
	tempNode = null;
	
	
}
static void removeNonListVariants() {
	
	VarNode node = FileRead.head.getNext();
	VarNode tempNode;
	Map.Entry<String, ArrayList<SampleNode>>  entry;
	int mutcount = 0;
	while(node != null ) {
		
		if(!node.inVarList || Main.drawCanvas.hideNode(node)) {			
				tempNode = node.getNext();
				node.removeNode();						
				node = tempNode;
		}
		else {	
			for(int v = 0; v<node.vars.size(); v++) {
				entry = node.vars.get(v);
				mutcount = 0;
					for(int m = 0; m<entry.getValue().size(); m++) {
						if(entry.getValue().get(m).alleles != null) {							
							break;
						}
						if(!Main.drawCanvas.hideVar(entry.getValue().get(m), entry.getKey().length() > 1)) {
							if(VariantHandler.onlyselected.isSelected()) {
								if(!entry.getValue().get(m).getSample().equals(Main.drawCanvas.selectedSample)) {
									entry.getValue().remove(m);
									m--;
									continue;
								}
							}							
							mutcount++;
						}
						else {							
							entry.getValue().remove(m);
							m--;
						}
				}
			}
			if(mutcount == 0) {
				tempNode = node.getNext();
				node.removeNode();						
				node = tempNode;
			}
			else {
				node = node.getNext();
			}
		}
	}
	node = null;
	tempNode = null;
	entry = null;
	Main.drawCanvas.current = FileRead.head;
	
	Draw.updatevars = true;
	Main.drawCanvas.repaint();
}
static void removeBedLinks() {
	
	for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
		if(Main.bedCanvas.bedTrack.get(i).used) {
			BedNode node = Main.bedCanvas.bedTrack.get(i).getHead();
			BedNode tempnode;
			try {
			while(node !=null) {
				if(node.varnodes != null) {
					for(int n = 0 ; n < node.varnodes.size(); n++) {
						
						if(Main.drawCanvas.hideNode(node.varnodes.get(n))) {
							node.varnodes.remove(n);
							n--;
							if(n == -1) {
								break;
							}
						}
						
					}
				}
				tempnode = node.getNext();
				if(!bigcalc) {
					node.putPrev(null);
					node.putNext(null);
				}
				node = tempnode;
			}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
	}
	
}

void varCalc() {
	try {
	contexts = null;
	Draw.calculateVars = false;
	Draw.variantcalculator = true;
	boolean isfreeze = VariantHandler.freeze.isSelected();
	if(VariantHandler.writetofile.isSelected()) {
		lastWriteVar = FileRead.head;
	}
	clearVariantsFromGenes();
	SplitClass split = Main.drawCanvas.splits.get(0);
	int chromcounter = 0;
	VariantHandler.outputStrings.clear();
	if(Main.drawCanvas.clusterNodes == null) {
		Main.drawCanvas.clusterNodes = new ArrayList<ClusterNode>();
		Main.drawCanvas.clusterId = 1;
	}
	else {
		Main.drawCanvas.clusterNodes.clear();
		Main.drawCanvas.clusterId = 1;
	}

	VariantHandler.table.variants = 0;
	VariantHandler.stattable.variants = 0;
	VariantHandler.clusterTable.variants = 0;
	for(int i = 0 ; i<VariantHandler.tables.size(); i++) {
		VariantHandler.tables.get(i).variants = 0;
	}

	cancelvarcount = false;
	
	VariantHandler.table.genearray.clear();
	VariantHandler.table.aminoarray.clear();	
	VariantHandler.table.controlarray = null;
	VariantHandler.stattable.sampleArray.clear();
	
	
	for(int i = 0; i<VariantHandler.tables.size(); i++) {
		VariantHandler.tables.get(i).bedarray =null;
		VariantHandler.tables.get(i).aminoarray=null;
		VariantHandler.tables.get(i).vararray=null;
		VariantHandler.tables.get(i).controlarray = null;
		VariantHandler.tables.get(i).variants = 0;		
	}	
	
	VarNode vardraw = Main.drawCanvas.current;	
	VariantHandler.freeze.setSelected(true);	
	Object[] addobject;

	for(int i=0; i<Main.samples; i++) {
		if(Main.drawCanvas.sampleList.get(i).getTabixFile() == null && !Main.drawCanvas.sampleList.get(i).calledvariants && !Main.drawCanvas.sampleList.get(i).multipart) {
			continue;
		}
		addobject = new  Object[VariantHandler.stattable.headerlengths.length];
		addobject[0] = Main.drawCanvas.sampleList.get(i);
		for(int j = 1; j<addobject.length; j++) {
			addobject[j] = 0;
		}
		Main.drawCanvas.sampleList.get(i).mutationTypes = new double[6];
		for(int j = 0; j<6; j++) {
			Main.drawCanvas.sampleList.get(i).mutationTypes[j] = 0.0;
		}
		VariantHandler.stattable.sampleArray.add(addobject);		
		Main.drawCanvas.sampleList.get(i).heterozygotes = 0;
		Main.drawCanvas.sampleList.get(i).homozygotes = 0;
		Main.drawCanvas.sampleList.get(i).varcount = 0;
		Main.drawCanvas.sampleList.get(i).indels = 0;
		Main.drawCanvas.sampleList.get(i).snvs = 0;
		Main.drawCanvas.sampleList.get(i).sitioRate = 0;
		Main.drawCanvas.sampleList.get(i).versioRate = 0;
		Main.drawCanvas.sampleList.get(i).ins = 0;
		Main.drawCanvas.sampleList.get(i).callrates = 0.0;
		Main.drawCanvas.sampleList.get(i).coding = 0;	
		
		Main.drawCanvas.sampleList.get(i).syn = 0;
		Main.drawCanvas.sampleList.get(i).nonsyn = 0;
		Main.drawCanvas.sampleList.get(i).missense = 0;
		Main.drawCanvas.sampleList.get(i).splice = 0;
		Main.drawCanvas.sampleList.get(i).nonsense = 0;
		Main.drawCanvas.sampleList.get(i).fshift = 0;
		Main.drawCanvas.sampleList.get(i).inframe = 0;
		
	}
	VariantHandler.stattable.setPreferredSize(new Dimension(VariantHandler.statsScroll.getViewport().getWidth(), (VariantHandler.stattable.sampleArray.size()+1)*15+2));
	VariantHandler.stattable.revalidate();
	VariantHandler.stattable.repaint();
	
	for(int g = 0 ; g<split.getGenes().size(); g++) {
		split.getGenes().get(g).mutations = 0;
		split.getGenes().get(g).missense = 0;
		split.getGenes().get(g).nonsense = 0;
		split.getGenes().get(g).synonymous = 0;
		split.getGenes().get(g).intronic = 0;
		split.getGenes().get(g).utr = 0;
		split.getGenes().get(g).samples.clear();
		split.getGenes().get(g).varnodes.clear();
		split.getGenes().get(g).intergenic = false;
		split.getGenes().get(g).transcriptString = new StringBuffer();		
	}
	
	if(VariantHandler.allChroms.isSelected()) {	
		if(VariantHandler.allChromsfrom.isSelected()) {
			chromcounter = Main.chromosomeDropdown.getSelectedIndex();
			Main.nothread = true;
			Main.chromosomeDropdown.setSelectedIndex(Main.chromosomeDropdown.getSelectedIndex());
			Main.nothread = false;
			vardraw = FileRead.head.getNext();		
			Main.drawCanvas.clusterNodes.clear();
			Main.drawCanvas.calcClusters(FileRead.head.getNext());			
		}		
		else if(Main.drawCanvas.splits.get(0).start != 1 || Main.chromosomeDropdown.getSelectedIndex() != 0) {
			Main.nothread = true;
			Main.chromosomeDropdown.setSelectedIndex(0);
			Main.nothread = false;
			vardraw = FileRead.head.getNext();		
			Main.drawCanvas.clusterNodes.clear();
			Main.drawCanvas.calcClusters(FileRead.head.getNext());			
		}
		else {
			vardraw = FileRead.head.getNext();
			Main.drawCanvas.clusterNodes.clear();
			Main.drawCanvas.calcClusters(FileRead.head.getNext());			
		}
	}	
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
	
	while(true) {
		
		if(cancelvarcount || !Main.drawCanvas.loading) {
			cancelFileRead();
			vardraw = null;
  			break;
  		}
		
		if(vardraw == null) {
			
			if(VariantHandler.writetofile.isSelected()) {
				 
				//flushVars(null);
				FileRead.lastpos = 0;
				Main.drawCanvas.clusterNodes.clear();				
			}
			try {
				if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {					
					VariantHandler.clusterTable.setPreferredSize(new Dimension(VariantHandler.tableScroll.getViewport().getWidth(), (Main.drawCanvas.clusterNodes.size()+1)*15));
					VariantHandler.clusterTable.revalidate();
					VariantHandler.clusterTable.repaint();
				}
				
				for(int i = 0; i<VariantHandler.stattable.sampleArray.size(); i++) {
					sample = (Sample)VariantHandler.stattable.sampleArray.get(i)[0];					
					VariantHandler.stattable.sampleArray.get(i)[1] = sample.varcount;
					VariantHandler.stattable.sampleArray.get(i)[2] = sample.snvs;
					VariantHandler.stattable.sampleArray.get(i)[3] = sample.indels;
					VariantHandler.stattable.sampleArray.get(i)[4] = sample.ins;
					VariantHandler.stattable.sampleArray.get(i)[5] = sample.coding;
					VariantHandler.stattable.sampleArray.get(i)[6] = MethodLibrary.round(sample.heterozygotes/(double)sample.homozygotes,2);
					VariantHandler.stattable.sampleArray.get(i)[7] = MethodLibrary.round((int)sample.sitioRate/(double)sample.versioRate,2);
					
					for(int j = 0; j<6; j++) {
						VariantHandler.stattable.sampleArray.get(i)[8+j] = MethodLibrary.round(sample.mutationTypes[j]/(double)sample.snvs, 2);
					}
					
					VariantHandler.stattable.sampleArray.get(i)[14] = MethodLibrary.round(sample.callrates / (double)sample.varcount,2);
					VariantHandler.stattable.repaint();
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		//	removeNonListVariants();
		//	removeBedLinks();
			
			if(VariantHandler.allChroms.isSelected()) {
				for(int i = 0; i<Main.bedCanvas.bedTrack.size(); i++) {
					removeNonListBeds( Main.bedCanvas.bedTrack.get(i).getHead(), Main.drawCanvas.splits.get(0).chromEnd);
				}
				
				if(VariantHandler.writetofile.isSelected()) {
					 FileRead.head.putNext(null);
					 nullifyVarNodes();
					clearVariantsFromGenes();
				}
				
				if(cancelvarcount || !Main.drawCanvas.loading) {
					cancelFileRead();
		  			break;
		  		}
				 FileRead.lastpos = 0;
				if(chromcounter < Main.chromosomeDropdown.getItemCount()) {
						
						chromcounter++;
						if(chromcounter == Main.chromosomeDropdown.getItemCount()) {
							break;
						}
						Main.nothread = true;					
						Main.chromosomeDropdown.setSelectedIndex(chromcounter);
						Main.nothread = false;				
						if(nobeds) {
							continue;
						}
					for(int g = 0 ; g<split.getGenes().size(); g++) {
						split.getGenes().get(g).mutations = 0;
						split.getGenes().get(g).missense = 0;
						split.getGenes().get(g).nonsense = 0;
						split.getGenes().get(g).synonymous = 0;
						split.getGenes().get(g).intronic = 0;
						split.getGenes().get(g).utr = 0;
						split.getGenes().get(g).samples.clear();
						split.getGenes().get(g).varnodes.clear();
						split.getGenes().get(g).intergenic = false;
						split.getGenes().get(g).transcriptString = new StringBuffer();
					
					}	
					Main.drawCanvas.calcClusters(FileRead.head.getNext());
					vardraw = FileRead.head.getNext();		
					FileRead.head.putNext(null);
					if(vardraw == null) {
						continue;
					}
					if(lastVar == null) {
						lastVar = FileRead.head;
					}
					if(lastWriteVar == null) {
						lastWriteVar = FileRead.head;
					}
					
					vardraw.putPrev(lastVar);
					lastVar.putNext(vardraw);
					
					VariantHandler.table.setPreferredSize(new Dimension(VariantHandler.tableScroll.getViewport().getWidth(), VariantHandler.table.getTableSize()*15));
					if(VariantHandler.tabs.getSelectedIndex() == 0) {
						VariantHandler.aminoCount.setText(VariantHandler.table.variants +" variants");
					}
					else if(VariantHandler.tabs.getSelectedIndex() == 1) {
						VariantHandler.aminoCount.setText(VariantHandler.stattable.variants +" variants");
					}
					else if(VariantHandler.tabs.getSelectedIndex() == 2) {
						VariantHandler.aminoCount.setText(VariantHandler.clusterTable.variants +" variants");
					}
					else {
						VariantHandler.aminoCount.setText(VariantHandler.tables.get(VariantHandler.tabs.getSelectedIndex()-3).variants +" variants");
					}
					VariantHandler.table.revalidate();
					VariantHandler.table.repaint();					
				}
				else {		
					Main.drawCanvas.ready("all");
					break;
				}
			}
			else {				
				break;
			}		
		}
		try {
			
		if(!VariantHandler.allChroms.isSelected()) {		
		
			if(vardraw != null && vardraw.getPosition() < Main.drawCanvas.splits.get(0).start) {				
				vardraw = vardraw.getNext();
				continue;
			}
			if(vardraw == null || vardraw.getPosition() > Main.drawCanvas.splits.get(0).end) {
				vardraw= null;
				continue;
			}
		}		
		//STATS
		
		/*if(vardraw == null) {
			
			continue;
		}*/
	//	Main.drawCanvas.loadbarAll = (int)((vardraw.getPosition()/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);
			
	//	System.out.println(vardraw.getPosition());
		//if(!VariantHandler.vcf.isSelected()) {
		
		vardraw = annotateVariant(vardraw);		
			
	/*	}
		else {			
			if(VariantHandler.writetofile.isSelected()) {
				VariantHandler.writeNodeToFile(vardraw,Main.chromosomeDropdown.getSelectedItem().toString(), output, outputgz);
			}
			else {
				annotateVariant(vardraw);		
			}
			
		}*/
	//	vardraw = vardraw.getNext();
		
		}
		catch(Exception ex) {
			ErrorLog.addError(ex.getStackTrace());
			ex.printStackTrace();
		}
	}
	
	/*if(VariantHandler.writetofile.isSelected()) {
		
		writeToFile(null, output);
	}*/
	nullifyVarNodes();
	vardraw =null;
	Draw.calculateVars = true;
	Draw.updatevars = true;
	Main.drawCanvas.repaint();
	for(int i = 0; i<VariantHandler.tables.size(); i++) {
		
		VariantHandler.tables.get(i).setPreferredSize(new Dimension(VariantHandler.tableScroll.getViewport().getWidth(), (VariantHandler.tables.get(i).getTableSize()+1)*15));
		VariantHandler.tables.get(i).revalidate();
		VariantHandler.tables.get(i).repaint();
	}

	VariantHandler.stattable.repaint();
	VariantHandler.table.setPreferredSize(new Dimension(VariantHandler.tableScroll.getViewport().getWidth(), (VariantHandler.table.getTableSize()+1)*15));
	
	if(VariantHandler.tabs.getSelectedIndex() == 0) {
		VariantHandler.aminoCount.setText(VariantHandler.table.variants +" variants");
	}
	else if(VariantHandler.tabs.getSelectedIndex() == 1) {
		VariantHandler.aminoCount.setText(VariantHandler.stattable.variants +" variants");
	}
	else if(VariantHandler.tabs.getSelectedIndex() == 2) {
	
		VariantHandler.aminoCount.setText(VariantHandler.clusterTable.variants +" variants");
	}
	else {
		VariantHandler.aminoCount.setText(VariantHandler.tables.get(VariantHandler.tabs.getSelectedIndex()-3).variants +" variants");
	}
	try {
		if(output != null) {
			
			output.close();
		}
		if(outputgz != null) {
			
			for(int i = 0 ; i<VariantHandler.outputStrings.size(); i++) {
				 outputgz.write(VariantHandler.outputStrings.get(i).getBytes());
				 
				 Feature vcf = VariantHandler.vcfCodec.decode(VariantHandler.outputStrings.get(i));
				
				 FileRead.indexCreator.addFeature(vcf, FileRead.filepointer);
				 FileRead.filepointer = outputgz.getFilePointer();
			 }
			VariantHandler.outputStrings.clear();
			outputgz.flush();
			 
		    Index index = indexCreator.finalizeIndex(outputgz.getFilePointer());
		    index.writeBasedOnFeatureFile(outFile);
			outputgz.close();
			
		}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	split = null;
	VariantHandler.table.revalidate();
	VariantHandler.table.repaint();
	if(!isfreeze) {
		VariantHandler.freeze.setSelected(false);
	}
	Main.drawCanvas.ready("all");
	if(!Main.drawCanvas.loading) {
		Draw.calculateVars = true;
	}
	Draw.variantcalculator = false;
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
	if(VariantHandler.allChroms.isSelected()) {
		Main.showError("Variant annotation ready!", "Note");
	}
	}
	catch(Exception e) {
		e.printStackTrace();
	}

	if(contexts != null) {
		
		Iterator<Entry<String, Integer[]>> iter = contexts.entrySet().iterator();
		System.out.print("#mutType");
		for(int i = 0 ; i<Main.varsamples; i++) {
			System.out.print("\t" +Main.drawCanvas.sampleList.get(i).getName());
		}
		System.out.println();
		while(iter.hasNext()) {
			Entry<String, Integer[]> entry = iter.next();			
			Integer[] samples = entry.getValue();			
			System.out.print(entry.getKey().substring(1)+">" +entry.getKey().charAt(1) +entry.getKey().charAt(0) +entry.getKey().charAt(3));
			for(int i = 0; i<samples.length; i++) {
				if(samples[i] == null) {
					System.out.print("\t0");
				}
				else {
					System.out.print("\t" +samples[i]);
				}				
			}
			System.out.println();			
		}	
	}
	boolean found = false;
	
	StringBuffer copy = new StringBuffer("");
	Clipboard clpbrd = Toolkit.getDefaultToolkit().getSystemClipboard();
	String type = "";
	for(int i = 0 ; i<VariantHandler.table.genearray.size(); i++) {
		if(VariantHandler.table.genearray.get(i).getName().equals(Main.searchField.getText().trim().toUpperCase())) {
			
			for(int j = 0; j<Main.drawCanvas.sampleList.size(); j++) {
				if(Main.drawCanvas.sampleList.get(j).multiVCF) {
					continue;
				}
				found = false;
				for(int s = 0; s < VariantHandler.table.genearray.get(i).samples.size(); s++) {
					if(VariantHandler.table.genearray.get(i).samples.get(s).getName().equals(Main.drawCanvas.sampleList.get(j).getName())) {
						found = true;
						break;
					}
				}
				if(found) {
					
					if(Main.drawCanvas.sampleList.get(j).fshift > 0) {
						type = "frameshift";
					}
					else if(Main.drawCanvas.sampleList.get(j).nonsense > 0) {
						type = "nonsense";
					}
					else if(Main.drawCanvas.sampleList.get(j).inframe > 0) {
						type = "inframe";
					}
					else if(Main.drawCanvas.sampleList.get(j).splice > 0) {
						type = "splice";
					}
					else if(Main.drawCanvas.sampleList.get(j).missense > 0) {
						type = "missense";
					}
					else if(Main.drawCanvas.sampleList.get(j).syn > 0) {
						type = "synonymous";
					}
					else {
						System.out.println(Main.drawCanvas.sampleList.get(j).getName() +" " +VariantHandler.table.genearray.get(i).getName());
					}
					copy.append(Main.drawCanvas.sampleList.get(j).getName() +"\t" +type +"\n" );
				}
				else {
					copy.append(Main.drawCanvas.sampleList.get(j).getName() +"\t0\n" );
				}
			}
			break;
		}
	}
	StringSelection stringSelection= new StringSelection(copy.toString());
	clpbrd.setContents(stringSelection, null);
}

static void clearVariantsFromBeds() {

	for(int b = 0; b< Main.bedCanvas.bedTrack.size(); b++) {
		BedNode node = Main.bedCanvas.bedTrack.get(b).getHead();
		while(node != null) {
			if(node.varnodes != null) {
				node.varnodes = null;
			}
			node = node.getNext();
		}
	}
	
}

static void nullifyVarNodes() {
	lastVar = null;
	lastWriteVar = null;
	returnnode = null;
}

void clearVariantsFromGenes() {
	SplitClass split = Main.drawCanvas.splits.get(0);

	for(int g = 0 ; g<split.getGenes().size(); g++) {
		split.getGenes().get(g).mutations = 0;
		split.getGenes().get(g).missense = 0;
		split.getGenes().get(g).nonsense = 0;
		split.getGenes().get(g).synonymous = 0;
		split.getGenes().get(g).intronic = 0;
		split.getGenes().get(g).utr = 0;
		split.getGenes().get(g).samples.clear();
		split.getGenes().get(g).varnodes.clear();
		split.getGenes().get(g).intergenic = false;
		split.getGenes().get(g).transcriptString = new StringBuffer();
	
	}
	split = null;
}

void varCalcBig() {
	bigcalc = true;
	boolean isfreeze = VariantHandler.freeze.isSelected();
	Draw.calculateVars = false;
	Draw.variantcalculator = true;
	int adder = Settings.windowSize, presearchpos = 1;
	clearVariantsFromGenes();
	int chromcounter = 0;
	clearVariantsFromBeds();
	if(Main.drawCanvas.clusterNodes == null) {
		Main.drawCanvas.clusterNodes = new ArrayList<ClusterNode>();
		Main.drawCanvas.clusterId = 1;
	}
	else {
		Main.drawCanvas.clusterNodes.clear();
		Main.drawCanvas.clusterId = 1;
	}
	
	VariantHandler.table.variants = 0;
	VariantHandler.stattable.variants = 0;
	VariantHandler.clusterTable.variants = 0;
	for(int i = 0 ; i<VariantHandler.tables.size(); i++) {
		VariantHandler.tables.get(i).variants = 0;
	}
	
	cancelvarcount = false;
	VariantHandler.table.genearray.clear();
	VariantHandler.table.aminoarray.clear();	
	VariantHandler.table.controlarray = null;
	VariantHandler.stattable.sampleArray.clear();


	for(int i = 0; i<VariantHandler.tables.size(); i++) {
		VariantHandler.tables.get(i).bedarray= null;
		VariantHandler.tables.get(i).aminoarray= null;
		VariantHandler.tables.get(i).vararray= null;
		VariantHandler.tables.get(i).controlarray= null;
		VariantHandler.tables.get(i).variants = 0;
		
	}	
	head.putNext(null);
	VarNode vardraw;
	
	VariantHandler.freeze.setSelected(true);
	
	Object[] addobject;
	
	for(int i=0; i<Main.samples; i++) {
		if(Main.drawCanvas.sampleList.get(i).getTabixFile() == null && !Main.drawCanvas.sampleList.get(i).multipart) {
			continue;
		}
		addobject = new  Object[VariantHandler.stattable.headerlengths.length];
		addobject[0] = Main.drawCanvas.sampleList.get(i);
		for(int j = 1; j<addobject.length; j++) {
			addobject[j] = 0;
		}
		Main.drawCanvas.sampleList.get(i).mutationTypes = new double[6];
		for(int j = 0; j<6; j++) {
			Main.drawCanvas.sampleList.get(i).mutationTypes[j] = 0.0;
		}
		VariantHandler.stattable.sampleArray.add(addobject);
		
		Main.drawCanvas.sampleList.get(i).heterozygotes = 0;
		Main.drawCanvas.sampleList.get(i).homozygotes = 0;
		Main.drawCanvas.sampleList.get(i).varcount = 0;
		Main.drawCanvas.sampleList.get(i).indels = 0;
		Main.drawCanvas.sampleList.get(i).snvs = 0;
		Main.drawCanvas.sampleList.get(i).sitioRate = 0;
		Main.drawCanvas.sampleList.get(i).versioRate = 0;
		Main.drawCanvas.sampleList.get(i).ins = 0;
		Main.drawCanvas.sampleList.get(i).callrates = 0.0;
		Main.drawCanvas.sampleList.get(i).coding = 0;
		
		Main.drawCanvas.sampleList.get(i).syn = 0;
		Main.drawCanvas.sampleList.get(i).nonsyn = 0;
		Main.drawCanvas.sampleList.get(i).missense = 0;
		Main.drawCanvas.sampleList.get(i).splice = 0;
		Main.drawCanvas.sampleList.get(i).nonsense = 0;
		Main.drawCanvas.sampleList.get(i).fshift = 0;
		Main.drawCanvas.sampleList.get(i).inframe = 0;
		
	}
	VariantHandler.stattable.setPreferredSize(new Dimension(VariantHandler.statsScroll.getViewport().getWidth(), (VariantHandler.stattable.sampleArray.size()+1)*15+2));
	VariantHandler.stattable.revalidate();
	VariantHandler.stattable.repaint();
	clearVariantsFromGenes();
	
	if(VariantHandler.allChroms.isSelected()) {		
		if(VariantHandler.allChromsfrom.isSelected()) {
			chromcounter = Main.chromosomeDropdown.getSelectedIndex();
			Main.nothread = true;
			search = true;
			searchStart = (int)Main.drawCanvas.splits.get(0).start;
			searchEnd = (int)Main.drawCanvas.splits.get(0).start+Settings.windowSize;
			Main.chromosomeDropdown.setSelectedIndex(Main.chromosomeDropdown.getSelectedIndex());
			Main.nothread = false;
				
			Main.drawCanvas.clusterNodes.clear();
			Main.drawCanvas.calcClusters(FileRead.head.getNext());	
			vardraw = FileRead.head.getNext();	
		}		
		else if(Main.drawCanvas.splits.get(0).start > 10 || Main.chromosomeDropdown.getSelectedIndex() != 0) {
			Main.nothread = true;
			search = true;			
			Main.chromosomeDropdown.setSelectedIndex(0);
			Main.nothread = false;
			vardraw = FileRead.head.getNext();		
			Main.drawCanvas.clusterNodes.clear();
			Main.drawCanvas.calcClusters(FileRead.head.getNext());
			
		}
		else {			
			searchStart = 0;
			searchEnd = searchStart+Settings.windowSize;			
			getVariantWindow(Main.drawCanvas.splits.get(0).chrom,searchStart,searchEnd);
			Main.drawCanvas.clusterNodes.clear();
			Main.drawCanvas.calcClusters(FileRead.head.getNext());
			vardraw = FileRead.head.getNext();
		}
		
	}	
	else {
	
		searchStart = (int)Main.drawCanvas.splits.get(0).start;
		
		searchEnd = getNextIntergenic(searchStart+adder);
		if(searchEnd >  (int)Main.drawCanvas.splits.get(0).end) {
			searchEnd = (int)Main.drawCanvas.splits.get(0).end;
					
		}
		
		
		getVariantWindow(Main.drawCanvas.splits.get(0).chrom,searchStart,searchEnd);
		vardraw = FileRead.head.getNext();
		
	}
	presearchpos = searchEnd;
	lastVar = FileRead.head;
	
	if(VariantHandler.writetofile.isSelected()) {
	
		lastWriteVar = FileRead.head;
	}
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
	while(true) {
		
		if(cancelvarcount || cancelfileread || !Main.drawCanvas.loading) {
			cancelFileRead();
			nullifyVarNodes();
			vardraw = null;
  			break;
  		}
		if(vardraw == null) {
			if(VariantHandler.writetofile.isSelected()) {
				
				if(VariantHandler.geneSlider.getValue() > 1) {
					 
					flushVars(vardraw);
					 
				 }
				else {
					clearVariantsFromGenes();
				}
				lastVar.putPrev(null);
				lastWriteVar.putPrev(null);
			}
			try {
				if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0 && Main.drawCanvas.clusterNodes.size() > 0) {					
					Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).width = Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.get(Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.size()-1).getPosition() -Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.get(0).getPosition() +1;
					VariantHandler.clusterTable.setPreferredSize(new Dimension(VariantHandler.tableScroll.getViewport().getWidth(), (Main.drawCanvas.clusterNodes.size()+1)*15));
					VariantHandler.clusterTable.revalidate();
					VariantHandler.clusterTable.repaint();
				}
				
				for(int i = 0; i<VariantHandler.stattable.sampleArray.size(); i++) {
					sample = (Sample)VariantHandler.stattable.sampleArray.get(i)[0];					
					VariantHandler.stattable.sampleArray.get(i)[1] = sample.varcount;
					VariantHandler.stattable.sampleArray.get(i)[2] = sample.snvs;
					VariantHandler.stattable.sampleArray.get(i)[3] = sample.indels;
					VariantHandler.stattable.sampleArray.get(i)[4] = sample.ins;
					VariantHandler.stattable.sampleArray.get(i)[5] = sample.coding;
					VariantHandler.stattable.sampleArray.get(i)[6] = MethodLibrary.round(sample.heterozygotes/(double)sample.homozygotes,2);
					VariantHandler.stattable.sampleArray.get(i)[7] = MethodLibrary.round((int)sample.sitioRate/(double)sample.versioRate,2);
					
					for(int j = 0; j<6; j++) {
						VariantHandler.stattable.sampleArray.get(i)[8+j] = MethodLibrary.round(sample.mutationTypes[j]/(double)sample.snvs, 2);
					}
					
					VariantHandler.stattable.sampleArray.get(i)[14] = MethodLibrary.round(sample.callrates / (double)sample.varcount,2);
					VariantHandler.stattable.repaint();
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			if(VariantHandler.allChroms.isSelected() && searchEnd > Main.drawCanvas.splits.get(0).chromEnd) {
				
				for(int i = 0; i<Main.bedCanvas.bedTrack.size(); i++) {
					removeNonListBeds( Main.bedCanvas.bedTrack.get(i).getHead(), Main.drawCanvas.splits.get(0).chromEnd);
				}
				
				if(cancelvarcount) {
		  			break;
		  		}
				if(Main.selectedChrom < Main.chromosomeDropdown.getItemCount()) {
					
					chromcounter++;
					if(chromcounter == Main.chromosomeDropdown.getItemCount()) {
						break;
					}
					//clearVariantsFromGenes();
					Main.nothread = true;					
					//search = true;
					
					Main.chromosomeDropdown.setSelectedIndex(chromcounter);
					searchStart = 1;
					searchEnd = adder;
					
					//search = false;
					Main.nothread = false;
					
					presearchpos = searchEnd;
					getVariantWindow(Main.drawCanvas.splits.get(0).chrom,searchStart,searchEnd);
					
					Main.drawCanvas.calcClusters(FileRead.head.getNext());
					vardraw = head.getNext();
					FileRead.head.putNext(null);
					if(vardraw == null) {
						continue;
					}
					
					vardraw.putPrev(lastVar);
					lastVar.putNext(vardraw);
					
					VariantHandler.table.setPreferredSize(new Dimension(VariantHandler.tableScroll.getViewport().getWidth(), VariantHandler.table.getTableSize()*15));
					if(VariantHandler.tabs.getSelectedIndex() == 0) {
						VariantHandler.aminoCount.setText(VariantHandler.table.variants +" variants");
					}
					else if(VariantHandler.tabs.getSelectedIndex() == 1) {
						VariantHandler.aminoCount.setText(VariantHandler.stattable.variants +" variants");
					}
					else if(VariantHandler.tabs.getSelectedIndex() == 2) {
						VariantHandler.aminoCount.setText(VariantHandler.clusterTable.variants +" variants");
					}
					else {
						VariantHandler.aminoCount.setText(VariantHandler.tables.get(VariantHandler.tabs.getSelectedIndex()-3).variants +" variants");
					}
					VariantHandler.table.revalidate();
					VariantHandler.table.repaint();		
					
				}
				
			}
			else {						
				
				
				VariantHandler.table.setPreferredSize(new Dimension(VariantHandler.tableScroll.getViewport().getWidth(), VariantHandler.table.getTableSize()*15));
				if(VariantHandler.tabs.getSelectedIndex() == 0) {
					VariantHandler.aminoCount.setText(VariantHandler.table.variants +" variants");
				}
				else if(VariantHandler.tabs.getSelectedIndex() == 1) {
					VariantHandler.aminoCount.setText(VariantHandler.stattable.variants +" variants");
				}
				else if(VariantHandler.tabs.getSelectedIndex() == 2) {
					VariantHandler.aminoCount.setText(VariantHandler.clusterTable.variants +" variants");
				}
				else {
					VariantHandler.aminoCount.setText(VariantHandler.tables.get(VariantHandler.tabs.getSelectedIndex()-3).variants +" variants");
				}
				
				VariantHandler.table.revalidate();
				VariantHandler.table.repaint();	
				
				Main.chromDraw.varnode = null;
				Main.chromDraw.vardraw = null;
			
				if(searchEnd >= (int)Main.drawCanvas.splits.get(0).end) {
					break;
				}
				
				searchStart = presearchpos;
				searchEnd = searchStart + adder; //getNextIntergenic(presearchpos+adder);
				presearchpos = searchEnd;
				for(int i = 0; i<Main.bedCanvas.bedTrack.size(); i++) {					
					removeNonListBeds( Main.bedCanvas.bedTrack.get(i).getHead(), searchStart);					
				}
				getVariantWindow(Main.drawCanvas.splits.get(0).chrom,searchStart,searchEnd);
				
				vardraw = head.getNext();
				head.putNext(null);
				if(vardraw == null) {
					continue;
				}
				vardraw.putPrev(lastVar);
				lastVar.putNext(vardraw);
			}		
		}
		try {
			
		if(!VariantHandler.allChroms.isSelected()) {		
		
			if(vardraw != null && vardraw.getPosition() < Main.drawCanvas.splits.get(0).start) {				
				vardraw = vardraw.getNext();
			}
			if(vardraw == null || vardraw.getPosition() > Main.drawCanvas.splits.get(0).end) {
				vardraw= null;
				continue;
			}
		}
		
		
		
		//STATS
	/*	if(vardraw == null) {
			continue;
		}
		*/
	//	Main.drawCanvas.loadbarAll = (int)((vardraw.getPosition()/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);
		
			
		
	/*	coding = false;
		if(vardraw.getExons() != null) {
			for(int e = 0; e < vardraw.getExons().size(); e++) {
				if(vardraw.getPosition()+1 > vardraw.getExons().get(e).getTranscript().getCodingStart() && vardraw.getPosition()+1 < vardraw.getExons().get(e).getTranscript().getCodingEnd()) {
					if(vardraw.getPosition()+1 >= vardraw.getExons().get(e).getStart()-2 && vardraw.getPosition()+1 < vardraw.getExons().get(e).getEnd()+2) {
						coding = true;
						break;
					}
				}
				
			}
		}
		*/
		
		
		//if(!VariantHandler.vcf.isSelected()) {
		
			vardraw = annotateVariant(vardraw);		
		/*}
		else {			
			if(VariantHandler.writetofile.isSelected()) {
				VariantHandler.writeNodeToFile(vardraw,Main.chromosomeDropdown.getSelectedItem().toString(), output, outputgz);
			}
			else {
				annotateVariant(vardraw);		
			}
		}
	*/
		
		//vardraw = vardraw.getNext();
		
		}
		catch(Exception ex) {
			ErrorLog.addError(ex.getStackTrace());
			ex.printStackTrace();
		}
	}
	Draw.calculateVars = true;
	Draw.updatevars = true;
	Main.drawCanvas.repaint();
	vardraw = null;
	for(int i = 0; i<VariantHandler.tables.size(); i++) {
		VariantHandler.tables.get(i).setPreferredSize(new Dimension(VariantHandler.tableScroll.getViewport().getWidth(), (VariantHandler.tables.get(i).getTableSize()+1)*15));
		VariantHandler.tables.get(i).revalidate();
		VariantHandler.tables.get(i).repaint();
	}
	VariantHandler.stattable.repaint();
	VariantHandler.table.setPreferredSize(new Dimension(VariantHandler.tableScroll.getViewport().getWidth(), (VariantHandler.table.getTableSize()+1)*15));
	
	if(VariantHandler.tabs.getSelectedIndex() == 0) {
		VariantHandler.aminoCount.setText(VariantHandler.table.variants +" variants");
	}
	else if(VariantHandler.tabs.getSelectedIndex() == 1) {
		VariantHandler.aminoCount.setText(VariantHandler.stattable.variants +" variants");
	}
	else if(VariantHandler.tabs.getSelectedIndex() == 2) {
	
		VariantHandler.aminoCount.setText(VariantHandler.clusterTable.variants +" variants");
	}
	else {
		VariantHandler.aminoCount.setText(VariantHandler.tables.get(VariantHandler.tabs.getSelectedIndex()-3).variants +" variants");
	}
	try {
		if(output != null) {
			output.close();
		}
		if(outputgz != null) {
			for(int i = 0 ; i<VariantHandler.outputStrings.size(); i++) {
				 outputgz.write(VariantHandler.outputStrings.get(i).getBytes());
				 
				 Feature vcf = VariantHandler.vcfCodec.decode(VariantHandler.outputStrings.get(i));
				 FileRead.indexCreator.addFeature(vcf, FileRead.filepointer);
				 FileRead.filepointer = outputgz.getFilePointer();
				 
			 }
			
			VariantHandler.outputStrings.clear();
			outputgz.flush();
			Index index = indexCreator.finalizeIndex(outputgz.getFilePointer());
		    index.writeBasedOnFeatureFile(outFile);
			outputgz.close();
		}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	VariantHandler.table.revalidate();
	VariantHandler.table.repaint();
	if(!isfreeze) {
		VariantHandler.freeze.setSelected(false);
	}
	
	Main.drawCanvas.current = FileRead.head;
	if(!Main.drawCanvas.loading) {
		Draw.calculateVars = true;
	}
	Main.drawCanvas.ready("all");
	nullifyVarNodes();
	bigcalc = false;
	Draw.variantcalculator = false;
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
	if(VariantHandler.allChroms.isSelected()) {
		Main.showError("Variant annotation ready!", "Note");
		 	
	}
}

int getNextIntergenic(int pos) {
	int returnpos = pos;
	int pointer = 0, gene = -1;
	
	while(true) {
		gene = MethodLibrary.getRegion(returnpos, Main.drawCanvas.splits.get(0), pointer);
		
		if(gene == -1) {
			break;
		}
		else {
			if(returnpos <  Main.drawCanvas.splits.get(0).getGenes().get(gene).getEnd()+1) {
				returnpos = Main.drawCanvas.splits.get(0).getGenes().get(gene).getEnd()+1;
			}
			pointer = gene;
		}
	}
	
	return returnpos;
}


void writeToFile(VarNode vardraw, BufferedWriter output, BlockCompressedOutputStream outputgz) {
	
	if(VariantHandler.table.genearray.size() > 0) {
		
		if(VariantHandler.tsv.isSelected() || VariantHandler.compactTsv.isSelected()) {
			
			for(int i = 0 ; i<VariantHandler.table.genearray.size(); i++) {
				if(vardraw == null || vardraw.getPosition() > VariantHandler.table.genearray.get(i).getEnd()) {
					
					VariantHandler.writeTranscriptToFile(VariantHandler.table.genearray.get(i), output);							
						
					if(VariantHandler.allChroms.isSelected() || bigcalc) {
						for(int s = 0; s<VariantHandler.table.genearray.get(i).varnodes.size(); s++) {
							if(VariantHandler.table.genearray.get(i).varnodes.get(s).getExons() != null) {
								for(int e = 0 ; e<VariantHandler.table.genearray.get(i).varnodes.get(s).getExons().size(); e++) {
									if(VariantHandler.table.genearray.get(i).varnodes.get(s).getExons().get(e).getTranscript().equals(VariantHandler.table.genearray.get(i))) {
										VariantHandler.table.genearray.get(i).varnodes.get(s).getExons().remove(e);
										e--;
									}
								}						
							}
							if(VariantHandler.table.genearray.get(i).varnodes.get(s).getTranscripts() != null) {
								for(int e = 0 ; e<VariantHandler.table.genearray.get(i).varnodes.get(s).getTranscripts().size(); e++) {
									if(VariantHandler.table.genearray.get(i).varnodes.get(s).getTranscripts().equals(VariantHandler.table.genearray.get(i))) {
										VariantHandler.table.genearray.get(i).varnodes.get(s).getTranscripts().remove(e);
										e--;
									}
									
								}
							}
						}
					}
					VariantHandler.table.genearray.get(i).samples.clear();
					VariantHandler.table.genearray.get(i).varnodes.clear();
					VariantHandler.table.genearray.remove(i);
					i--;
				}				
			}
		}
		else if(VariantHandler.vcf.isSelected()){
			
			/* ArrayList<VarNode> nodes = new ArrayList<VarNode>();			
			 String lastChrom;
			 
			 for(int v=0; v<VariantHandler.table.genearray.get(0).varnodes.size(); v++) {		    					
				 if(!nodes.contains(VariantHandler.table.genearray.get(0).varnodes.get(v))) {		    						
					 nodes.add(VariantHandler.table.genearray.get(0).varnodes.get(v));	
					
				 }		    					
			 }	  
			 if(lastpos < vardraw.getPosition()) {
				 lastpos = vardraw.getPosition();
			 }
			 else {
				 return;
			 }

			 lastChrom = Main.chromosomeDropdown.getSelectedItem().toString();
			 */
			/* if(bigcalc) {
				 Main.drawCanvas.loadbarAll = (int)((lastpos/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);
			 }
			 else {
				 Main.drawCanvas.loadbarAll = (int)((lastpos/(double)(Main.drawCanvas.variantsEnd))*100);
			 }
			 Main.drawCanvas.loadBarSample = Main.drawCanvas.loadbarAll;	
			 */
			if(VariantHandler.geneSlider.getValue() > 1) {
				
			}
			else {				
				
				VariantHandler.writeNodeToFile(vardraw,vardraw.getChrom(), output, outputgz);	
				lastWriteVar = vardraw;
			}			 
			
			 for(int i = 0 ; i< VariantHandler.table.genearray.size() ; i++) {
				 if(vardraw == null || VariantHandler.table.genearray.get(i).getEnd() < vardraw.getPosition()) {					 
					 
					 if(VariantHandler.allChroms.isSelected() || bigcalc) {
							for(int s = 0; s<VariantHandler.table.genearray.get(i).varnodes.size(); s++) {
								if(VariantHandler.table.genearray.get(i).varnodes.get(s).getExons() != null) {
									for(int e = 0 ; e<VariantHandler.table.genearray.get(i).varnodes.get(s).getExons().size(); e++) {
										if(VariantHandler.table.genearray.get(i).varnodes.get(s).getExons().get(e).getTranscript().equals(VariantHandler.table.genearray.get(i))) {
											VariantHandler.table.genearray.get(i).varnodes.get(s).getExons().remove(e);
											e--;
										}
									}						
								}
								if(VariantHandler.table.genearray.get(i).varnodes.get(s).getTranscripts() != null) {
									for(int e = 0 ; e<VariantHandler.table.genearray.get(i).varnodes.get(s).getTranscripts().size(); e++) {
										if(VariantHandler.table.genearray.get(i).varnodes.get(s).getTranscripts().equals(VariantHandler.table.genearray.get(i))) {
											VariantHandler.table.genearray.get(i).varnodes.get(s).getTranscripts().remove(e);
											e--;
										}										
									}
								}
							}
						}					 
						 if(VariantHandler.geneSlider.getValue() == 1) {
							VariantHandler.table.genearray.get(i).samples.clear();
						 }
						VariantHandler.table.genearray.get(i).varnodes.clear();
						VariantHandler.table.genearray.remove(i);
						i--;
				 }
			 }		
			 
			 if(VariantHandler.geneSlider.getValue() > 1) {
				 
				 if(VariantHandler.table.genearray.size() == 0) {
					
					flushVars(vardraw);
				 }
			 }			 
		}
		
	}	
}

void flushVars(VarNode vardraw) {		
		if(VariantHandler.vcf.isSelected() || VariantHandler.oncodrive.isSelected()) {
			if(lastWriteVar != null) {
				if(lastWriteVar.getNext()== null) {
					found = false;
					
					if(lastWriteVar.isInGene()) {
						if(lastWriteVar.getExons() != null) {
							for(int i = 0 ; i<lastWriteVar.getExons().size(); i++) {
								if(lastWriteVar.getExons().get(i).getTranscript().getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
									
									found = true;
								}
								else {
									lastWriteVar.getExons().remove(i);
									i--;
								}
							}
						}
						else {
							for(int i = 0 ; i<lastWriteVar.getTranscripts().size(); i++) {
								if(lastWriteVar.getTranscripts().get(i).getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
									found = true;
								}
								else {
									lastWriteVar.getTranscripts().remove(i);
									i--;
								}
							}
						}	
					}
					if(found) {					
						VariantHandler.writeNodeToFile(lastWriteVar,lastWriteVar.getChrom(), output, outputgz);					
					}
				}
				else {
					lastWriteVar = lastWriteVar.getNext();
				}
			}
			if(vardraw == null) {
				
				if(VariantHandler.table.genearray.size() > 0) {
					
					while(lastWriteVar != null) {	
						
							found = false;
							if(lastWriteVar.isInGene()) {
								if(lastWriteVar.getExons() != null) {
									for(int i = 0 ; i<lastWriteVar.getExons().size(); i++) {
										if(lastWriteVar.getExons().get(i).getTranscript().getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
											
											found = true;
											
										}
										else {
											lastWriteVar.getExons().remove(i);
											i--;
										}
									}
								}
								else {
									for(int i = 0 ; i<lastWriteVar.getTranscripts().size(); i++) {
										if(lastWriteVar.getTranscripts().get(i).getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
											found = true;
											
										}
										
										else {
											lastWriteVar.getTranscripts().remove(i);
											i--;
										}
									}
								}	
							}
							if(found) {
								
								VariantHandler.writeNodeToFile(lastWriteVar,lastWriteVar.getChrom(), output, outputgz);		
								
							}
							/*if(lastWriteVar.getNext() == null) {
								break;
							}*/
							
							lastWriteVar = lastWriteVar.getNext();
					}
				}
			}
			else {
				
			 while(!lastWriteVar.equals(vardraw)) {
				 
				 found = false;
					if(lastWriteVar.isInGene()) {
						if(lastWriteVar.getExons() != null) {
							for(int i = 0 ; i<lastWriteVar.getExons().size(); i++) {
								if(lastWriteVar.getExons().get(i).getTranscript().getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
									
									found = true;
								}
								else {
									lastWriteVar.getExons().remove(i);
									i--;
								}
							}
						}
						else {
							for(int i = 0 ; i<lastWriteVar.getTranscripts().size(); i++) {
								if(lastWriteVar.getTranscripts().get(i).getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
									found = true;
								}
								else {
									lastWriteVar.getTranscripts().remove(i);
									i--;
								}
							}
						}	
					}
					if(found) {
						
						VariantHandler.writeNodeToFile(lastWriteVar,lastWriteVar.getChrom(), output, outputgz);		
						
					}
				 
					lastWriteVar = lastWriteVar.getNext();
			}
			
			}
		 
		}
		else {
			if(vardraw == null && VariantHandler.table.genearray.size() > 0) {
				for(int i =0;i< VariantHandler.table.genearray.size(); i++) {
					VariantHandler.writeTranscriptToFile(VariantHandler.table.genearray.get(i), output);
				}
				VariantHandler.table.genearray.clear();
			}
		}
}

VarNode annotateVariant(VarNode vardraw) {
	//System.out.println(vardraw.getPosition());
	if(Main.drawCanvas.hideNode(vardraw)) {	
		
	    returnnode = vardraw.getNext();
	    if(VariantHandler.allChroms.isSelected()) {
	    	vardraw.removeNode();		
	    }
		return returnnode;					
	}	
	
	Map.Entry<String, ArrayList<SampleNode>> entry = null;
	String base = null, amino = null;
	int pretrack = -1, preI = -1;
	
	Main.drawCanvas.loadbarAll = (int)((vardraw.getPosition()/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);	
	Main.drawCanvas.loadBarSample = Main.drawCanvas.loadbarAll;
	
	for(int v = 0; v<vardraw.vars.size(); v++) {
	 	entry = vardraw.vars.get(v);
	 	
		if(Main.drawCanvas.hideNodeVar(vardraw, entry)) {
			if(VariantHandler.allChroms.isSelected()) {
				vardraw.vars.remove(v);
				v--;
			}
			continue;
		}			
		
		base = entry.getKey();
		mutcount = 0;
	
		for(int m = 0; m<entry.getValue().size(); m++) {
				if(entry.getValue().get(m).alleles != null) {
					break;
				}
			
				if(Main.drawCanvas.hideVar(entry.getValue().get(m), entry.getKey().length() > 1)) {
					if(VariantHandler.allChroms.isSelected()) {
						entry.getValue().remove(m);
						m--;
					}
					continue;
				}
				
				sample = entry.getValue().get(m).getSample();
				if(VariantHandler.onlyselected.isSelected()) {
					if(!sample.equals(Main.drawCanvas.selectedSample)) {
						continue;
					}
				}
				sample.varcount++;
				mutcount++;
				VariantHandler.stattable.variants++;
				if(vardraw.coding) {
					sample.coding++;
				}
				sample.callrates += entry.getValue().get(m).getAlleleFraction();
				
				if(entry.getKey().length() > 1) {						
					if(entry.getKey().contains("del")) {
						sample.indels++;
					}
					else {
						sample.ins++;
					}						
				}
			
				if(entry.getValue().get(m).isHomozygous()) {
					sample.homozygotes++;
				}
				else {
					sample.heterozygotes++;
				}							
			
				if(entry.getKey().length() < 2) {
					if(VariantHandler.onlyStats.isSelected()) {
						String context = Main.chromDraw.getSeq(Main.drawCanvas.splits.get(0).chrom,vardraw.getPosition()-1, vardraw.getPosition()+2, Main.referenceFile).toString();
						
						if(contexts == null) {
							contexts = new HashMap<String, Integer[]>();
						}
					/*	if(base.equals("C") && context.equals("ATA")) {
							System.out.println(sample.getName() +" " +vardraw.getPosition());
						}*/
						basecontext = base+context;							
						
						
						if(contexts.containsKey(base+context)) {
							if(contexts.get(base+context)[sample.getIndex()] == null) {
								contexts.get(base+context)[sample.getIndex()] = 1; 
							}
							else {
								contexts.get(base+context)[sample.getIndex()]++;
							}
						}
						else {
							if(contexts.containsKey(MethodLibrary.reverseComplement(base)+MethodLibrary.reverseComplement(context))) {
								if(contexts.get(MethodLibrary.reverseComplement(base)+MethodLibrary.reverseComplement(context))[sample.getIndex()] == null) {
									contexts.get(MethodLibrary.reverseComplement(base)+MethodLibrary.reverseComplement(context))[sample.getIndex()] = 1; 
								}
								else {
									contexts.get(MethodLibrary.reverseComplement(base)+MethodLibrary.reverseComplement(context))[sample.getIndex()]++;
								}
							}
							else {
								contexts.put(base+context, new Integer[Main.samples]);
								contexts.get(base+context)[sample.getIndex()] = 1; 
							}
						}
												
					}
					
					sample.snvs++;
					try {
						
						sample.mutationTypes[Main.mutTypes.get(Main.getBase.get(vardraw.getRefBase()) +entry.getKey())]++;
					}
					catch(Exception e) {
						
						System.out.println(vardraw.getPosition() +" " +(char)vardraw.getRefBase() +entry.getKey());
						e.printStackTrace();
					}
					if(((char)vardraw.getRefBase() == 'A' && entry.getKey().equals("G")) || ((char)vardraw.getRefBase() == 'G' && entry.getKey().equals("A")) || ((char)vardraw.getRefBase() == 'C' && entry.getKey().equals("T")) || ((char)vardraw.getRefBase() == 'T' && entry.getKey().equals("C"))) {
						sample.sitioRate++;
					}
					else {
						sample.versioRate++;
					}						
				}
		}
	if(mutcount == 0) {
		continue;
	}
	// INTRONIC
	
	if(VariantHandler.intronic.isSelected() && vardraw.isInGene() && vardraw.getTranscripts() != null && vardraw.getExons() == null) {
		
		for(int t = 0; t<vardraw.getTranscripts().size(); t++) {				
		
			for(int i = 0; i<entry.getValue().size(); i++) {
				if(entry.getValue().get(i).alleles != null) {
					break;
				}
				if(Main.drawCanvas.hideVar(entry.getValue().get(i), entry.getKey().length() > 1)) {
					if(VariantHandler.allChroms.isSelected()) {
						entry.getValue().remove(i);
						i--;
					}
					continue;
				}
				if(VariantHandler.onlyselected.isSelected()) {
					if(!entry.getValue().get(i).getSample().equals(Main.drawCanvas.selectedSample)) {
						continue;
					}
				}
				if(!vardraw.getTranscripts().get(t).getGene().samples.contains(entry.getValue().get(i).getSample())) {
					vardraw.getTranscripts().get(t).getGene().samples.add(entry.getValue().get(i).getSample());
				}				
			}
			if(!VariantHandler.onlyStats.isSelected()) {
				
				if(vardraw.getTranscripts().get(t).getGene().mutations == 0 && vardraw.getTranscripts().get(t).getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
					VariantHandler.table.addEntry(vardraw.getTranscripts().get(t).getGene());			
					VariantHandler.table.revalidate();
					VariantHandler.table.repaint();	
				}
			}
			vardraw.getTranscripts().get(t).getGene().intronic+=mutcount;				
			vardraw.inVarList = true;							
		
			if(vardraw.inVarList) {
				
				if(!vardraw.getTranscripts().get(t).getGene().varnodes.contains(vardraw)) {
					if(!VariantHandler.onlyStats.isSelected()) {
						vardraw.getTranscripts().get(t).getGene().varnodes.add(vardraw);	
					}
					preI = -1;
				}
				if(v != preI) {
					VariantHandler.table.variants += mutcount;
					vardraw.getTranscripts().get(t).getGene().mutations+=mutcount;
					preI = v;
				}
				
			}
		}	
	}
	
	
	if(vardraw != null && vardraw.getExons() != null) {	
		ArrayList<Gene> calcgene = new ArrayList<Gene>();
		for(int exon = 0; exon<vardraw.getExons().size(); exon++) {	
				
				amino = Main.chromDraw.getAminoChange(vardraw,base,vardraw.getExons().get(exon));
				
				if(amino.contains("UTR")) {
					
						if(VariantHandler.utr.isSelected()) {
							if(!VariantHandler.table.genearray.contains(vardraw.getExons().get(exon).getTranscript().getGene()) && vardraw.getExons().get(exon).getTranscript().getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
								if(!VariantHandler.onlyStats.isSelected()) {
									VariantHandler.table.addEntry(vardraw.getExons().get(exon).getTranscript().getGene());
								}
							}	
							vardraw.inVarList = true;						
							vardraw.getExons().get(exon).getTranscript().getGene().utr +=mutcount;
						
						}						
						else {
							continue;
						}
				}
				
				
				
				if(VariantHandler.nonsense.isSelected()) {
					if(!MethodLibrary.aminoEffect(amino).contains("nonsense")) {							
						continue;
					}						
				}
				if(VariantHandler.synonymous.isSelected()) {
					
					if(MethodLibrary.aminoEffect(amino).contains("synonymous") ) {	
						
						continue;
					}							
				}			
				
					for(int i = 0; i<entry.getValue().size(); i++) {
						if(entry.getValue().get(i).alleles != null) {
							break;
						}
						if(Main.drawCanvas.hideVar(entry.getValue().get(i), entry.getKey().length() > 1)) {
							continue;
						}
						if(VariantHandler.onlyselected.isSelected()) {
							if(!entry.getValue().get(i).getSample().equals(Main.drawCanvas.selectedSample)) {
								continue;
							}
						}
						if(!vardraw.getExons().get(exon).getTranscript().getGene().samples.contains(entry.getValue().get(i).getSample())) {
							if(!VariantHandler.onlyStats.isSelected()) {
								vardraw.getExons().get(exon).getTranscript().getGene().samples.add(entry.getValue().get(i).getSample());		
							}
						}
					}		
					if(!VariantHandler.onlyStats.isSelected()) {
						
						if(!VariantHandler.table.genearray.contains(vardraw.getExons().get(exon).getTranscript().getGene()) && vardraw.getExons().get(exon).getTranscript().getGene().samples.size() >= VariantHandler.geneSlider.getValue() ) {
							
							VariantHandler.table.addEntry(vardraw.getExons().get(exon).getTranscript().getGene());								
						}				
					}
					if(MethodLibrary.aminoEffect(amino).contains("nonsense")) {
						if(!calcgene.contains(vardraw.getExons().get(exon).getTranscript().getGene())) {
							for(int i = 0; i<entry.getValue().size(); i++) {
								if(entry.getValue().get(i).alleles != null) {
									break;
								}
								if(Main.drawCanvas.hideVar(entry.getValue().get(i), entry.getKey().length() > 1)) {
									continue;
								}
								entry.getValue().get(i).getSample().nonsense++;
								if(entry.getKey().length() > 1) {
									entry.getValue().get(i).getSample().fshift++;
								}
								if(amino.contains("spl")) {
									entry.getValue().get(i).getSample().splice++;
								}
								else {
									entry.getValue().get(i).getSample().nonsyn++;
								}
							}
							
							
						}
						vardraw.getExons().get(exon).getTranscript().getGene().nonsense +=mutcount;
						
					}
					else if(MethodLibrary.aminoEffect(amino).contains("missense")) {		
						if(!calcgene.contains(vardraw.getExons().get(exon).getTranscript().getGene())) {
							for(int i = 0; i<entry.getValue().size(); i++) {
								if(entry.getValue().get(i).alleles != null) {
									break;
								}
								if(Main.drawCanvas.hideVar(entry.getValue().get(i), entry.getKey().length() > 1)) {
									continue;
								}
								if(entry.getKey().length() > 1) {
									entry.getValue().get(i).getSample().inframe++;
								}
								entry.getValue().get(i).getSample().missense++;
								entry.getValue().get(i).getSample().nonsyn++;
							}						
						}
						vardraw.getExons().get(exon).getTranscript().getGene().missense +=mutcount;
					}
					else if(MethodLibrary.aminoEffect(amino).contains("synonymous")) {
						if(!calcgene.contains(vardraw.getExons().get(exon).getTranscript().getGene())) {
							for(int i = 0; i<entry.getValue().size(); i++) {
								if(entry.getValue().get(i).alleles != null) {
									break;
								}
								if(Main.drawCanvas.hideVar(entry.getValue().get(i), entry.getKey().length() > 1)) {
									continue;
								}
								entry.getValue().get(i).getSample().syn++;
							}						
						}
						vardraw.getExons().get(exon).getTranscript().getGene().synonymous +=mutcount;
					}					
				
					
				vardraw.inVarList = true;
				if(!vardraw.getExons().get(exon).getTranscript().getGene().varnodes.contains(vardraw)) {	
					if(!VariantHandler.onlyStats.isSelected()) {
						vardraw.getExons().get(exon).getTranscript().getGene().varnodes.add(vardraw);
					}
					preI = -1;
					
				}	
				if(v != preI) {
					vardraw.getExons().get(exon).getTranscript().getGene().mutations+=mutcount;	
					VariantHandler.table.variants +=mutcount;
					preI = v;
				}
				
				VariantHandler.table.revalidate();
				VariantHandler.table.repaint();
				if(!calcgene.contains(vardraw.getExons().get(exon).getTranscript().getGene())) {
					calcgene.add(vardraw.getExons().get(exon).getTranscript().getGene());
				}
			}	
		}	
		preI = v;
	
	
	if(!vardraw.isInGene() && VariantHandler.intergenic.isSelected()) {
		
	/*
			for(int v = 0; v<vardraw.vars.size(); v++) {
			 	entry = vardraw.vars.get(v);
				if(Main.drawCanvas.hideNodeVar(vardraw, entry)) {
					continue;
				}			
				
				base = entry.getKey();
				mutcount = 0;
				
				for(int m = 0; m<entry.getValue().size(); m++) {
						if(entry.getValue().get(m).alleles != null) {
							break;
						}
					
						if(Main.drawCanvas.hideVar(entry.getValue().get(m))) {
							continue;
						}
						mutcount++;
						
				}				
		*/	
		
			if(mutcount > 0) {
				
				if(VariantHandler.table.genearray.size() == 0  || !VariantHandler.table.genearray.get(VariantHandler.table.genearray.size()-1).intergenic || !VariantHandler.table.genearray.get(VariantHandler.table.genearray.size()-1).getName().equals(vardraw.getTranscripts().get(0).getGene().getName())) {
					Gene addgene =null;
					try {
						addgene = new Gene(vardraw.getTranscripts().get(0).getGene());
						
					}
					catch(Exception e) {
						e.printStackTrace();
						//Main.drawCanvas.ready("all");
						break;
					}
					addgene.intergenic = true;
					addgene.mutations = mutcount;
					VariantHandler.table.variants += mutcount;
					if(!VariantHandler.onlyStats.isSelected()) {
					addgene.varnodes.add(vardraw);
					VariantHandler.table.addEntry(addgene);
					}
					for(int m = 0; m<entry.getValue().size(); m++) {
						if(entry.getValue().get(m).alleles != null) {
							break;
						}
					
						if(Main.drawCanvas.hideVar(entry.getValue().get(m), entry.getKey().length() > 1)) {
							continue;
						}
						if(!addgene.samples.contains(entry.getValue().get(m).getSample())) {
							if(!VariantHandler.onlyStats.isSelected()) {
								addgene.samples.add(entry.getValue().get(m).getSample());
							}
						}
					}		
				}	
				else {
					if(!VariantHandler.onlyStats.isSelected()) {
						VariantHandler.table.genearray.get(VariantHandler.table.genearray.size()-1).mutations +=mutcount;					
						VariantHandler.table.genearray.get(VariantHandler.table.genearray.size()-1).varnodes.add(vardraw);
					}
					VariantHandler.table.variants += mutcount;
					for(int m = 0; m<entry.getValue().size(); m++) {
						if(entry.getValue().get(m).alleles != null) {
							break;
						}
					
						if(Main.drawCanvas.hideVar(entry.getValue().get(m), entry.getKey().length() > 1)) {
							continue;
						}
						if(!VariantHandler.table.genearray.get(VariantHandler.table.genearray.size()-1).samples.contains(entry.getValue().get(m).getSample())) {
							if(!VariantHandler.onlyStats.isSelected()) {
								VariantHandler.table.genearray.get(VariantHandler.table.genearray.size()-1).samples.add(entry.getValue().get(m).getSample());
							}
						}
					}		
		//		}
				
			
			}	
				vardraw.inVarList = true;
		}		
			
	}
		//	System.out.println(vardraw.getTranscripts().get(0).getGenename() +"..." +vardraw.getTranscripts().get(1).getGenename() );
		}
		if(!vardraw.inVarList) {			
			returnnode = vardraw.getNext();
			if(VariantHandler.allChroms.isSelected()) {
				vardraw.removeNode();
			}
			return returnnode;		
				
		}
		else {
			mutcount = 0;
			for(int v = 0; v<vardraw.vars.size(); v++) {
			 	entry = vardraw.vars.get(v);
			 	
				if(Main.drawCanvas.hideNodeVar(vardraw, entry)) {
					if(VariantHandler.allChroms.isSelected()) {
						vardraw.vars.remove(v);
						v--;
					}
					continue;
				}			
				
				base = entry.getKey();
				
			
				for(int m = 0; m<entry.getValue().size(); m++) {
						if(entry.getValue().get(m).alleles != null) {
							break;
						}
					
						if(Main.drawCanvas.hideVar(entry.getValue().get(m), entry.getKey().length() > 1)) {
							if(VariantHandler.allChroms.isSelected()) {
								entry.getValue().remove(m);
								m--;
							}
							continue;
						}
						
						mutcount++;
				}
			}
				
			if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {
				
				if(Main.drawCanvas.clusterNodes.size() == 0) {
					
					ClusterNode cluster = new ClusterNode();
					cluster.nodecount=mutcount;
					cluster.ID = vardraw.clusterId;
					vardraw.clusterNode = cluster;
					cluster.varnodes.add(vardraw);			
					cluster.width = 1;
					Main.drawCanvas.clusterNodes.add(cluster);
					
				}
				else if(Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).ID != vardraw.clusterId) {
					
					Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).width = Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.get(Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.size()-1).getPosition() -Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.get(0).getPosition() +1;
					ClusterNode cluster = new ClusterNode();
					cluster.nodecount+= mutcount;
					cluster.ID = vardraw.clusterId;
					vardraw.clusterNode = cluster;
					cluster.varnodes.add(vardraw);			
					cluster.width = 1;
					Main.drawCanvas.clusterNodes.add(cluster);
				}
				else {
					
					vardraw.clusterNode = Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1);
					Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).nodecount+=mutcount;
					Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.add(vardraw);
				}
				
				VariantHandler.clusterTable.variants+=mutcount;
				
			}
			if(!VariantHandler.writetofile.isSelected()) {
				if(vardraw != null && vardraw.inVarList && vardraw.bedhit && vardraw.getBedHits() != null) {		
						
						pretrack = -1;
						for(int i = 0; i<vardraw.getBedHits().size(); i++) {
							
								vardraw.getBedHits().get(i).inVarlist = true;
							
							if(pretrack != vardraw.getBedHits().get(i).getTrack().trackIndex) {
								vardraw.getBedHits().get(i).getTrack().getTable().variants += mutcount;
								pretrack = vardraw.getBedHits().get(i).getTrack().trackIndex;
							}
							if(vardraw.getBedHits().get(i).getTrack().getTable().bedarray == null) {
								vardraw.getBedHits().get(i).getTrack().getTable().bedarray = new ArrayList<BedNode>();
							}
							if(!vardraw.getBedHits().get(i).getTrack().getTable().bedarray.contains(vardraw.getBedHits().get(i))) {
								
								vardraw.getBedHits().get(i).mutations+=mutcount;
								
									vardraw.getBedHits().get(i).getTrack().getTable().bedarray.add(vardraw.getBedHits().get(i));
									vardraw.getBedHits().get(i).getTrack().getTable().setPreferredSize(new Dimension(vardraw.getBedHits().get(i).getTrack().getTable().tablescroll.getViewport().getWidth(), (vardraw.getBedHits().get(i).getTrack().getTable().getTableSize()+2+samplecount)*vardraw.getBedHits().get(i).getTrack().getTable().rowHeight));										
									vardraw.getBedHits().get(i).getTrack().getTable().revalidate();		
								
							//	vardraw.inVarList = true;
							}
							else {
								vardraw.getBedHits().get(i).mutations+=mutcount;
							//	vardraw.inVarList = true;
								
							}
						}
				}
			}
		}
		if(VariantHandler.writetofile.isSelected()) {
			
			writeToFile(vardraw, output, outputgz);
		}
		lastVar = vardraw;
		
		return vardraw.getNext();
//	}
}

/*
static void fetchEnsembl2() {
	try {
	MysqlDataSource dataSource = new MysqlDataSource();
	dataSource.setUser("genome");
	dataSource.setPassword("");
	dataSource.setServerName("genome-mysql.cse.ucsc.edu");
	
	dataSource.setDatabaseName("hg19");
	
	Connection conn = dataSource.getConnection();
	
	Statement stmt = conn.createStatement();
	
	

	
//	ResultSet rs = stmt.executeQuery("select S.*,X.*,G.* from ensemblSource as S,ensGene as G,knownToEnsembl as KE, kgXref as X where S.name = G.name and G.name=KE.value and KE.name=X.kgID");
	ResultSet rs = stmt.executeQuery("select E.chrom, E.txStart, E.txEnd, N.value, E.exonCount, E.strand, E.name2, E.name, X.spID, C.transcript, S.source, E.cdsStart, E.cdsEnd, E.exonStarts, E.exonEnds, E.exonFrames, X.description from "
			+ "knownToEnsembl as K left outer join knownCanonical as C on K.name = C.transcript "
			+ "left outer join ensGene as E on K.value = E.name "
			+ "left outer join kgXref as X on K.name = X.kgID "
			+ "left outer join ensemblSource as S on S.name = E.name "
			+ "left outer join ensemblToGeneName as N on E.name = N.name");
	
	

	BufferedWriter writer = new BufferedWriter(new FileWriter("exons.txt"));
	String header = "#Chrom\tGeneStart\tGeneEnd\tName\tExonCount\tStrand\tENSG\tENST\tUniProt\tCanonical\tBiotype\tCodingStart\tCodingEnd\tExonStarts\tExonEnds\tStartPhases\tDescription";
	writer.write(header +"\n");
	
	while(rs.next()) {
		
		
	
		writer.write(rs.getString(1).substring(3) +"\t");
	
		for(int i = 2; i <= rs.getMetaData().getColumnCount(); i++) {
				if(rs.getString(i) == null || rs.getString(i).equals("")) {
			
					writer.write("-\t");
				}
				else if(i == 10) {
					if(rs.getString(i).equals("null")) {
					
						writer.write("0\t");
					}
					else {
				
						writer.write("1\t");
					}
				}
				else {
			
					writer.write(rs.getString(i) +"\t");
				}
				
		}
	
		writer.write("\n");
		
	}
	writer.close();
	rs.close();
	
	
	
}
catch(Exception e) {
	ErrorLog.addError(e.getStackTrace());
	e.printStackTrace();
}
}
static void fetchEnsembl() {
	try {
	MysqlDataSource dataSource = new MysqlDataSource();
	dataSource.setUser("anonymous");
	dataSource.setPassword("");
	dataSource.setServerName("ensembldb.ensembl.org");
	dataSource.setPort(3337);
	dataSource.setDatabaseName("homo_sapiens_core_84_37");
	
	Connection conn = dataSource.getConnection();
	
	Statement stmt = conn.createStatement();
	
	

	
//	ResultSet rs = stmt.executeQuery("select S.*,X.*,G.* from ensemblSource as S,ensGene as G,knownToEnsembl as KE, kgXref as X where S.name = G.name and G.name=KE.value and KE.name=X.kgID");
	ResultSet rs = stmt.executeQuery("SELECT transcript.stable_id, seq_region.name, transcript.seq_region_strand, transcript.biotype, start_exon.stable_id AS start_exon_id, translation.seq_start, end_exon.stable_id AS end_exon_id ,translation.seq_end " +
            "FROM translation JOIN transcript ON translation.transcript_id=transcript.transcript_id " +
            "JOIN exon AS start_exon ON translation.start_exon_id = start_exon.exon_id "+
            "JOIN exon AS end_exon ON translation.end_exon_id = end_exon.exon_id " +
            "JOIN gene ON transcript.gene_id = gene.gene_id " +
            "JOIN seq_region ON gene.seq_region_id = seq_region.seq_region_id " +
            "WHERE seq_region.coord_system_id = 2");
	
	
	

	
	BufferedWriter writer = new BufferedWriter(new FileWriter("exons.txt"));
	String header = "#Chrom\tGeneStart\tGeneEnd\tName\tExonCount\tStrand\tENSG\tENST\tUniProt\tCanonical\tBiotype\tCodingStart\tCodingEnd\tExonStarts\tExonEnds\tStartPhases\tDescription";
	writer.write(header +"\n");
	
	while(rs.next()) {
		
		
	
		writer.write(rs.getString(1).substring(3) +"\t");
		System.out.print(rs.getString(1).substring(3) +"\t");
		for(int i = 2; i <= rs.getMetaData().getColumnCount(); i++) {
				if(rs.getString(i) == null || rs.getString(i).equals("")) {
					
					writer.write("-\t");
					System.out.print("-\t");
				}
				else if(i == 10) {
					if(rs.getString(i).equals("null")) {
				
						writer.write("0\t");
						System.out.print("0\t");
					}
					else {
					
						writer.write("1\t");
						System.out.print("1\t");
					}
				}
				else {
				
					writer.write(rs.getString(i) +"\t");
					System.out.print(rs.getString(i) +"\t");
				}
				
		}
		System.out.println();
		writer.write("\n");
		
	}
	writer.close();
	rs.close();
	
	
	
}
catch(Exception e) {
	ErrorLog.addError(e.getStackTrace());
	e.printStackTrace();
}
}
static void fetchAnnotation() {
	try {
		
		MysqlDataSource dataSource = new MysqlDataSource();
		dataSource.setUser("genome");
		dataSource.setPassword("");
		dataSource.setServerName("genome-mysql.cse.ucsc.edu");
		
	//	dataSource.setDatabaseName("hg19");
		dataSource.setDatabaseName("mm10");
		Connection conn = dataSource.getConnection();
		
		Statement stmt = conn.createStatement();
				
	//	ResultSet rs = conn.getMetaData().getCatalogs();
	
	//	ResultSet rs = stmt.executeQuery("select S.*,X.*,G.* from ensemblSource as S,ensGene as G,knownToEnsembl as KE, kgXref as X where S.name = G.name and G.name=KE.value and KE.name=X.kgID");
		System.out.println("Fetching...");
		ResultSet rs = stmt.executeQuery("select E.chrom, E.txStart, E.txEnd, N.value, E.exonCount, E.strand, E.name2, E.name, E.name, E.name, S.source, E.cdsStart, E.cdsEnd, E.exonStarts, E.exonEnds, E.exonFrames, E.name from "
				//+ "knownToEnsembl as K left outer join knownCanonical as C on K.name = C.transcript "
				+ "ensGene as E left outer join ensemblToGeneName as N on E.name = N.name "
				+ "left outer join ensemblSource as S on S.name = E.name");				
				
		System.out.println("Fecthed");
	
		BufferedWriter writer = new BufferedWriter(new FileWriter("exons.txt"));
		String header = "#Chrom\tGeneStart\tGeneEnd\tName\tExonCount\tStrand\tENSG\tENST\tUniProt\tCanonical\tBiotype\tCodingStart\tCodingEnd\tExonStarts\tExonEnds\tStartPhases\tDescription";
		writer.write(header +"\n");
		
		while(rs.next()) {
			
			//ResultSetMetaData joo = rs.getMetaData();
			writer.write(rs.getString(1).substring(3) +"\t");
		//	System.out.print(rs.getString(1).substring(3) +"\t");
			for(int i = 2; i <= rs.getMetaData().getColumnCount(); i++) {
					if(rs.getString(i) == null || rs.getString(i).equals("")) {
				//		System.out.print("-\t");
						writer.write("-\t");
					}
					else if(i == 10) {
						if(rs.getString(i).equals("null")) {
						//	System.out.print("0\t");
							writer.write("0\t");
						}
						else {
						//	System.out.print("1\t");
							writer.write("1\t");
						}
					}
					else {
					//	System.out.print(rs.getString(i) +"\t");
						writer.write(rs.getString(i) +"\t");
					}
					
			}
		//	System.out.println();
			writer.write("\n");
			
		}
		writer.close();
		rs.close();
		
		
		
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
}*/
public int[][] variantCaller(String chrom, int startpos, int endpos, Reads readClass, int minBaseQuality, int minReadQuality, ReferenceSeq reference) {
	
		int[][] coverages = new int[endpos-startpos +300][8];		
		Iterator<SAMRecord> bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos, endpos);	
		
		try {
			readClass.setCoverageStart(startpos);
			readClass.startpos = startpos;
			readClass.endpos = endpos;
			readClass.setReadStart(startpos);
			readClass.setReadEnd(endpos);			
					
			readClass.setCoverageEnd(startpos + coverages.length);		
		
			
			while(bamIterator != null && bamIterator.hasNext()) {	
				if(!Main.drawCanvas.loading) {
					Main.drawCanvas.ready("all");
					break;
				}
				try {	
					
				try {
					samRecord = bamIterator.next(); 	
					
				}
				catch(htsjdk.samtools.SAMFormatException ex) {
					ex.printStackTrace();					
				}
				if(samRecord.getMappingQuality() < minReadQuality) {
					continue;
				}
				if(samRecord.getReadUnmappedFlag()) {					
					continue;
				}
				
				//TODO jos on pienempi ku vika unclipped start		
				
				
				if(samRecord.getUnclippedEnd() < startpos) { //this.readSeqStart+1) {
					
					continue;
				}
				
				if(samRecord.getUnclippedStart() >= endpos){
					
					
					break;
				}		
				
				if(readClass.sample.longestRead <samRecord.getCigar().getReferenceLength()) {
					readClass.sample.longestRead = samRecord.getCigar().getReferenceLength();
					
				}
				if(readClass.sample.getComplete() == null) {
					if(samRecord.getReadName().startsWith("GS")) {
						readClass.sample.setcomplete(true);
					}
					else {
						readClass.sample.setcomplete(false);
					}
				}
				
			//	java.util.ArrayList<java.util.Map.Entry<Integer,Integer>> insertions = null;
			
				if(samRecord.getCigarLength() > 1) {			
					readstart = samRecord.getUnclippedStart();			
					readpos= 0;
					mispos= 0;
					/*if(readClass.getCoverageStart()>readstart) {
						return null;
					}*/
					
					for(int i = 0; i<samRecord.getCigarLength(); i++) {
						element = samRecord.getCigar().getCigarElement(i);
						if(element.getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH)== 0) {
							for(int r = readpos; r<readpos+element.getLength(); r++) {																	
								if(((readstart+r)-readClass.getCoverageStart()) < 0) {
									continue;
								}
								if(samRecord.getReadBases()[mispos] !=  reference.getSeq()[((readstart+r)-reference.getStartPos()-1)]) {	
									if((int)samRecord.getBaseQualityString().charAt(mispos)-readClass.sample.phred >= minBaseQuality) {
										
										coverages[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get((byte)samRecord.getReadBases()[mispos])]++;
									}									
									
								}
							
								coverages[((readstart+r)-readClass.getCoverageStart())][0]++; 
								
								if(coverages[((readstart+r)-readClass.getCoverageStart())][0] > readClass.getMaxcoverage()) {
									readClass.setMaxcoverage(coverages[((readstart+r)-readClass.getCoverageStart())][0]);
								}
								mispos++;
							}
							
							readpos+=element.getLength();
						}
						else if(element.getOperator().compareTo(CigarOperator.DELETION)== 0) {
					//		readWidth+=element.getLength();						
						/*	for(int d = 0; d<element.getLength(); d++) {
								coverages[(readstart+readpos+d)- readClass.getCoverageStart()][Main.baseMap.get((byte)'D')]++;
								
							}*/
							readpos+=element.getLength();
						}
						else if(element.getOperator().compareTo(CigarOperator.INSERTION)== 0) {										
						//	coverages[(readstart+readpos)- readClass.getCoverageStart()][Main.baseMap.get((byte)'I')]++;
							mispos+=element.getLength();					
						
						}
						else if(element.getOperator().compareTo(CigarOperator.SOFT_CLIP)== 0 ) {
							/*	if(Settings.softClips == 1) {
									for(int r = readpos; r<readpos+element.getLength(); r++) {
																				
										if(samRecord.getReadBases()[mispos] != reference.getSeq()[((readstart+r)-reference.getStartPos()-1)]) {													
											coverages[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get(samRecord.getReadBases()[mispos])]++;
																				
										}
									
										coverages[((readstart+r)-readClass.getCoverageStart())][0]++; 
										
										mispos++;
									}
									readpos+=element.getLength();
								}
								else {*/
									if(i == 0) {							
										readstart = samRecord.getAlignmentStart();							
										mispos+=element.getLength();
									}
								//}							
						}
						else if (element.getOperator().compareTo(CigarOperator.HARD_CLIP)== 0) {
							if(i == 0) {
								readstart = samRecord.getAlignmentStart();
							}
						}
						else if(element.getOperator().compareTo(CigarOperator.SKIPPED_REGION)== 0) {										
							readpos+=element.getLength();				
						}
					}						
					
				}
				else {
					readstart = samRecord.getUnclippedStart();
					
					for(int r = 0; r<samRecord.getReadLength(); r++) {
						try {									
							
							if(samRecord.getReadBases()[r] != reference.getSeq()[((readstart+r)-reference.getStartPos()-1)]) {									
								
								if((readstart+r)-readClass.getCoverageStart() < coverages.length-1 && (readstart+r)- readClass.getCoverageStart() > -1) {
									
									if((int)samRecord.getBaseQualityString().charAt(r)-readClass.sample.phred >= minBaseQuality) {
										
										coverages[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get(samRecord.getReadBases()[r])]++;	
									}
								}									
							}										
					//		try {
								if((readstart+r)-readClass.getCoverageStart() < 0) {
									continue;
								}
								coverages[(readstart+r)-readClass.getCoverageStart()][0]++; 
						/*	}
							catch(Exception e) {							
								e.printStackTrace();
							}
					/*		if(coverages[(samRecord.getUnclippedStart()-readClass.getCoverageStart())+r][0] > readClass.getMaxcoverage()) {
								readClass.setMaxcoverage(coverages[(samRecord.getUnclippedStart()-readClass.getCoverageStart())+r][0]);
							}*/
						}
						catch(Exception e) {						
							ErrorLog.addError(e.getStackTrace());
							e.printStackTrace();							
							break;
						}
					}							
				
				}				
				}
				catch(Exception e) {				
					e.printStackTrace();
					break;
				}				
			}		
		}
		catch(Exception ex) {			
			ex.printStackTrace();	
		
			return null;	
		}		
		
		return coverages;
	}
	public static void main(String args[]) {
		/*String cram = "X:/cg7/Heikki/CRAMtest/cram3/Y_crc47_1_13-0818_cram3_test.cram";
		String crai = cram + ".crai";
		File ref = new File("C:/HY-Data/RKATAINE/BasePlayer/BasePlayer/genomes/Homo_sapiens_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa");
		*/
		
		try {
		/*	String file = "X:/cg8/Riku/PlatinumGenomes/Ohutsuoli9_1_15-1064TP_T_160129_K00110_0062_AH55TWBBXX_1.snps_indels.ug.vcf.idx";
			BufferedReader reader = new BufferedReader(new FileReader(file.replace(".idx", "")));
			
		//	VCFFileReader reader = new VCFFileReader(new File(file));
			Index indexfile = IndexFactory.loadIndex(file);
			List<Block> list = indexfile.getBlocks("1", 100000000,100000000);
			
			System.out.println(list.get(0).getStartPosition());
			
			reader.skip(list.get(0).getStartPosition());
			int count = 0;
			while(count < 20) {
				System.out.println(reader.readLine());
				
				count++;
			}
			reader.close();
			
			/*
			
			while(iter.hasNext()) {
				
				System.out.println(iter.next().getStart());
				
			}
			read.close();
			*/
	/*		RandomAccessFile random = new RandomAccessFile(ref, "r");
		CRAMFileReader reader = new CRAMFileReader(new File(cram), new File(crai),new ReferenceSource(ref), random, ValidationStringency.SILENT);
		QueryInterval[] interval = { new QueryInterval(0, 1000000, 1010000) };
		
		Iterator<SAMRecord> iter = reader.query(interval, true);
		System.out.println(iter.hasNext());
	//	while(iter.hasNext()) {
	//		SAMRecord rec = iter.next();
		//	System.out.println(rec.getAlignmentStart() +" " +rec.getStringAttribute("MD"));
	//	}
		/*SAMRecord record = iter.next();
		do {
			System.out.println(record.getAlignmentStart());
			record= iter.next();
		
		}
		while(record != null);
		*/
		
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
	}
}
