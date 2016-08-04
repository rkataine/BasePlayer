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
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import htsjdk.tribble.readers.TabixReader;

import htsjdk.tribble.index.Block;

import htsjdk.tribble.index.tabix.TabixIndex;

import htsjdk.samtools.CRAMFileReader;
import htsjdk.samtools.cram.ref.ReferenceSource;
//import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;

import java.awt.Dimension;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.InputStreamReader;
import java.net.URL;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.SwingWorker;

import com.mysql.jdbc.jdbc2.optional.MysqlDataSource;

public class FileRead extends SwingWorker< String, Object >  {
	
	Boolean readVCF = false, changeChrom = false, readBAM = false,searchInsSites = false, varcalc = false, getreads = false;
	File[] files;
	String chrom = "0";
	int pos;
	Short calls1, calls2;	
	
	private boolean found, noref = false;
	private int xpos;	
	private ReadNode mundane = null;
	static int searchwindow = 1000;
	private boolean left;
	boolean stop = false;
	private ReadNode lastAdded;	
	private ReadNode addNode;
	TabixReader tabixreader;
	//VCFFileReader vcfreader;
	SamReader samFileReader;	
	CRAMFileReader CRAMReader = null;
	Iterator<SAMRecord> bamIterator = null;	
	SAMRecord samRecord, samRecordBuffer;
	private String[] info;
	private String[] coverages;
	private Sample sample;
	private short quality;
	static boolean novars = false;
	int start, end, viewLength;
	double pixel;
//	private TabixReader.Iterator iterator;
//	private Iterator<VariantContext> vcfIterator;
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
	private int foundexon;
	public boolean statcalc;
	private boolean genotype;
	private int mutcount, linecounter = 0;
	private short firstallele;
	private short secondallele;
//	ArrayList<SAMRecord> samBuffer = new ArrayList<SAMRecord>();	
	private String altbase;
	private short refallele;
	private Short refcalls;
	private short altallele;
	private Short altcalls;
	private String altbase2;
	private String[] headersplit;
	private boolean first;
	private int samplecount;
//	private int index;
//	private int[] tempadder;
	private int middle;
//	private ReadBuffer readHead;
	private ReadBuffer currentread;
//	private Integer hashindex;
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
	private int dpindex;
	private int previewLength = Integer.MAX_VALUE;
	static boolean bigcalc=false;
	static boolean changing = false;
//	static int readcount = 0;	
	static boolean readFiles;	
	public static boolean search;
	public static int searchStart;
	public static int searchEnd;
	public static boolean cancelvarcount= false;
	public static boolean cancelfileread = false;
	public static boolean cancelreadread= false;
	public static BufferedWriter output = null;
	public static String outputName = "";
	
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
			if(FileRead.head.getNext() == null) {
				if(VariantHandler.allChroms.isSelected() && Main.selectedChrom != 0) {
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
void setVCFFileStart(String chrom, int start, int end, Sample sample) {
	try {
		
			
		
			TabixIndex index = new TabixIndex(new File(sample.getTabixFile()+".tbi"));
			java.util.List<Block> blocks = null;
			
			if(index.containsChromosome(sample.vcfchr +chrom)) {
				
				try {
					
					blocks = index.getBlocks(sample.vcfchr +chrom, start, end);
					
				}
				catch(Exception e) {
					sample.vcfEndPos = 0;
					
					return;
				}
			
				if(blocks.size() > 0) {
					sample.getVCFInput().seek(0);
					sample.getVCFInput().seek(blocks.get(0).getStartPosition());
					sample.vcfEndPos = blocks.get(blocks.size()-1).getEndPosition();
				
					
				}
				else {
					sample.getVCFInput().seek(0);
				}
				
			}
			else {
				if(index.containsChromosome(sample.vcfchr +(Main.chromosomeDropdown.getSelectedIndex()+1))) {
					try {
						blocks = index.getBlocks(sample.vcfchr+(Main.chromosomeDropdown.getSelectedIndex()+1), start, end);
					}
					catch(Exception e) {
						sample.vcfEndPos = 0;
						return;
					}
					sample.getVCFInput().seek(0);
			
					if(blocks.size() > 0) {
						sample.getVCFInput().seek(blocks.get(0).getStartPosition());
					}
					else {
						sample.getVCFInput().seek(0);
					}
				}
				else {
					sample.vcfEndPos = 0;
				}
			}						
		
	}
	catch(Exception e) {
		e.printStackTrace();

		
	}
	
}
void cancelFileRead() {
	changing = false;
	FileRead.search = false;
	bigcalc = false;
	Main.opensamples.setText("Open samples");
	head.putNext(null);
	current = null;
	Main.drawCanvas.current = FileRead.head;
	Draw.updatevars = true;
	Main.drawCanvas.repaint();
	tabixreader = null;		
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
	removeNonListVariants();
	removeBedLinks();
	for (int i = 0; i<Control.controlData.fileArray.size(); i++) {
	  Control.controlData.fileArray.get(i).controlled = false;
	}
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
	
	if(Main.bedCanvas.bedOn) {
		for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
			if(Main.bedCanvas.bedTrack.get(i).intersect && !Main.bedCanvas.bedTrack.get(i).loading) {
				Main.bedCanvas.annotate(Main.bedCanvas.bedTrack.get(i).getHead());
			}		
			else if(Main.bedCanvas.bedTrack.get(i).intersect && Main.bedCanvas.bedTrack.get(i).loading) {
				Main.bedCanvas.bedTrack.get(i).waiting = true;
			}
		}
	}

	if(Control.controlData.controlsOn) {		    	
		Control.applyControl();
	}
	if(bigcalc) {
		Main.drawCanvas.calcClusters(FileRead.head);
	}
	else {
		Main.drawCanvas.calcClusters(FileRead.head,1);
	}
	
}

void getVariants(String chrom, int start, int end, Sample sample) {
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
	
		setVCFFileStart(chrom, start, end, sample);

	sample.maxCoverage = 0;
	
	current = head;
	line = "";
	first = true;
	stop = false;
		while(true) {	    
			if(stop) {
				
				break;
			}
			if(cancelfileread || cancelvarcount) {
				cancelFileRead();
				
				break;
	  		}
			
				try {
					
					line = sample.getVCFInput().readLine();
					
					if(line == null) {
						break;
					}
				}
				catch(Exception ex) {							
					JOptionPane.showMessageDialog(Main.chromDraw, ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
					ErrorLog.addError(ex.getStackTrace());
					ex.printStackTrace();
					Main.cancel();
					changing = false;
					break;
				}												
			
				split = line.split("\\t+");
				if(sample.getVCFInput().getFilePointer() > sample.vcfEndPos ) {			
					break;
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
			
			
		}
		
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
		JOptionPane.showMessageDialog(Main.chromDraw, exc.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
		System.out.println(sample.getName());
		exc.printStackTrace();
		ErrorLog.addError(exc.getStackTrace());
		changing = false;
		
	}
	
}
	
	public void changeChrom(String chrom) {
		try {
			
			
			cancelfileread = false;
			
		try {	
			Main.drawCanvas.loading("Loading annotation...");
		//	if(Main.drawCanvas.splits.get(0).getTranscripts().size() == 0 || !Main.drawCanvas.splits.get(0).getTranscripts().get(0).getChrom().equals(chrom)) {			
			Main.drawCanvas.splits.get(0).setGenes(getExons(chrom));						
			Main.chromDraw.updateExons = true;
			Main.chromDraw.repaint();
		//	}
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
					
					if(Main.bedCanvas.bedTrack.get(i).small) {
						
						if(Main.bedCanvas.bedTrack.get(i).getHead().getNext() == null) {
							
							Main.bedCanvas.getBEDfeatures(Main.bedCanvas.bedTrack.get(i), 1, Main.drawCanvas.splits.get(0).chromEnd);
							
						}
					}
					else {
						if(searchEnd- searchStart < 1000000) {
							Main.bedCanvas.getBEDfeatures(Main.bedCanvas.bedTrack.get(i), searchStart, searchEnd);
						}
					}
				}				
			}
			else {
				
				for(int i = 0; i< Main.bedCanvas.bedTrack.size(); i++) {
					if(Main.bedCanvas.bedTrack.get(i).small) {						
						Main.bedCanvas.getBEDfeatures(Main.bedCanvas.bedTrack.get(i), 1, Main.drawCanvas.splits.get(0).chromEnd);
					}
				}
			}
			Main.drawCanvas.ready("Loading BED-files...");
			
		}
		if(novars) {
			Main.drawCanvas.variantsStart= 0;
			Main.drawCanvas.variantsEnd = 0;
			
		}
		else {
			changing = true;
		}
		 
		if(Main.varsamples > 0 && !novars) {			
			removeNonListVariants();
			removeBedLinks();
			
			Main.drawCanvas.loading("Loading variants...");
			head.putNext(null);
			current = FileRead.head;			
			Main.drawCanvas.current = head;		
			Main.chromDraw.varnode = null;
			Main.chromDraw.vardraw = null;
			for(int i = 0; i<Main.samples; i++) {
				if(cancelfileread) {
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
					getVariants(chrom, 1, Main.drawCanvas.splits.get(0).chromEnd, Main.drawCanvas.sampleList.get(i));
				}
				
			}
			
			annotate();
			
			if(Main.bedCanvas.bedOn) {
				for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
					if(Main.bedCanvas.bedTrack.get(i).intersect && !Main.bedCanvas.bedTrack.get(i).loading) {
						Main.bedCanvas.annotate(Main.bedCanvas.bedTrack.get(i).getHead());
					}		
					else if(Main.bedCanvas.bedTrack.get(i).intersect && Main.bedCanvas.bedTrack.get(i).loading) {
						Main.bedCanvas.bedTrack.get(i).waiting = true;
					}
				}
			}
	
	//		Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.samples);
	//		Main.drawCanvas.checkSampleZoom();
		  	Draw.updatevars = true;
		  	
		  	Main.drawCanvas.ready("Loading variants...");
		  	if(!Main.drawCanvas.loading) {
				Draw.calculateVars = true;
			}
		  	Main.drawCanvas.loading("Applying controls...");
			if(Control.controlData.controlsOn) {		    	
			    Control.applyControl();
			}
			Main.drawCanvas.ready("Applying controls...");
		}
		if(novars) {
			Main.drawCanvas.variantsStart= 1;
			Main.drawCanvas.variantsEnd = Main.drawCanvas.splits.get(0).chromEnd;
		}
		FileRead.novars = false;
		FileRead.search = false;
		readFiles = false;
		
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
		
	}
	
	public ArrayList<Gene> getExons(String chrom) {		
		
		ArrayList<Gene> transcriptsTemp = new ArrayList<Gene>();
		
		try {			
			if(Main.genomehash.get(Main.defaultGenome).size() == 0) {
				return new ArrayList<Gene>();
			}
			
			ChromDraw.exonReader = new TabixReader(Main.genomehash.get(Main.defaultGenome).get(Main.annotation).getCanonicalPath());
			if(chrom == null) {
				return null;
			}
			TabixReader.Iterator exonIterator = null;
				try {				
					
					exonIterator = ChromDraw.exonReader.query(chrom);
					
				}
				catch(Exception e) {
					try {
						
						if(chrom.matches("\\w+")) {
							if(chrom.contains("X")) {
								
								exonIterator = ChromDraw.exonReader.query("23");
							}
							else if(chrom.contains("Y")) {
								exonIterator = ChromDraw.exonReader.query("24");
							}
							else if(chrom.contains("M")) {
								exonIterator = ChromDraw.exonReader.query("25");
							}
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
					addtrans.setGene(addgene);
					transcriptsTemp.add(addgene);														
				}
				else {
					setGene = genes.get(exonSplit[6]);
					setGene.addTranscript(addtrans);
					addtrans.setGene(setGene);					
				}					
			}
			genes.clear();
			ChromDraw.exonReader.close();
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
		}
		if(transcriptsTemp != null && transcriptsTemp.size() > 0 && transcriptsTemp.get(0).getCanonical() == null) {
			VariantHandler.allIsoforms.setEnabled(false);
		}
		else {
			VariantHandler.allIsoforms.setEnabled(true);
		}
		return transcriptsTemp;
	}
	private void readBAM(File[] files) {
		try {
		 
		 File addFile =null;
	  	 File[] addDir;
	  	 Sample currentSample = null;
	  	 
	
	  	 for(int i = 0; i<files.length; i++) { 	       
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
		      	  currentSample = new Sample(addFile.getName(), (Short)Main.samples, null);	      	 
		      	  Main.drawCanvas.sampleList.add(currentSample);	      	 
		      	  currentSample.samFile = addFile;
		      	  currentSample.resetreadHash();
		      	 
		      
		      try {
		    	  if(currentSample.samFile.getName().endsWith(".cram")) {
		    		 CRAMReader = new CRAMFileReader(currentSample.samFile, new File(currentSample.samFile.getCanonicalFile() +".crai"),new ReferenceSource(Main.ref), ValidationStringency.SILENT);
		    		 if(CRAMReader.getFileHeader().getTextHeader().contains("SN:chr")) {
		    			 currentSample.chr = "chr";
		    		 }
		    		
		    			CRAMReader.close();
		    	  }
		    	  else {
		    		  samFileReader = SamReaderFactory.makeDefault().open(currentSample.samFile);
			      	  if(samFileReader.getFileHeader().getTextHeader().contains("SN:chr")) {
			      		  currentSample.chr = "chr";
			      	  }
			      	 samFileReader.close();
		    	  }
		      	 }
		      	 catch(Exception e) {
		      		ErrorLog.addError(e.getStackTrace());
		      		 e.printStackTrace();
		      	 }
				  Main.samples++;
		         }   
	  		 
	  	 }
	   	 Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.samples);
	  	 Main.drawCanvas.checkSampleZoom();
	  	 
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
		 Main.drawCanvas.splits.get(0).updateReads = true;
		 Draw.updatevars = true;
		 Main.drawCanvas.repaint();
}	  	 
	
private void readVCF(File[] files) {
		try {	 
	  
  	  String line;
  	  String[] split;  	  
  	  Main.drawCanvas.loading("Loading variants...");
  	  File[] addDir;
  	  int sampletemp = Main.samples;
  	
	   	Sample addSample = null;
	  	StringBuffer sampleString = new StringBuffer(); 
	  	int fileindex = -1;
	  	readFiles = true;
	  	cancelfileread = false;
	  	
	 if(Control.controlData.controlsOn) {	    	
	    Control.dismissControls(head);
	 }
  	
  	 for(int fi = 0; fi<files.length; fi++) { 
  	    
  		if(Main.cancel) {
  			current = null;
  			FileRead.head.putNext(null);
  			return;
  		}
  	    addDir = null;  	   
  	   
  	    if(files[fi].isDirectory()) {  	    	
  	    	
  	    	addDir = files[fi].listFiles(new FilenameFilter() {
  	    	     public boolean accept(File dir, String name) {
  	    	        return name.toLowerCase().endsWith(".vcf.gz");
  	    	     }
  	    	});
  	    	
  	    	for(int f= 0; f<addDir.length; f++) {
  	    		if(cancelfileread) {
  	    			current = null;
  	    			FileRead.head.putNext(null);
  	    			return;
  	    		}
  	    		 addSample = new Sample(addDir[f].getName(), (Short)Main.samples, addDir[f].getCanonicalPath());      
  	    		 Main.drawCanvas.sampleList.add(addSample);
  	        	 Main.varsamples++;
  	        	 Main.samples++;
  	        	 VariantHandler.commonSlider.setMaximum(Main.varsamples);
  	        	 VariantHandler.commonSlider.setUpperValue(Main.varsamples);
  	        	 VariantHandler.geneSlider.setMaximum(Main.varsamples);
  	        	 sampleString.append(addSample.getName() +";");
  	    		
  	    		}
  	    		File[] bamfiles = files[fi].listFiles(new FilenameFilter() {
  	    		 
  	    	     public boolean accept(File dir, String name) {
  	    	        return name.toLowerCase().endsWith(".bam");
  	    	     }
  	    	   	 });
  	    	  	 if(bamfiles.length > 0) {
  	    	  		
  	    	  		int index = -1, sampleindex;
  	    	  		 for(int i = 0; i<bamfiles.length; i++) {
  	    		  		 sampleindex = 0;
  	    		  		 index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".bam")));
  	    		  	 	 if (index < 0) continue; 
  	    		  	 	 
  	    		  	 	 for(char letter : sampleString.substring(0, index).toString().toCharArray()) {
  	    		  	 		 if(letter == ';') sampleindex++;
  	    		  	 	 }
  	    		  	 	 Main.drawCanvas.sampleList.get(sampleindex+sampletemp).resetreadHash();
  	    		  	 	 Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile = new File(bamfiles[i].getAbsolutePath());	 
 	    			  	  samFileReader = SamReaderFactory.makeDefault().open(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
  	    		      	  if(samFileReader.getFileHeader().getTextHeader().contains("SN:chr")) {
  	    		      		Main.drawCanvas.sampleList.get(sampleindex+sampletemp).chr = "chr";
  	    		      	  }  	    		      	 
  	    		  	 }  	
  	    	  	 }
  	    	  	 else {
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
  	    	  	 }
  	    	  	
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
  	    	fileindex =fi;
  	    	
        	addSample = new Sample(files[fi].getName(), (Short)Main.samples, files[fi].getCanonicalPath());        	
        	
        	 Main.drawCanvas.sampleList.add(addSample);
          	 Main.varsamples++;
          	 Main.samples++;
          	 VariantHandler.commonSlider.setMaximum(Main.varsamples);
          	 VariantHandler.commonSlider.setUpperValue(Main.varsamples);
          	 VariantHandler.geneSlider.setMaximum(Main.varsamples);
          	 sampleString.append(addSample.getName() +";");
       	  }     	
  	 }
  	
  	 if(fileindex > -1) {
     File[] bamfiles = files[fileindex].getParentFile().listFiles(new FilenameFilter() {
     public boolean accept(File dir, String name) {
        return name.toLowerCase().endsWith(".bam");
     }
   	 });
  	 if(bamfiles.length > 0) {
  		
  		int index = -1, sampleindex;
  		 for(int i = 0; i<bamfiles.length; i++) {
	  		 sampleindex = 0;
	  		 index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".")));
	  	 	 if (index < 0) continue; 
	  	 	 
	  	 	 for(char letter : sampleString.substring(0, index).toString().toCharArray()) {
	  	 		 if(letter == ';') sampleindex++;
	  	 	 }
	  	 	
	  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile = new File(bamfiles[i].getAbsolutePath());	 
	  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).resetreadHash();
	  	 	samFileReader = SamReaderFactory.makeDefault().open(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
      	    
	  	    if(samFileReader.getFileHeader().getTextHeader().contains("SN:chr")) {
      	    	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).chr = "chr";
      	    }
	      	samFileReader.close();
	  	 }  	
  	 }
  	
  		 
	  	 File[] listfiles = files[0].getParentFile().listFiles(new FilenameFilter() {
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
  	 
  	 	sampleString = new StringBuffer();
  	 }
  	
  	 files = null;
  
  //	Main.drawCanvas.drawVariables.visibleend = (short)(Main.drawCanvas.sampleList.size()-1);
  	Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.drawCanvas.sampleList.size());
  	Main.drawCanvas.checkSampleZoom();
  	Main.drawCanvas.resizeCanvas(Main.drawScroll.getViewport().getWidth(),Main.drawScroll.getViewport().getHeight());
  	
  	 int loading = 0;
  	if(!(VariantHandler.hideIndels.isSelected() && VariantHandler.hideSNVs.isSelected())) { 
  		
  	 for(int i = sampletemp; i < Main.drawCanvas.sampleList.size(); i++) {
  		 linecounter = 0;
  		if(cancelfileread) {
			cancelFileRead();
			break;
		}
  		 if(( Main.drawCanvas.sampleList.get(i).getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight < Main.drawScroll.getViewport().getHeight()+Main.drawCanvas.drawVariables.sampleHeight) {
 		//	Main.drawCanvas.drawVariables.visibleend = (short)(Main.drawCanvas.sampleList.get(i).getIndex());		
 			Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.drawCanvas.sampleList.get(i).getIndex()+1);	
 			Main.drawCanvas.checkSampleZoom();
 		}
	      Main.drawCanvas.loadbarAll = (int)((loading/(double)(Main.drawCanvas.sampleList.size()-sampletemp))*100);    
	   //   if(!Main.drawCanvas.sampleList.get(i).multipart) {	    	  
	    	  Main.drawCanvas.sampleList.get(i).maxCoverage = 0;  
	      try {
	    	 
	   // 	  vcfreader = new VCFFileReader(new File(Main.drawCanvas.sampleList.get(i).getTabixFile()));
		      tabixreader = new TabixReader(Main.drawCanvas.sampleList.get(i).getTabixFile());
		//      iterator=null;   
		      GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(Main.drawCanvas.sampleList.get(i).getTabixFile()));
		      BufferedReader reader = new BufferedReader(new InputStreamReader(gzip));		
		      line = reader.readLine();
		     
			if(!Main.drawCanvas.sampleList.get(i).multipart && line != null) {
				while(line.startsWith("#")) {
				
					if(line.startsWith("#CHROM")) {
						headersplit = line.split("\t");
						if(headersplit.length -9 > 2) {
							Main.drawCanvas.sampleList.get(i).multiVCF = true;
							Main.varsamples--;
							for(int h = 9; h<headersplit.length; h++) {
								 addSample = new Sample(headersplit[h], (Short)Main.samples, null);        	
					        	 addSample.multipart = true;
					        	 Main.drawCanvas.sampleList.add(addSample);				          	
					          	 Main.samples++;
					          	 Main.varsamples++;
					          
							}
							 VariantHandler.commonSlider.setMaximum(Main.varsamples);
				          	 VariantHandler.commonSlider.setUpperValue(Main.varsamples);
				          	 VariantHandler.geneSlider.setMaximum(Main.varsamples);
						//	Main.drawCanvas.drawVariables.visibleend = (short)(Main.drawCanvas.sampleList.size()-1);
							Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.drawCanvas.sampleList.size());	
							Main.drawCanvas.checkSampleZoom();
						  	Main.drawCanvas.resizeCanvas(Main.drawScroll.getViewport().getWidth(),Main.drawScroll.getViewport().getHeight());
						}			
						
					}					
					
					line = reader.readLine();
				}
			line = reader.readLine();
			if(line != null && line.startsWith("chr")) {
				Main.drawCanvas.sampleList.get(i).vcfchr = "chr";
			}
			line = null;
			reader.close();
			gzip.close();
			}
			else {
				reader.close();
				gzip.close();
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
			  Main.drawCanvas.variantsStart = 1;
			  Main.drawCanvas.variantsEnd = Main.drawCanvas.splits.get(0).chromEnd;
			  getVariants(Main.chromosomeDropdown.getSelectedItem().toString(), Main.drawCanvas.variantsStart,Main.drawCanvas.variantsEnd, Main.drawCanvas.sampleList.get(i));
				 
		  }		 
	  }
	  catch(Exception e) {
		  JOptionPane.showMessageDialog(Main.chromDraw, e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
		  ErrorLog.addError(e.getStackTrace());
		  e.printStackTrace();
		  continue;
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
	Main.opensamples.setText("Open samples");
  	annotate();
  	 tabixreader = null;
    readFiles = false;
 //   Main.drawCanvas.clusterCalc = true;
    if(Control.controlData.controlsOn) {
    	Control.control(head);
    }
    if(Main.drawCanvas.splits.get(0).viewLength < 2000000) {
		
		Main.chromDraw.updateExons = true;
		Main.chromDraw.repaint();
	}
//	Main.drawCanvas.drawVariables.visibleend = Main.samples;
	Main.drawCanvas.drawVariables.visiblesamples = Main.samples;
	Main.drawCanvas.checkSampleZoom();
	Main.drawCanvas.current = head;
  //	Draw.updatevars = true;
//	Main.drawCanvas.repaint();	
	current = null;
	Main.drawCanvas.ready("Loading variants...");
	}
	catch(Exception e)  {
		e.printStackTrace();
		ErrorLog.addError(e.getStackTrace());
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
				
				while(current != null && current.getPosition() < gene.getStart()) {
					if(!current.isInGene()) {
						if(current.getTranscripts() == null) {
							current.setTranscripts();
						}					
						current.getTranscripts().add(prevGene.getTranscripts().get(0));
						current.getTranscripts().add(gene.getTranscripts().get(0));
					}
					current = current.getNext();
				}
				if(gene.getEnd() > prevGene.getEnd()) {
					prevGene = gene;
				}
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
			}						
			while(current != null) {
				if(!current.isInGene()) {
					if(current.getTranscripts() == null) {
						current.setTranscripts();
					}					
					current.getTranscripts().add(prevGene.getTranscripts().get(0));
					current.getTranscripts().add(prevGene.getTranscripts().get(0));
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
		
		pos = Integer.parseInt(split[1])-1;		
		if(pos < Main.drawCanvas.variantsStart) {
			return;
		}
		else if(pos >=Main.drawCanvas.variantsEnd) {
			stop = true;
			return;
		}
		else if(sample.getVCFInput().getFilePointer() > sample.vcfEndPos) {
			return;
		}
		if(linecounter % 100 == 0) {
		//	if(search)  {
				 Main.drawCanvas.loadBarSample = (int)(((pos-Main.drawCanvas.variantsStart)/(double)(Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart))*100);     
				 Draw.updatevars = true;
				 Main.drawCanvas.repaint();
		/*	}
			else {
				 Main.drawCanvas.loadBarSample = (int)((pos/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);     
		
				Draw.updatevars = true;
				 Main.drawCanvas.repaint();
			}*/
		}
		
		linecounter++;
		first = true;
		while(current != null && current.getNext() != null && current.getPosition() < pos ){			
			current = current.getNext();		      				
		}
		for(int s = 9; s < split.length;s++) {
			try {	
				samplecount++;
				currentSample = Main.drawCanvas.sampleList.get(sample.getIndex() +samplecount);			
				info = split[s].split(":");								
				noref = false;
				
				if(split[8].contains("GT")) {	
					
					if(info[split[8].indexOf("GT")/3].contains(".") || info[split[8].indexOf("GT")/3].length() < 3) {
						continue;
					}
					firstallele = Short.parseShort(""+info[split[8].indexOf("GT")/3].charAt(0));
					secondallele = Short.parseShort(""+info[split[8].indexOf("GT")/3].charAt(2));			
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
				if(altbase.contains("*")) {
					continue;
				}
				
				if(split[8].contains("Q")) {
					if(split[8].contains("GQ") && !info[split[8].indexOf("GQ")/3].equals(".")) {						
						quality = (short)Double.parseDouble(info[split[8].indexOf("GQ")/3]);						
					}
					else if(split[8].contains(":BQ:") && !info[split[8].indexOf("BQ")/3].equals(".")) {						
						quality = (short)Double.parseDouble(info[split[8].indexOf("BQ")/3]);						
					}					
					else if(split[7].contains("SSC")) {
						quality = Short.parseShort(split[7].substring(split[7].indexOf("SSC") +4).split(";")[0]);							
					}
					else {						
						quality = (short)20;						
					}
				}
				else if(split[5].matches("\\d+.?.*")) {
					quality = (short)Double.parseDouble(split[5]);
				}
				else {		
					quality = (short)20;					
				}	
				if(VariantHandler.freeze.isSelected()) {
					if(refcalls+altcalls < VariantHandler.coverageSlider.getValue()) {
						continue;
					}
					if(quality < VariantHandler.qualitySlider.getValue()) {
						continue;
					}
					if(altcalls/(double)(refcalls+altcalls) < VariantHandler.callSlider.getValue()/100.0) {
						continue;
					}
					if(VariantHandler.hideSNVs.isSelected() && split[3].length() == 1 && split[4].length() == 1) {
						continue;
					}
					if(VariantHandler.hideIndels.isSelected() && split[3].length() > 1 || split[4].length() > 1) {
						continue;
					}
					if(VariantHandler.rscode.isSelected() && !split[2].equals(".")) {
						continue;
					}					
				}				
				
				if(currentSample.maxCoverage < calls1+calls2) {
					currentSample.maxCoverage = (short)(calls1+calls2);
				}
					
					if(first && current.getNext() == null) {							
									
							if(!split[2].equals(".")) {
								current.putNext(new VarNode(pos, ""+split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, split[2], currentSample, current, null));		
							}
							else {
								current.putNext(new VarNode(pos, ""+split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, null, currentSample, current, null));		
							}
							if(noref ) {
								current.getNext().addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, currentSample);
							}						
							
							current = current.getNext();	
							first = false;
								
					}						
					else if(pos == current.getPosition()){		
						
							current.addSample(altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, currentSample);
							if(noref ) {
								current.addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, currentSample);
							}
							if(current.isRscode() == null && !split[2].equals(".")) {							
								current.setRscode(split[2]);
							}					
							
							
					}
				
					else if(current.getPosition() > pos) {
							if(!split[2].equals(".")) {
								current.getPrev().putNext(new VarNode(pos, ""+split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, split[2], currentSample, current.getPrev(), current));
							}
							else {
								current.getPrev().putNext(new VarNode(pos, ""+split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, null, currentSample, current.getPrev(), current));
								
							}
							if(noref ) {
								current.getPrev().addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, currentSample);
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
	
	public void readLine(String[] split, Sample sample) {		
		try {
	
		if(split.length < 3) {			
			return;
		}	
		
		if(split[4].equals("*")) {
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
		if(linecounter % 100 == 0) {
			//	if(search)  {
					 Main.drawCanvas.loadBarSample = (int)(((pos-Main.drawCanvas.variantsStart)/(double)(Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart))*100);     
					 Draw.updatevars = true;
					 Main.drawCanvas.repaint();
			/*	}
				else {
					 Main.drawCanvas.loadBarSample = (int)((pos/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);     
			
					Draw.updatevars = true;
					 Main.drawCanvas.repaint();
				}*/
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
		
		*/	info = split[split.length-1].split(":");
		//}
		noref = false;
	
		if(split[8].contains("GT")) {	
			
			if(info[split[8].indexOf("GT")/3].contains(".")) {
				return;
			}
			
			firstallele = Short.parseShort(""+info[split[8].indexOf("GT")/3].charAt(0));
			secondallele = Short.parseShort(""+info[split[8].indexOf("GT")/3].charAt(2));			
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
				dpindex = split[7].indexOf("DP4");
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
		
		
		if(split[8].contains("Q")) {
			
			if(split[8].contains("GQ") && !info[split[8].indexOf("GQ")/3].equals(".")) {
				
				quality = (short)Double.parseDouble(info[split[8].indexOf("GQ")/3]);
				
			}
			else if(split[8].contains(":BQ:") && !info[split[8].indexOf("BQ")/3].equals(".")) {
				
				quality = (short)Double.parseDouble(info[split[8].indexOf("BQ")/3]);
				
			}
			
			else if(split[7].contains("SSC")) {
				quality = Short.parseShort(split[7].substring(split[7].indexOf("SSC") +4).split(";")[0]);				
					
			}
			else {
				
				quality = (short)20;
				//quality = (short)Double.parseDouble(info[split[8].indexOf("Q")/3]);
			}
		}
		else if(split[5].matches("\\d+.?.*")) {
			quality = (short)Double.parseDouble(split[5]);
			
		}
		else {		
			quality = (short)20;
			
		}	
		
		if(VariantHandler.freeze.isSelected()) {
			if(refcalls+altcalls < VariantHandler.coverageSlider.getValue()) {				
				return;
			}
			if(quality < VariantHandler.qualitySlider.getValue()) {				
				return;
			}
			if(altcalls/(double)(refcalls+altcalls) < VariantHandler.callSlider.getValue()/100.0) {
				return;
			}
			
			if(VariantHandler.hideSNVs.isSelected() && (split[3].length() == 1 && split[4].length() == 1)) {
				return;
			}
			
			if(VariantHandler.hideIndels.isSelected() && (split[3].length() > 1 || split[4].length() > 1)) {
				return;
			}
			
			if(VariantHandler.rscode.isSelected() && !split[2].equals(".")) {
				return;
			}			
		}
		
			
		
		if(linecounter % 100 == 0) {
			if(search)  {
			
				 Main.drawCanvas.loadBarSample = (int)(((pos-searchStart)/(double)(searchEnd-searchStart))*100);     
				 Draw.updatevars = true;
				 if(!Main.drawCanvas.loading) {
						Draw.calculateVars = true;
					}
				 Main.drawCanvas.repaint();
			}
			else {
				 Main.drawCanvas.loadBarSample = (int)((pos/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);     
				 Draw.updatevars = true;
				 if(!Main.drawCanvas.loading) {
						Draw.calculateVars = true;
					}
				 Main.drawCanvas.repaint();
			}
		}
		
		linecounter++;
		
		if(sample.maxCoverage < refcalls+altcalls) {
			sample.maxCoverage = (short)(refcalls+altcalls);
			
		}		
		
		while(current.getNext() != null && current != null && current.getPosition() < pos ){	
		
			current = current.getNext();		      				
		}		
		
		if(current.getNext() == null && current.getPosition() < pos) {	
			
			try {				
				if(!split[2].equals(".")) {
					current.putNext(new VarNode(pos, ""+split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, split[2], sample, current, null));		
				}
				else {
					current.putNext(new VarNode(pos, ""+split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, null, sample, current, null));		
					
				}
				if(noref ) {
					
					current.getNext().addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, sample);
				}
				
		
			}
			catch(Exception e) {
				ErrorLog.addError(e.getStackTrace());
				e.printStackTrace();
				
			}
			
		}		
		
		else if(current.getPosition() == pos) {
				
				current.addSample(altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, sample);
				
				if(noref ) {
					current.addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, sample);
				}
				if(current.isRscode() == null && !split[2].equals(".")) {							
					current.setRscode(split[2]);
				}
			
		}
		else if(current.getPosition() > pos) {
			
				if(!split[2].equals(".")) {
					current.getPrev().putNext(new VarNode(pos, ""+split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, split[2], sample, current.getPrev(), current));
				}
				else {
					current.getPrev().putNext(new VarNode(pos, ""+split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, null, sample, current.getPrev(), current));
					
				}
				
				current.putPrev(current.getPrev().getNext());
				current = current.getPrev();
				if(noref ) {
					current.addSample(altbase2, (short)(refcalls+altcalls), refcalls, genotype, quality, sample);
				}
		}
		
		}
		catch(Exception ex) {
			System.out.println(split[8] +" " +split[10] +" " +split[3] +" " +split[4]);
			ErrorLog.addError(ex.getStackTrace());
			ex.printStackTrace();
			
		}
		
	}
	
	File getListFile(File listfile, String chrom) {
		
		try {
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
			else if(alt.length() == 2) {
					return "ins" +alt.substring(1);					
			}
			else if(ref.length() > 2) {				
					return "del" +(ref.length()-1);						
			}
			else {				
					return "ins" +(alt.length()-1);				
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
				if(alt.length() -ref.length() == 1) {
					return "ins" +alt.substring(1, alt.length()-(ref.length()-1));
				}
				else {
					return "ins" +(alt.length() -ref.length());
				}
				
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
				File file = getListFile(readClass.sample.samFile,readClass.sample.chr +  chrom);
				if(file != null) {
					samFileReader = SamReaderFactory.make().open(file);
				}
			}
			else {
				if(readClass.sample.samFile.getName().endsWith("cram")) {					
					
					 CRAMReader = new CRAMFileReader(readClass.sample.samFile, new File(readClass.sample.samFile.getCanonicalFile() +".crai"),new ReferenceSource(Main.ref),ValidationStringency.SILENT);
					
				}
				else {
					samFileReader = SamReaderFactory.make().open(readClass.sample.samFile);
				}
			}
			if(samFileReader != null && !samFileReader.hasIndex()) {				
				JOptionPane.showMessageDialog(null, "Index file is missing (.bai) for " +readClass.sample.samFile.getName(), "Error", JOptionPane.ERROR_MESSAGE);		    		
				return null;
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}		
		try {			
			if(readClass.sample.samFile.getName().endsWith("cram")) {
				QueryInterval[] interval = { new QueryInterval(CRAMReader.getFileHeader().getSequence(chrom).getSequenceIndex(), startpos, endpos) };
				return CRAMReader.query(interval, false);
			}
			else {
				return samFileReader.queryOverlapping(chrom, startpos, endpos);	
			}
		}
		catch(Exception e) {			
			try {
				if(chrom.contains("M")) {
					return samFileReader.queryOverlapping(readClass.sample.chr + "M", startpos, endpos);	
				}
			}
			catch(Exception ex) {
				ex.printStackTrace();
				ErrorLog.addError(e.getStackTrace());
				return null;
			}
			
		}
		return null;
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
		
		if(Main.drawCanvas.drawVariables.sampleHeight > 100) {
			if(viewLength > Settings.readDrawDistance) {
				
				splitIndex.readSequence = null;
				
				if(firstCov || readClass.getCoverageStart() == 0 || (readClass.getCoverageStart() < startpos && readClass.getCoverageEnd()  > endpos) || (readClass.getCoverageStart() > startpos && readClass.getCoverageEnd() < endpos) ) {
					
					double[][] coverages = new double[(int)Main.frame.getWidth()][8];
				//	System.out.println(coverages.length);
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
					
					try {
						bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos-searchwindow,  endpos+searchwindow  );
						
					}
					catch(Exception e) {
						try {
							
						}
						catch(java.lang.IllegalArgumentException ex) {
							readClass.loading = false;
							System.out.println("Chromosome " +chrom +" not found in " +readClass.sample.samFile.getName());
							ErrorLog.addError(ex.getStackTrace());
							return null;
						}
						
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
			
				bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos-searchwindow,  endpos+searchwindow  );
								
				startpos = readClass.getReadEnd();				
				endpos = endpos+searchwindow;
				addlength = 1000;
				
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
	
				
				bamIterator = getBamIterator(readClass,readClass.sample.chr + chrom,startpos-searchwindow,  endpos+searchwindow  );
				endpos = readClass.getReadStart();
				
				startpos = startpos-searchwindow;
				
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
				readClass.setCoverages(coverages);				
			}			

			
			
			}
			
		//	System.out.println(splitIndex.readSeqStart +" " +splitIndex.readSequence.length());
			coverageArray = readClass.getCoverages();			
			
			while(bamIterator != null && bamIterator.hasNext()) {	
			
				try {	
					
					if(cancelreadread) {
		//				readrunner.end = true;						
			  			return null;
			  		}
					try {
						samRecord = bamIterator.next(); 	
						
					}
					catch(htsjdk.samtools.SAMFormatException ex) {
					
						
					}
					
				if(samRecord.getReadUnmappedFlag()) {					
					continue;
				}
				//TODO jos on pienempi ku vika unclipped start		
				
				if(left && samRecord.getUnclippedStart() >= endpos) {
					Main.drawCanvas.loadBarSample = 0;  
					Main.drawCanvas.loadbarAll  =0;
					setLevels(lastAdded, readClass.sample, readClass);					
					lastAdded = null;
					break;
				}
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
				
			if(viewLength < Settings.readDrawDistance) {
				 if(splitIndex.readSequence != null && samRecord.getUnclippedEnd() > splitIndex.readSeqStart+splitIndex.readSequence.length()-1) {
					splitIndex.readSequence.append(Main.chromDraw.getSeq(chrom, splitIndex.readSeqStart+splitIndex.readSequence.length(), samRecord.getUnclippedEnd()+(int)splitIndex.viewLength, Main.referenceFile));
				 
				 }
				
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
					java.util.ArrayList<java.util.Map.Entry<Integer,String>> mismatches = null;
					
					if(samRecord.getMappingQuality() > Settings.getMappingQuality()) {
						
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
							
							mismatches = getMismatches(samRecord, readClass, coverageArray);
							
						}
					}
					
					if(viewLength <= Settings.readDrawDistance) {	
						
						if(right && readClass.getLastRead() != null && (samRecord.getUnclippedStart() <= startpos /* || samRecord.getReadName().equals(readClass.getLastRead().getName())*/)) {
							
							continue;
						}
						
						if(right && samRecord.getUnclippedStart() > endpos) {
							right = false;
							
						}
						
						if(left) {
						
							if(readClass.getFirstRead().getPosition() > samRecord.getUnclippedStart()) {							
								
								ReadNode addNode = new ReadNode(samRecord, readClass.complete, chrom, splitIndex.readSequence, splitIndex.readSeqStart, readClass.sample, splitIndex,readClass, mismatches);
								readClass.setFirstRead(addNode);
								startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-(readClass.readHeight+2)));					
						//		addNode.setRectBounds((int)((addNode.getPosition()-start)*pixel), startY, (int)(samRecord.getReadLength()*pixel), Main.drawCanvas.readHeight);
														
								addNode.setPrev(null);
								addNode.setNext(readClass.getHeadAndTail().get(0)[headnode]);
								readClass.getHeadAndTail().get(0)[headnode].setPrev(addNode);
								lastAdded = addNode;
						
							}
							else {
								
								ReadNode addNode = new ReadNode(samRecord, readClass.complete, chrom, splitIndex.readSequence, splitIndex.readSeqStart, readClass.sample, splitIndex,readClass, mismatches);
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
			//					readrunner.end = true;
								
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
		
	//		current = null;
	//		readHead = null;
	//		readrunner.end = true;
			
		}
		
		}
		catch(Exception ex) {	
			if(readClass.getReadStart() > startpos) {
				readClass.setReadStart(startpos);
			}
			if(readClass.getReadEnd() < endpos) {
				readClass.setReadEnd(endpos);
			}
			
			
			ErrorLog.addError(ex.getStackTrace());
			ex.printStackTrace();		
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
	
	readClass.setCoverageEnd(readClass.getCoverageStart()+coverages.length);
	readClass.setCoverages(coverages);
	return coverages;
	
}
	
	
java.util.ArrayList<java.util.Map.Entry<Integer,String>> getMismatches(SAMRecord samRecord, Reads readClass, double[][] coverageArray) {
	
	
	java.util.ArrayList<java.util.Map.Entry<Integer,String>> mismatches = null;
	String MD = samRecord.getStringAttribute("MD");	
	String SA = samRecord.getStringAttribute("SA");
	
	try {
	
	if(MD == null) {
		
		if(readClass.sample.MD) {
			readClass.sample.MD = false;
		}
	}
	if((readClass.sample.MD || MD!=null) && ((!samRecord.getCigarString().contains("S") && SA ==null) || SA !=null)) {		
		
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
						readClass.getCoverages()[(readstart+readpos+d)- readClass.getCoverageStart()][Main.baseMap.get("D")]++;
						
					}
					readpos+=element.getLength();
				}
				else if(element.getOperator().compareTo(CigarOperator.INSERTION)== 0) {	
					if(insertions == null) {
						insertions = new java.util.ArrayList<java.util.Map.Entry<Integer,Integer>>();						
					}
					java.util.Map.Entry<Integer, Integer> ins = new java.util.AbstractMap.SimpleEntry<>(readpos, element.getLength());
					insertions.add(ins);
					readClass.getCoverages()[(readstart+readpos)- readClass.getCoverageStart()][Main.baseMap.get("I")]++;
		
					mispos+=element.getLength();
					
				
				}
				else if(element.getOperator().compareTo(CigarOperator.SOFT_CLIP)== 0 ) {
					
						if(i==0 && SA == null) {						
								softstart = element.getLength();
						}				
					
						if(Settings.softClips) {
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
							mismatches = new java.util.ArrayList<java.util.Map.Entry<Integer,String>>();
						}
						
						readClass.getCoverages()[(readstart+readpos)- readClass.getCoverageStart()][Main.baseMap.get(""+bases.charAt(mispos))]++;
						
							if(samRecord.getBaseQualityString().length() != 1 && (int)qualities.charAt(mispos)-readClass.sample.phred < Settings.getBaseQuality() ) {
								
								java.util.Map.Entry<Integer, String> mismatch = new java.util.AbstractMap.SimpleEntry<>(readpos, ""+Character.toLowerCase(bases.charAt(mispos)));
								mismatches.add(mismatch);
							}
							else {
								
								java.util.Map.Entry<Integer, String> mismatch = new java.util.AbstractMap.SimpleEntry<>(readpos, ""+bases.charAt(mispos));
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
		if(splitIndex.readSequence == null) {
			splitIndex.readSequence = Main.chromDraw.getSeq(splitIndex.chrom, (int)(splitIndex.start-searchwindow-1), (int)(splitIndex.end+splitIndex.viewLength+200), Main.referenceFile);
			splitIndex.readSeqStart = (int)(splitIndex.start-searchwindow-1);	
		}
		if(samRecord.getCigarLength() > 1) {			
			readstart = samRecord.getUnclippedStart();			
			readpos= 0;
			mispos= 0;
			if(readClass.getCoverageStart()>readstart) {
				return null;
			}
			if(mismatches == null) {
				mismatches = new java.util.ArrayList<java.util.Map.Entry<Integer,String>>();
			}
			for(int i = 0; i<samRecord.getCigarLength(); i++) {
				element = samRecord.getCigar().getCigarElement(i);
				if(element.getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH)== 0) {
					for(int r = readpos; r<readpos+element.getLength(); r++) {																	
						if((readstart+r)- splitIndex.readSeqStart-1 < 0) {
							while((readstart+r)- splitIndex.readSeqStart-1 < 0) {
								splitIndex.readSequence.insert(0, Main.chromDraw.getSeq(splitIndex.chrom, (int)(splitIndex.readSeqStart-1000), splitIndex.readSeqStart, Main.referenceFile));
								splitIndex.readSeqStart =  (int)(splitIndex.readSeqStart-1000);				
							}
						}
					
						if(samRecord.getReadString().charAt(mispos) !=  splitIndex.readSequence.charAt((readstart+r)-splitIndex.readSeqStart-1)) {	
							readClass.getCoverages()[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get(""+samRecord.getReadString().charAt(mispos))]++;
							
							if(samRecord.getBaseQualityString().length() != 1 && (int)samRecord.getBaseQualityString().charAt(mispos)-readClass.sample.phred < Settings.getBaseQuality() ) {
								
								java.util.Map.Entry<Integer, String> mismatch = new java.util.AbstractMap.SimpleEntry<>(r, ""+Character.toLowerCase(samRecord.getReadString().charAt(mispos)));
								mismatches.add(mismatch);
							}
							else {
								
								java.util.Map.Entry<Integer, String> mismatch = new java.util.AbstractMap.SimpleEntry<>(r, ""+samRecord.getReadString().charAt(mispos));
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
						readClass.getCoverages()[(readstart+readpos+d)- readClass.getCoverageStart()][Main.baseMap.get("D")]++;
						
					}
					readpos+=element.getLength();
				}
				else if(element.getOperator().compareTo(CigarOperator.INSERTION)== 0) {										
					readClass.getCoverages()[(readstart+readpos)- readClass.getCoverageStart()][Main.baseMap.get("I")]++;
			//		readClass.getCoverages()[(readstart+mispos)- readClass.getCoverageStart()][Main.baseMap.get("I")]++;
				/*	if(mismatches == null) {
						mismatches = new java.util.ArrayList<java.util.Map.Entry<Integer,String>>();
					}
					java.util.Map.Entry<Integer, String> mismatch = new java.util.AbstractMap.SimpleEntry<>(readpos, "I " +element.getLength());
					
					mismatches.add(mismatch);
					*/
					mispos+=element.getLength();
					
				
				}
				else if(element.getOperator().compareTo(CigarOperator.SOFT_CLIP)== 0 ) {
						if(Settings.softClips) {
							for(int r = readpos; r<readpos+element.getLength(); r++) {
								
								if(samRecord.getReadString().charAt(mispos) !=  splitIndex.readSequence.charAt((readstart+r)- splitIndex.readSeqStart-1)) {													
	
									if(mismatches == null) {
										mismatches = new java.util.ArrayList<java.util.Map.Entry<Integer,String>>();
									}
									java.util.Map.Entry<Integer, String> mismatch = new java.util.AbstractMap.SimpleEntry<>(r, ""+samRecord.getReadString().charAt(mispos));													
									mismatches.add(mismatch);												
								}
								if(SA == null) {
									coverageArray[((readstart+r)-readClass.getCoverageStart())][0]++; 
									readClass.getCoverages()[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get(""+samRecord.getReadString().charAt(mispos))]++;
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
			for(int i = 0; i<samRecord.getReadLength(); i++) {
				try {									
					if((readstart+i)- splitIndex.readSeqStart-1 < 0) {
						while((readstart+i)- splitIndex.readSeqStart-1 < 0) {
							splitIndex.readSequence.insert(0, Main.chromDraw.getSeq(splitIndex.chrom, (int)(splitIndex.readSeqStart-1000), splitIndex.readSeqStart, Main.referenceFile));
							splitIndex.readSeqStart =  (int)(splitIndex.readSeqStart-1000);	
						
						}
					}
					if((readstart+i)- splitIndex.readSeqStart-1 > splitIndex.readSequence.length()-1) {
						while((readstart+i)- splitIndex.readSeqStart-1 > splitIndex.readSequence.length()-1) {
							splitIndex.readSequence.append(Main.chromDraw.getSeq(splitIndex.chrom, splitIndex.readSeqStart+splitIndex.readSequence.length(), splitIndex.readSeqStart+splitIndex.readSequence.length()+1000, Main.referenceFile));
							
						
						}
					}
					if(samRecord.getReadString().charAt(i) != splitIndex.readSequence.charAt((readstart+i)- splitIndex.readSeqStart-1)) {									
						if((readstart+i)-readClass.getCoverageStart() < readClass.getCoverages().length-1) {
							readClass.getCoverages()[(readstart+i)- readClass.getCoverageStart()][Main.baseMap.get(""+samRecord.getReadString().charAt(i))]++;
							if(mismatches == null) {
								mismatches = new java.util.ArrayList<java.util.Map.Entry<Integer,String>>();
							}
							java.util.Map.Entry<Integer, String> mismatch = new java.util.AbstractMap.SimpleEntry<>(i, ""+samRecord.getReadString().charAt(i));												
							mismatches.add(mismatch);
						}									
					}										
					try {
						coverageArray[(readstart+i)-readClass.getCoverageStart()][0]++; 
					}
					catch(Exception e) {
					//	System.out.println(readClass.getCoverageStart());
						e.printStackTrace();
					}
					if(coverageArray[(samRecord.getUnclippedStart()-readClass.getCoverageStart())+i][0] > readClass.getMaxcoverage()) {
						readClass.setMaxcoverage(coverageArray[(samRecord.getUnclippedStart()-readClass.getCoverageStart())+i][0]);
					}
				}
				catch(Exception e) {
				//	System.out.println(splitIndex.readSeqStart +" " +splitIndex.readSequence.length());
					ErrorLog.addError(e.getStackTrace());
					e.printStackTrace();
					FileRead.cancelreadread = true;
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
		 
		  if(!Main.trackPane.isVisible()) {
			  Main.trackPane.setVisible(true);
			  Main.varpane.setDividerSize(4);	  
			  Main.varpane.setResizeWeight(files.length *0.1);
			  Main.varpane.setDividerLocation(files.length *0.1);			  
			  Main.varpane.revalidate();
		  }
		  else {
			  Main.varpane.setDividerLocation(Main.varpane.getDividerLocation()+100);
			  Main.varpane.revalidate();
			  Main.trackPane.setDividerLocation(Main.varpane.getDividerLocation()/2);			  
			  
			  if(Main.controlScroll.isVisible()) {
				  Main.trackPane.setDividerSize(4);
			  }
		  }
		  Main.bedScroll.setVisible(true);
		  Main.bedCanvas.setVisible(true);
		  
		 for(int i = 0; i<files.length; i++) {
		   
			  bedfile = files[i];		   
		      BedTrack addTrack = new BedTrack(bedfile,Main.bedCanvas.bedTrack.size());	
		      if(bedfile.getName().toLowerCase().endsWith(".bedgraph") || bedfile.getName().toLowerCase().endsWith(".bedgraph.gz") || bedfile.getName().toLowerCase().endsWith("bigwig") || bedfile.getName().toLowerCase().endsWith("bw")) {
		    	  addTrack.graph = true;
		      }
		      Main.bedCanvas.bedTrack.add(addTrack);
		      Main.bedCanvas.trackDivider.add(0.0);
		      if(bedfile.length() /1048576 < 200) {
		    	  addTrack.small = true;		    	
		    	  Main.bedCanvas.getBEDfeatures(addTrack, 1, Main.drawCanvas.splits.get(0).chromEnd);		    		    	  
		    	 
		      }	
		      setTable(addTrack);
		 }		
		
		 for(int i = 0 ; i<Main.bedCanvas.trackDivider.size(); i++) {
			 Main.bedCanvas.trackDivider.set(i, ((i+1)*(Main.varpane.getDividerLocation()/(double)Main.bedCanvas.trackDivider.size())));
		 }		
		}
		catch(Exception e) {			
			e.printStackTrace();
			ErrorLog.addError(e.getStackTrace());
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
static void readGFF(File infile, String outfile) {
	
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
		line = reader.readLine();
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
		    		if(!lineHash.get("type").equals("gene")) {
		    			
		    			continue;
		    		}
		    		
		    		
		    		genes.put(getInfoValue(lineHash, "id"),new Gene(chrom, lineHash));
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
	    MethodLibrary.blockCompressAndIndex(geneArray, outfile, false);	  
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

static void readGFF2(File infile, String outfile) {	
	
		try {
			BufferedReader reader = null;
			GZIPInputStream gzip = null;
		
			if(infile.getName().endsWith(".gz")) {
				gzip = new GZIPInputStream(new FileInputStream(infile));
				reader = new BufferedReader(new InputStreamReader(gzip));	
				
			}
			else {
				reader = new BufferedReader(new FileReader(infile));	
			}
	//		String header = "#Chrom\tGeneStart\tGeneEnd\tName\tExonCount\tStrand\tENSG\tENST\tUniProt\tCanonical\tBiotype\tCodingStart\tCodingEnd\tExonStarts\tExonEnds\tStartPhases\tDescription\n";
	//		System.out.println(header);
			
			String line, chrom = "-1";
			StringBuffer exonStarts, exonEnds, exonPhases;
			ArrayList<Integer> eStarts = new ArrayList<Integer>(), eEnds = new ArrayList<Integer>(), ePhases = new ArrayList<Integer>();
			HashMap<String, String> parentHash, childHash, tempChildHash, subParentHash;
			ArrayList<String[]> rows = new ArrayList<String[]>();
			HashMap<String,Gene> genes = new HashMap<String, Gene>();
			HashMap<String,Transcript> transcripts = new HashMap<String, Transcript>();
			String[] split, infosplit, resultTable;
			int count = 0, index, exonpointer;			
			line = reader.readLine();
			boolean first = true;
			String[] firstrow = {"0", "0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0"};
			rows.add(firstrow);
			
	    while(true) { 
	    	
	    	if(line == null) {
	    		break;
	    	}
	    /*	count++;
	    	if(count > 300) {
	    		
	    		break;
	    		//System.out.println(count);
	    		
	    	}*/
	    	if(line.startsWith("#")) {
	    		line = reader.readLine();
	    		continue;
	    	}
	    	
	    	parentHash = makeHash(line.split("\t"));
	    	if(parentHash == null) {
	    		line = reader.readLine();
	    		continue;
	    	}
	    	if(parentHash.get("name") == null) {
	    		line = reader.readLine();
	    		continue;
	    	}
	    	if(parentHash.get("type").startsWith("region")) {
	    		if(parentHash.get("type").startsWith("region") && parentHash.get("chromosome") != null) {
	    			if(parentHash.get("name") != null) {
	    				chrom = parentHash.get("name").replace("chr", "");	    				
	    			}	    			
	    		}
	    		
	    		line = reader.readLine();
	    		continue;
	    	}
	    
	    	if(!parentHash.containsKey("parent")) {
	    		line = reader.readLine();
	    		
	    		tempChildHash = makeHash(line.split("\t"));
	    		if(tempChildHash != null && !tempChildHash.containsKey("parent")) {
	    			childHash = tempChildHash;
	    			resultTable = new String[17];
	    			if(!chrom.equals("-1")) {
	    				resultTable[0] = chrom;
	    			}
	    			else {
	    				resultTable[0] = getInfoValue(parentHash,"seqid");	    			
	    			}
	    			resultTable[1] = getInfoValue(childHash,"start");
					resultTable[2] = getInfoValue(childHash,"end");
					resultTable[3] = getInfoValue(parentHash,"name");
					resultTable[4] = "1";
					resultTable[5] = getInfoValue(parentHash,"strand");
					resultTable[6] = getInfoValue(parentHash,"id");
					
					resultTable[7] = getInfoValue(childHash,"id");
					resultTable[8] = getInfoValue(childHash,"uniprot");
					resultTable[9] =  getInfoValue(childHash,"canonical");;
					resultTable[10] = getInfoValue(childHash,"biotype");
					resultTable[11] = getInfoValue(childHash,"end");
					resultTable[12]	= getInfoValue(childHash,"end");
					resultTable[13] = resultTable[1]+",";
					resultTable[14] = resultTable[2]+",";
					if(!parentHash.containsKey("description")) {
						resultTable[16] = getInfoValue(parentHash,"note").replace("%20", " ");
					}
					else {
						resultTable[16] = getInfoValue(parentHash,"description").replace("%20", " ");
					}
					if(childHash.get("phase").equals(".")) {
	    				
						resultTable[15] = "-1";
	    			}
	    			else {
	    				if(childHash.get("phase").equals("0")) {
	    					resultTable[15] = "0";
    					}
    					else if(childHash.get("phase").equals("1")) {
    						resultTable[15] = "2";
    					}
    					else {
    						resultTable[15] = "1";
    					}
    					
	    			}
					
					rows.add(resultTable);
				/*	String writeString;
					for(int i = 0; i<17; i++) {
						writeString = resultTable[i] +"\t";
						writer.write(writeString.getBytes());
					
					}
					writer.write("\n".getBytes());*/
			//		System.out.println();
	    			continue;
	    		}
	    		childHash = tempChildHash;
	    		while(childHash !=null && childHash.containsKey("parent")) {
	    		
	    			if(childHash.get("parent").contains(parentHash.get("id"))) {
	    			
	    				resultTable = new String[17];
	    				if(!chrom.equals("-1")) {
		    				resultTable[0] = chrom;
		    			}
		    			else {
		    				resultTable[0] = getInfoValue(parentHash,"seqid");	    			
		    			}
						resultTable[1] = getInfoValue(childHash,"start");
						resultTable[2] = getInfoValue(childHash,"end");
						resultTable[3] = getInfoValue(parentHash,"name");
						resultTable[4] = "1";
						resultTable[5] = getInfoValue(parentHash,"strand");
						resultTable[6] = getInfoValue(parentHash,"id");
						
						resultTable[7] = getInfoValue(childHash,"id");
						resultTable[8] = getInfoValue(childHash,"uniprot");
						resultTable[9] =  getInfoValue(childHash,"canonical");;
						resultTable[10] = getInfoValue(childHash,"biotype");
						resultTable[11] = getInfoValue(childHash,"end");
						resultTable[12]	= getInfoValue(childHash,"end");
						resultTable[13] = childHash.get("start");
						resultTable[14] = childHash.get("end");
						resultTable[15] = "-1,";
						if(!parentHash.containsKey("description")) {
							resultTable[16] = getInfoValue(parentHash,"note").replace("%20", " ");
						}
						else {
							resultTable[16] = getInfoValue(parentHash,"description").replace("%20", " ");
						}
						if(childHash.containsKey("id")) {
							subParentHash = childHash;					
						}
						else {
							subParentHash = parentHash;
						}
						line = reader.readLine();
						
			    		tempChildHash = makeHash(line.split("\t"));
			    		if(tempChildHash != null && tempChildHash.containsKey("parent") && !tempChildHash.get("parent").contains(subParentHash.get("id"))) {
			    			
			    			resultTable[13] = childHash.get("start");
							resultTable[14] = childHash.get("end");
							
			    			if(childHash.get("phase").equals(".")) {
			    				
								resultTable[15] = "-1";
			    			}
			    			else {
			    				if(childHash.get("phase").equals("0")) {
			    					resultTable[15] = "0";
		    					}
		    					else if(childHash.get("phase").equals("1")) {
		    						resultTable[15] = "2";
		    					}
		    					else {
		    						resultTable[15] = "1";
		    					}
			    				resultTable[11] = getInfoValue(childHash,"start");
			    				
			    			}
			    			rows.add(resultTable);
			    			/*
			    			String writeString;
							for(int i = 0; i<17; i++) {
								writeString = resultTable[i] +"\t";
								writer.write(writeString.getBytes());
							
							}
							writer.write("\n".getBytes());
			    			*/
			    			break;
			    		}
			    		childHash = tempChildHash;
			    		first = true;
			    		eStarts.clear();
			    		eEnds.clear();
			    		ePhases.clear();
			    		exonStarts = new StringBuffer();
			    		exonEnds = new StringBuffer();
			    		exonPhases = new StringBuffer();
			    		
			    		while(childHash != null && childHash.containsKey("parent") && childHash.get("parent").contains(subParentHash.get("id"))) {
			    			
			    			if(childHash.get("type").contains("UTR")) {
			    				line = reader.readLine();
					    		childHash = makeHash(line.split("\t"));
					    		continue;
			    			}
			    			if(childHash.get("phase").equals(".")) {
			    				if(eStarts.size() > 0 && eEnds.get(eEnds.size()-1) > Integer.parseInt(childHash.get("start"))) {
			    					if( eStarts.get(eStarts.size()-1) > Integer.parseInt(childHash.get("start"))) {
			    						eStarts.set(eStarts.size()-1, Integer.parseInt(childHash.get("start")));
			    					}
			    					if(eEnds.get(eEnds.size()-1) < Integer.parseInt(childHash.get("end"))) {
			    						eEnds.set(eEnds.size()-1, Integer.parseInt(childHash.get("end")));
			    					}
			    				}
			    				else {
				    				eStarts.add(Integer.parseInt(childHash.get("start")));
				    				eEnds.add(Integer.parseInt(childHash.get("end")));
				    				ePhases.add(-1);	
				    				
			    				}
			    			}
			    			else {
			    				if(first) {
			    					first = false;
			    					resultTable[11] = childHash.get("start");			    					
			    				}
			    				resultTable[12] = childHash.get("end");
			    				
			    				if(eStarts.size() > 0 && eEnds.get(eEnds.size()-1) > Integer.parseInt(childHash.get("start"))) {
			    					if(childHash.get("phase").equals("0")) {
			    						ePhases.set(ePhases.size()-1, 0);
			    					}
			    					else if(childHash.get("phase").equals("1")) {
			    						ePhases.set(ePhases.size()-1, 2);
			    					}
			    					else if(childHash.get("phase").equals("2")) {
			    						ePhases.set(ePhases.size()-1, 1);
			    					}
			    					else {
			    				//		System.out.println(childHash.get("phase"));
			    						ePhases.set(ePhases.size()-1, -1);
			    					}
			    				}
			    				else {
				    				eStarts.add(Integer.parseInt(childHash.get("start")));
				    				eEnds.add(Integer.parseInt(childHash.get("end")));
				    				if(childHash.get("phase").equals("0")) {
			    						ePhases.add(0);
			    					}
			    					else if(childHash.get("phase").equals("1")) {
			    						ePhases.add(2);
			    					}
			    					else if(childHash.get("phase").equals("2")) {
			    						ePhases.add(1);
			    					}
			    					else {
			    			//			System.out.println(childHash.get("phase"));
			    						ePhases.add(-1);
			    					}
			    				}
			    				
			    			}
			    			line = reader.readLine();
			    			
				    		childHash = makeHash(line.split("\t"));
			    		}	    				
			    		for(int i=0; i<eStarts.size(); i++) {
			    			exonStarts.append(eStarts.get(i) +",");
			    			exonEnds.append(eEnds.get(i) +",");
			    			exonPhases.append(ePhases.get(i)+",");
			    		}
			    		resultTable[4] = ""+eStarts.size();
			    		
			    		resultTable[13] = exonStarts.toString();
						resultTable[14] = exonEnds.toString();
						resultTable[15] = exonPhases.toString();
						rows.add(resultTable);
						/*String writeString;
						for(int i = 0; i<17; i++) {
							writeString = resultTable[i] +"\t";
							writer.write(writeString.getBytes());
						
						}
						writer.write("\n".getBytes());
	    			*/
	    			}    			
	    			
	    			line = reader.readLine();
	    			if(line != null) {
	    				childHash = makeHash(line.split("\t"));
	    			}
	    			else {
	    			
	    				childHash = null;
	    			}
	    		}
	    	}
	    	else {
	    		line = reader.readLine();
	    		continue;
	    	}
	    	
		
	    	
	    }
	   
	    rows.remove(0);
	//    writer.write(header.getBytes());
	    gffSorter gffsorter = new gffSorter();
	    Collections.sort(rows, gffsorter);
	   
	   /* 
	    for(int i = 0 ; i< rows.size(); i++) {
	    	for(int j = 0 ; j<rows.get(i).length; j++) {
	    		writer.write(rows.get(i)[j].getBytes());
	    		if(j < rows.get(i).length-1) {
	    			writer.write('\n');
	    		}
	    	}
	    	if(i < rows.size()-1) {
	    		writer.write('\n');
	    	}
	    }
	    writer.flush();
	    writer.close();
	 //  	TabixIndexCreator indexer = new TabixIndexCreator(TabixFormat.PSLTBL);
	   */
	    
	   MethodLibrary.blockCompressAndIndex(rows, outfile, false);
	    rows.clear();
	    
	   
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
		  Main.varpane.setDividerSize(4);	  
		  Main.varpane.setResizeWeight(0.1);
		  Main.varpane.setDividerLocation(0.1);			  
		 
	  }
	  else {
		  Main.varpane.setDividerLocation(Main.varpane.getDividerLocation()+100);
		  Main.varpane.revalidate();
		  Main.trackPane.setDividerLocation(Main.varpane.getDividerLocation()/2);			  
		  
		  if(Main.controlScroll.isVisible()) {
			  Main.trackPane.setDividerSize(4);
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
	 if(addTrack != null) {
		 if(url.toLowerCase().endsWith(".bedgraph") || url.toLowerCase().endsWith("bigwig")  || url.toLowerCase().endsWith("bw")) {
	    	  addTrack.graph = true;
	      }
		Main.bedCanvas.trackDivider.add(0.0);
		Main.bedCanvas.bedTrack.add(addTrack);
	 }
	 
	}
	catch(Exception e) {			
		e.printStackTrace();
		ErrorLog.addError(e.getStackTrace());
	}  	 
}
static void setTable(BedTrack track) {
	/*if(!track.small) {
		return;
	}*/
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
	  
	  if(track.file.getName().length() > 10) {
		  VariantHandler.tabs.add(track.file.getName().substring(0, 10) +"...", addpane);
	  }
	  else {
		  VariantHandler.tabs.add(track.file.getName(), addpane);
	  }
	  VariantHandler.clusterTable.addHeaderRow(track.file.getName());
	  VariantHandler.tabs.revalidate();
}
static void removeTable(BedTrack track) {
	try {
		
	if(VariantHandler.tables.indexOf(track.getTable()) > -1) {
		
		 VariantHandler.tabs.remove(VariantHandler.tables.indexOf(track.getTable())+3);
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

void readSam(String chrom, Reads readClass, SAMRecord samRecordBuffer, java.util.ArrayList<java.util.Map.Entry<Integer,String>> mismatches) {
	try {
		
		if(readClass.getReads().isEmpty()) {
			
			addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, splitIndex.readSequence, splitIndex.readSeqStart, readClass.sample, splitIndex,readClass, mismatches);			
			readClass.setFirstRead(addNode);			
			startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-(readClass.readHeight+2)));					
			addNode.setRectBounds((int)((addNode.getPosition()-start)*pixel), startY, (int)(samRecordBuffer.getReadLength()*pixel), readClass.readHeight);		
			readClass.getReads().add(addNode);			
			ReadNode[] addList = new ReadNode[2];
			addList[headnode] = addNode;
			addList[tailnode] = addNode;
			readClass.getHeadAndTail().add(addList);
			readClass.setLastRead(addNode);
			
		/*	int[] adder = new int[2];
			adder[0] = 0;
			adder[1] = addNode.getPosition()+addNode.getWidth()+2;
			readClass.getHelpArray().add(adder);
			readClass.getHelpHash().put(0, 0);*/
		}
		else {		
			/*
			if(readClass.getHelpArray().get(0)[1]>=samRecordBuffer.getUnclippedStart()) {
				if(readClass.getHelpArray().size() < 100) {
					xpos = (int)((samRecordBuffer.getUnclippedStart()-start)*pixel);					
					addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, readSequence, readSeqStart, readClass.sample, splitIndex);	
					
					int[] adder = new int[2];
					adder[0] = readClass.getHelpArray().size();
					adder[1] = addNode.getPosition()+addNode.getWidth()+2;
					readClass.getHelpArray().add(adder);
					readClass.getHelpHash().put(readClass.getHelpArray().size(), readClass.getHelpArray().size());
					ReadNode[] addList = new ReadNode[2];
					addList[headnode] = addNode;
					addList[tailnode] = addNode;
					addNode.setPrev(null);
					addNode.setNext(null);
					readClass.getHeadAndTail().add(addList);
					readClass.getReads().add(addNode);
					startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((readClass.getReads().size())*(Main.drawCanvas.readHeight+2))));					
					addNode.setRectBounds(xpos, startY, (int)(addNode.getWidth()*pixel), Main.drawCanvas.readHeight);
				}
				
			}
			else {
				/*for(int i = 0; i<readClass.getHelpArray().size(); i++) {
					if(samRecordBuffer.getUnclippedStart() >readClass.getHelpArray().get(i)[1]) {
						
						index = i;
						break;
					}
				}*/
				
				//Binary search
			/*	xpos = (int)((samRecordBuffer.getUnclippedStart()-start)*pixel);	
				addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, readSequence, readSeqStart, readClass.sample, splitIndex);	
				readClass.setLastRead(addNode);	
				
				for(int i = 0; i<readClass.getHeadAndTail().size(); i++) {
					if(readClass.getHeadAndTail().get(i)[tailnode].getPosition() +readClass.getHeadAndTail().get(i)[tailnode].getWidth() +2 < samRecordBuffer.getUnclippedStart()) {
						try {
							//readClass.setLastRead(addNode);	
							startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((i+1)*(Main.drawCanvas.readHeight+2))));							
							addNode.setRectBounds(xpos, startY, (int)(addNode.getWidth()*pixel), Main.drawCanvas.readHeight);
							addNode.setPrev(readClass.getHeadAndTail().get(i)[tailnode]);						
							readClass.getHeadAndTail().get(i)[tailnode].setNext(addNode);
							readClass.getHeadAndTail().get(i)[tailnode] = addNode;	
							hashindex = readClass.getHelpHash().get(i);
							readClass.getHelpArray().get(hashindex)[1] = addNode.getPosition()+addNode.getWidth()+2;
							
							if(readClass.getHelpArray().size() > 1) {
								
								tempadder = readClass.getHelpArray().get(hashindex);
								
								readClass.getHelpArray().remove(hashindex);	
								index = searchIndex(readClass.getHelpArray(), tempadder[1], 0, readClass.getHelpArray().size());
								
								if(index == readClass.getHelpArray().size() || index == -1) {
									
									readClass.getHelpArray().add(tempadder);
									readClass.getHelpHash().put(i, readClass.getHelpArray().size()-1);
								}
								else {
									
									readClass.getHelpArray().add(index, tempadder);
									readClass.getHelpHash().put(i, index);
								}
							}
							
							break;
						}
						catch(Exception e) {
							ErrorLog.addError(e.getStackTrace());
							e.printStackTrace();
						}
						
					}
				}
				
			/*	startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((readClass.getHelpArray().get(0)[0]+1)*(Main.drawCanvas.readHeight+2))));							
				addNode.setRectBounds(xpos, startY, (int)(addNode.getWidth()*pixel), Main.drawCanvas.readHeight);
				addNode.setPrev(readClass.getHeadAndTail().get(readClass.getHelpArray().get(0)[0])[tailnode]);						
				readClass.getHeadAndTail().get(readClass.getHelpArray().get(0)[0])[tailnode].setNext(addNode);
				readClass.getHeadAndTail().get(readClass.getHelpArray().get(0)[0])[tailnode] = addNode;			
				
				
			}
			*/
			
			mundane = null;
			found = false;			
			
			isDiscordant = MethodLibrary.isDiscordant(samRecordBuffer, readClass.complete);		
			
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
											addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, splitIndex.readSequence, splitIndex.readSeqStart, readClass.sample, splitIndex,readClass, mismatches);	
											readClass.setLastRead(addNode);	
										}									
										
										
										startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-((i+1)*(readClass.readHeight+2))));							
										addNode.setRectBounds(xpos, startY, (int)(addNode.getWidth()*pixel), readClass.readHeight);
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
								addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, splitIndex.readSequence, splitIndex.readSeqStart, readClass.sample, splitIndex,readClass, mismatches);	
								
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
								addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, splitIndex.readSequence, splitIndex.readSeqStart, readClass.sample, splitIndex,readClass, mismatches);	
								
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
						addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, splitIndex.readSequence, splitIndex.readSeqStart, readClass.sample, splitIndex,readClass, mismatches);	
						
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

void removeNonListVariants() {
	
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
						if(!Main.drawCanvas.hideVar(entry.getValue().get(m))) {
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
void removeBedLinks() {
	
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
				node.putPrev(null);
				node.putNext(null);				
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
	Draw.calculateVars = false;
	Draw.variantcalculator = true;
	SplitClass split = Main.drawCanvas.splits.get(0);
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
	//boolean coding = false;
	cancelvarcount = false;
	VariantHandler.table.transarray.clear();
	VariantHandler.table.aminoarray.clear();	
	VariantHandler.table.controlarray = null;
	VariantHandler.stattable.sampleArray.clear();


	for(int i = 0; i<VariantHandler.tables.size(); i++) {
		VariantHandler.tables.get(i).bedarray.clear();
		VariantHandler.tables.get(i).aminoarray.clear();
		VariantHandler.tables.get(i).vararray.clear();
		VariantHandler.tables.get(i).controlarray = null;
		VariantHandler.tables.get(i).variants = 0;
		
	}	
	
	VarNode vardraw = Main.drawCanvas.current;
	
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
		
	}
	VariantHandler.stattable.setPreferredSize(new Dimension(VariantHandler.statsScroll.getViewport().getWidth(), (VariantHandler.stattable.sampleArray.size()+1)*15+2));
	VariantHandler.stattable.revalidate();
	VariantHandler.stattable.repaint();
	
	for(int g = 0 ; g<split.getGenes().size(); g++) {
		for(int i = 0; i<split.getGenes().get(g).getTranscripts().size(); i++) {
			split.getGenes().get(g).getTranscripts().get(i).mutations = 0;
			split.getGenes().get(g).getTranscripts().get(i).missense = 0;
			split.getGenes().get(g).getTranscripts().get(i).nonsense = 0;
			split.getGenes().get(g).getTranscripts().get(i).synonymous = 0;
			split.getGenes().get(g).getTranscripts().get(i).intronic = 0;
			split.getGenes().get(g).getTranscripts().get(i).utr = 0;
			split.getGenes().get(g).getTranscripts().get(i).samples.clear();
			split.getGenes().get(g).getTranscripts().get(i).varnodes.clear();
		}
	}
	
	if(VariantHandler.allChroms.isSelected()) {		
		if(Main.drawCanvas.splits.get(0).start != 1 || Main.chromosomeDropdown.getSelectedIndex() != 0) {
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
	while(true) {
		
		if(cancelvarcount) {
			cancelFileRead();
			vardraw = null;
  			break;
  		}
		
		if(vardraw == null) {
			
			if(VariantHandler.writetofile.isSelected()) {
				if(VariantHandler.table.transarray.size() > 0) {
					for(int i = 0 ; i<VariantHandler.table.transarray.size(); i++) {
						
							VariantHandler.writeTranscriptToFile(VariantHandler.table.transarray.get(i), output);
							for(int s = 0; s<VariantHandler.table.transarray.get(i).varnodes.size(); s++) {
								if(VariantHandler.table.transarray.get(i).varnodes.get(s).getExons() != null) {
									VariantHandler.table.transarray.get(i).varnodes.get(s).getExons().clear();
								}
								if(VariantHandler.table.transarray.get(i).varnodes.get(s).getTranscripts() != null) {
									VariantHandler.table.transarray.get(i).varnodes.get(s).getTranscripts().clear();
								}
							}							
							VariantHandler.table.transarray.get(i).samples.clear();
							VariantHandler.table.transarray.get(i).varnodes.clear();
							VariantHandler.table.transarray.remove(i);
							i--;						
					}
				}
				
			}
			try {
				if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {					
		//			Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).width = Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.get(Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.size()-1).getPosition() -Main.drawCanvas.clusterNodes.get(Main.drawCanvas.clusterNodes.size()-1).varnodes.get(0).getPosition() +1;
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
			if(VariantHandler.allChroms.isSelected()) {
							
				if(cancelvarcount) {
					cancelFileRead();
		  			break;
		  		}
				if(Main.selectedChrom < Main.chromosomeDropdown.getItemCount()) {
					Main.nothread = true;					
					Main.chromosomeDropdown.setSelectedIndex(Main.selectedChrom+1);
					Main.nothread = false;
					
					if(Main.drawCanvas.splits.get(0).getGenes().size() == 0) {
						break;
					}
					
					for(int g = 0 ; g<split.getGenes().size(); g++) {
						for(int i = 0; i<split.getGenes().get(g).getTranscripts().size(); i++) {
							split.getGenes().get(g).getTranscripts().get(i).mutations = 0;
							split.getGenes().get(g).getTranscripts().get(i).missense = 0;
							split.getGenes().get(g).getTranscripts().get(i).nonsense = 0;
							split.getGenes().get(g).getTranscripts().get(i).synonymous = 0;
							split.getGenes().get(g).getTranscripts().get(i).intronic = 0;
							split.getGenes().get(g).getTranscripts().get(i).utr = 0;
							split.getGenes().get(g).getTranscripts().get(i).samples.clear();
							split.getGenes().get(g).getTranscripts().get(i).varnodes.clear();
						}
					}	
					Main.drawCanvas.calcClusters(FileRead.head.getNext());
					vardraw = FileRead.head.getNext();		
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
			}
			if(vardraw == null || vardraw.getPosition() > Main.drawCanvas.splits.get(0).end) {
				vardraw= null;
				continue;
			}
		}		
		//STATS
		if(vardraw == null) {
			continue;
		}
	//	Main.drawCanvas.loadbarAll = (int)((vardraw.getPosition()/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);
		if(Main.drawCanvas.hideNode(vardraw)) {			
			vardraw = vardraw.getNext();			
			continue;			
		}		
	
		annotateVariant(vardraw);		
		vardraw = vardraw.getNext();
		
		}
		catch(Exception ex) {
			ErrorLog.addError(ex.getStackTrace());
			ex.printStackTrace();
		}
	}
	
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
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	split = null;
	VariantHandler.table.revalidate();
	VariantHandler.table.repaint();
	VariantHandler.freeze.setSelected(false);
	Main.drawCanvas.current = FileRead.head;
	if(!Main.drawCanvas.loading) {
		Draw.calculateVars = true;
	}
	Draw.variantcalculator = false;
}

void varCalcBig() {
	bigcalc = false;
	Draw.calculateVars = false;
	Draw.variantcalculator = true;
	int adder = 2000000, presearchpos = 1;
	
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
	VariantHandler.table.transarray.clear();
	VariantHandler.table.aminoarray.clear();	
	VariantHandler.table.controlarray = null;
	VariantHandler.stattable.sampleArray.clear();


	for(int i = 0; i<VariantHandler.tables.size(); i++) {
		VariantHandler.tables.get(i).bedarray.clear();
		VariantHandler.tables.get(i).aminoarray.clear();
		VariantHandler.tables.get(i).vararray.clear();
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
		
	}
	VariantHandler.stattable.setPreferredSize(new Dimension(VariantHandler.statsScroll.getViewport().getWidth(), (VariantHandler.stattable.sampleArray.size()+1)*15+2));
	VariantHandler.stattable.revalidate();
	VariantHandler.stattable.repaint();
	SplitClass split = Main.drawCanvas.splits.get(0);
	for(int g = 0 ; g<split.getGenes().size(); g++) {
		for(int i = 0; i<split.getGenes().get(g).getTranscripts().size(); i++) {
			split.getGenes().get(g).getTranscripts().get(i).mutations = 0;
			split.getGenes().get(g).getTranscripts().get(i).missense = 0;
			split.getGenes().get(g).getTranscripts().get(i).nonsense = 0;
			split.getGenes().get(g).getTranscripts().get(i).synonymous = 0;
			split.getGenes().get(g).getTranscripts().get(i).intronic = 0;
			split.getGenes().get(g).getTranscripts().get(i).utr = 0;
			split.getGenes().get(g).getTranscripts().get(i).samples.clear();
			split.getGenes().get(g).getTranscripts().get(i).varnodes.clear();
		}
	}
	split = null;
	if(VariantHandler.allChroms.isSelected()) {		
		if(Main.drawCanvas.splits.get(0).start != 1 || Main.chromosomeDropdown.getSelectedIndex() != 0) {
			Main.nothread = true;
			search = true;
			
			Main.chromosomeDropdown.setSelectedIndex(0);
			Main.nothread = false;
			vardraw = FileRead.head.getNext();		
			Main.drawCanvas.clusterNodes.clear();
			Main.drawCanvas.calcClusters(FileRead.head.getNext());
			
		}
		else {
			Main.nothread = true;
			search = true;
			searchStart = (int)Main.drawCanvas.splits.get(0).start;
			searchEnd = (int)Main.drawCanvas.splits.get(0).start+1000000;
			Main.chromosomeDropdown.setSelectedIndex(0);
			Main.nothread = false;
			
			Main.drawCanvas.clusterNodes.clear();
			Main.drawCanvas.calcClusters(FileRead.head.getNext());
			vardraw = FileRead.head.getNext();
		}
		
	}	
	else {
		//Main.nothread = true;
		//search = true;
	
		searchStart = 1;
		
		searchEnd = getNextIntergenic(presearchpos+adder);
		presearchpos = searchEnd;
		
		getVariantWindow(Main.drawCanvas.splits.get(0).chrom,searchStart,searchEnd);
		vardraw = FileRead.head.getNext();
		
	}
	bigcalc = true;
	Main.drawCanvas.loadbarAll = 0;
	while(true) {
		
		if(cancelvarcount || cancelfileread) {
			cancelFileRead();
			vardraw = null;
  			break;
  		}
		if(vardraw == null) {
			if(VariantHandler.writetofile.isSelected()) {
				if(VariantHandler.table.transarray.size() > 0) {
					for(int i = 0 ; i<VariantHandler.table.transarray.size(); i++) {
						
							VariantHandler.writeTranscriptToFile(VariantHandler.table.transarray.get(i), output);
							for(int s = 0; s<VariantHandler.table.transarray.get(i).varnodes.size(); s++) {
								if(VariantHandler.table.transarray.get(i).varnodes.get(s).getExons() != null) {
									VariantHandler.table.transarray.get(i).varnodes.get(s).getExons().clear();
								}
								if(VariantHandler.table.transarray.get(i).varnodes.get(s).getTranscripts() != null) {
									VariantHandler.table.transarray.get(i).varnodes.get(s).getTranscripts().clear();
								}
							}
							
							VariantHandler.table.transarray.get(i).samples.clear();
							VariantHandler.table.transarray.get(i).varnodes.clear();
							VariantHandler.table.transarray.remove(i);
							i--;						
					}
				}				
			}
			try {
				if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {					
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
				presearchpos = 1;
				if(cancelvarcount) {
		  			break;
		  		}
				if(Main.selectedChrom < Main.chromosomeDropdown.getItemCount()) {
					Main.nothread = true;					
					search = true;
					searchStart = 1;
					searchEnd = adder;
					Main.chromosomeDropdown.setSelectedIndex(Main.selectedChrom+1);
					search = false;
					Main.nothread = false;
					
					if(Main.drawCanvas.splits.get(0).getGenes().size() == 0) {
						break;
					}
					for(int g = 0 ; g<Main.drawCanvas.splits.get(0).getGenes().size(); g++) {
						for(int i = 0; i<Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().size(); i++) {
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).mutations = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).missense = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).nonsense = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).synonymous = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).intronic = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).utr = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).samples.clear();
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).varnodes.clear();
						}
					}	
					Main.drawCanvas.calcClusters(FileRead.head.getNext());
					vardraw = head.getNext();		
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
					break;
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
				if(searchEnd > Main.drawCanvas.splits.get(0).end) {
					break;
				}
				searchStart = presearchpos;
				searchEnd = getNextIntergenic(presearchpos+adder);
				presearchpos = searchEnd;
				
				getVariantWindow(Main.drawCanvas.splits.get(0).chrom,searchStart,searchEnd);
			/*	for(int i = 0; i<Main.samples; i++) {
					if( Main.drawCanvas.sampleList.get(i).getTabixFile() == null ||  Main.drawCanvas.sampleList.get(i).multipart) {
						continue;
					}					
					getVariants(Main.drawCanvas.splits.get(0).chrom,searchStart,searchEnd, Main.drawCanvas.sampleList.get(i));					
					
				}				
				
				annotate();
				Main.drawCanvas.current = head;		
				*/
				vardraw = head.getNext();
				
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
		if(vardraw == null) {
			continue;
		}
	//	Main.drawCanvas.loadbarAll = (int)((vardraw.getPosition()/(double)(Main.drawCanvas.splits.get(0).chromEnd))*100);
		if(Main.drawCanvas.hideNode(vardraw)) {
			
			vardraw = vardraw.getNext();
			
			continue;			
		}		
		
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
		annotateVariant(vardraw);
	
	
		
		vardraw = vardraw.getNext();
		
		}
		catch(Exception ex) {
			ErrorLog.addError(ex.getStackTrace());
			ex.printStackTrace();
		}
	}
	
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
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	VariantHandler.table.revalidate();
	VariantHandler.table.repaint();
	VariantHandler.freeze.setSelected(false);
	Main.drawCanvas.current = FileRead.head;
	if(!Main.drawCanvas.loading) {
		Draw.calculateVars = true;
	}
	bigcalc = false;
	Draw.variantcalculator = false;
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
			returnpos = Main.drawCanvas.splits.get(0).getGenes().get(gene).getEnd()+1;
			pointer = gene;
		}
	}
	
	return returnpos;
}

void annotateVariant(VarNode vardraw) {
	if(VariantHandler.writetofile.isSelected()) {
		if(VariantHandler.table.transarray.size() > 0) {
			for(int i = 0 ; i<VariantHandler.table.transarray.size(); i++) {
				if(vardraw.getPosition() > VariantHandler.table.transarray.get(i).getEnd()) {
					
					VariantHandler.writeTranscriptToFile(VariantHandler.table.transarray.get(i), output);
					if(VariantHandler.allChroms.isSelected() || bigcalc) {
						for(int s = 0; s<VariantHandler.table.transarray.get(i).varnodes.size(); s++) {
							if(VariantHandler.table.transarray.get(i).varnodes.get(s).getExons() != null) {
								for(int e = 0 ; e<VariantHandler.table.transarray.get(i).varnodes.get(s).getExons().size(); e++) {
									if(VariantHandler.table.transarray.get(i).varnodes.get(s).getExons().get(e).getTranscript().equals(VariantHandler.table.transarray.get(i))) {
										VariantHandler.table.transarray.get(i).varnodes.get(s).getExons().remove(e);
										e--;
									}
								}						
							}
							if(VariantHandler.table.transarray.get(i).varnodes.get(s).getTranscripts() != null) {
								for(int e = 0 ; e<VariantHandler.table.transarray.get(i).varnodes.get(s).getTranscripts().size(); e++) {
									if(VariantHandler.table.transarray.get(i).varnodes.get(s).getTranscripts().equals(VariantHandler.table.transarray.get(i))) {
										VariantHandler.table.transarray.get(i).varnodes.get(s).getTranscripts().remove(e);
										e--;
									}
									
								}
							}
						}
					}
					VariantHandler.table.transarray.get(i).samples.clear();
					VariantHandler.table.transarray.get(i).varnodes.clear();
					VariantHandler.table.transarray.remove(i);
					i--;
				}
			}
		}
		
	}
	Map.Entry<String, ArrayList<SampleNode>> entry = null;
	String base = null, amino = null;
	int pretrack = -1;
	mutcount = 0;
	for(int v = 0; v<vardraw.vars.size(); v++) {
	 	entry = vardraw.vars.get(v);
		if(Main.drawCanvas.hideNodeVar(vardraw, entry)) {
			continue;
		}			
		
		base = entry.getKey();
		for(int m = 0; m<entry.getValue().size(); m++) {
			if(!Main.drawCanvas.hideVar(entry.getValue().get(m))) {
				sample = entry.getValue().get(m).getSample();
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
					sample.snvs++;
					try {
						sample.mutationTypes[Main.mutTypes.get(vardraw.getRefBase() +entry.getKey())]++;
					}
					catch(Exception e) {
						System.out.println(vardraw.getPosition() +" " +vardraw.getRefBase() +entry.getKey());
						e.printStackTrace();
					}
					if((vardraw.getRefBase().equals("A") && entry.getKey().equals("G")) || (vardraw.getRefBase().equals("G") && entry.getKey().equals("A")) || (vardraw.getRefBase().equals("C") && entry.getKey().equals("T")) || (vardraw.getRefBase().equals("T") && entry.getKey().equals("C"))) {
						sample.sitioRate++;
					}
					else {
						sample.versioRate++;
					}						
				}										
			}
		}
	}
	
	if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {
		if(Main.drawCanvas.clusterNodes.size() == 0) {
			ClusterNode cluster = new ClusterNode();
			cluster.nodecount= mutcount;
			cluster.ID = vardraw.clusterId;
			cluster.varnodes.add(vardraw);			
			cluster.width = 1;
			Main.drawCanvas.clusterNodes.add(cluster);
		}
		else if(vardraw.clusterId-1 >= Main.drawCanvas.clusterNodes.size()) {
			Main.drawCanvas.clusterNodes.get(vardraw.clusterId-2).width = Main.drawCanvas.clusterNodes.get(vardraw.clusterId-2).varnodes.get(Main.drawCanvas.clusterNodes.get(vardraw.clusterId-2).varnodes.size()-1).getPosition() -Main.drawCanvas.clusterNodes.get(vardraw.clusterId-2).varnodes.get(0).getPosition() +1;
			
			ClusterNode cluster = new ClusterNode();
			cluster.nodecount= mutcount;
			cluster.ID = vardraw.clusterId;
			cluster.varnodes.add(vardraw);			
			cluster.width = 1;
			Main.drawCanvas.clusterNodes.add(cluster);
		}
		else {
			Main.drawCanvas.clusterNodes.get(vardraw.clusterId-1).nodecount+=mutcount;
			Main.drawCanvas.clusterNodes.get(vardraw.clusterId-1).varnodes.add(vardraw);
		}
		VariantHandler.clusterTable.variants+=mutcount;
	}
	
	if(vardraw != null && vardraw.bedhit) {
		
		for(int v = 0; v<vardraw.vars.size(); v++) {
			mutcount = 0;
		 	entry = vardraw.vars.get(v);
			if(Main.drawCanvas.hideNodeVar(vardraw, entry)) {
				continue;
			}
							
			for(int m = 0; m<entry.getValue().size(); m++) {
				if(!Main.drawCanvas.hideVar(entry.getValue().get(m))) {
					if(VariantHandler.onlyselected.isSelected()) {
						if(!entry.getValue().get(m).getSample().equals(Main.drawCanvas.selectedSample)) {
							continue;
						}
					}
					mutcount++;
				}
			}
			if(mutcount == 0) {
				continue;
			}
			pretrack = -1;
			for(int i = 0; i<vardraw.getBedHits().size(); i++) {
				if(pretrack != vardraw.getBedHits().get(i).getTrack().trackIndex) {
					vardraw.getBedHits().get(i).getTrack().getTable().variants += mutcount;
					pretrack = vardraw.getBedHits().get(i).getTrack().trackIndex;
				}
				if(!vardraw.getBedHits().get(i).getTrack().getTable().bedarray.contains(vardraw.getBedHits().get(i))) {
					vardraw.getBedHits().get(i).mutations= mutcount;
					vardraw.getBedHits().get(i).getTrack().getTable().bedarray.add(vardraw.getBedHits().get(i));
					vardraw.getBedHits().get(i).getTrack().getTable().setPreferredSize(new Dimension(vardraw.getBedHits().get(i).getTrack().getTable().tablescroll.getViewport().getWidth(), (vardraw.getBedHits().get(i).getTrack().getTable().getTableSize()+2+samplecount)*vardraw.getBedHits().get(i).getTrack().getTable().rowHeight));										
					vardraw.getBedHits().get(i).getTrack().getTable().revalidate();		
					vardraw.inVarList = true;
				}
				else {
					vardraw.getBedHits().get(i).mutations+= mutcount;
					vardraw.inVarList = true;
				}
			}
		}
		
	}
	// INTRONIC
	
	if(VariantHandler.intronic.isSelected() && vardraw.isInGene() && vardraw.getTranscripts() != null ) {
		
		for(int t = 0; t<vardraw.getTranscripts().size(); t++) {
			if(!VariantHandler.allIsoforms.isSelected() && vardraw.getTranscripts().get(t).getGene().getCanonical() != null && !vardraw.getTranscripts().get(t).isCanonical()) {
				continue;
			}
			
			for(int v = 0; v<vardraw.vars.size(); v++) {
				 	entry = vardraw.vars.get(v);
					if(Main.drawCanvas.hideNodeVar(vardraw, entry)) {
						continue;
					}
					
					mutcount = 0;
					base = entry.getKey();
					for(int m = 0; m<entry.getValue().size(); m++) {
						if(!Main.drawCanvas.hideVar(entry.getValue().get(m))) {
							if(VariantHandler.onlyselected.isSelected()) {
								if(!entry.getValue().get(m).getSample().equals(Main.drawCanvas.selectedSample)) {
									continue;
								}
							}
							mutcount++;
						}
					}
					if(mutcount == 0) {
						continue;
					}
				
					for(int i = 0; i<entry.getValue().size(); i++) {
						if(entry.getValue().get(i).alleles != null) {
							break;
						}
						if(Main.drawCanvas.hideVar(entry.getValue().get(i))) {
							continue;
						}
						if(VariantHandler.onlyselected.isSelected()) {
							if(!entry.getValue().get(i).getSample().equals(Main.drawCanvas.selectedSample)) {
								continue;
							}
						}
						if(!vardraw.getTranscripts().get(t).samples.contains(entry.getValue().get(i).getSample())) {
							vardraw.getTranscripts().get(t).samples.add(entry.getValue().get(i).getSample());
						}
					}
					if(vardraw.getTranscripts().get(t).mutations == 0 && vardraw.getTranscripts().get(t).samples.size() >= VariantHandler.geneSlider.getValue()) {
						VariantHandler.table.addEntry(vardraw.getTranscripts().get(t));
						
					}
					
					vardraw.inVarList = true;
					
					vardraw.getTranscripts().get(t).mutations+=mutcount;
					vardraw.getTranscripts().get(t).intronic+=mutcount;						
					VariantHandler.table.variants += mutcount;
					VariantHandler.table.revalidate();
					VariantHandler.table.repaint();						
				}	
			if(vardraw.inVarList) {
				vardraw.getTranscripts().get(t).varnodes.add(vardraw);
			}
		}	
	}
	
	if(vardraw != null && vardraw.getExons() != null) {	
		
		for(int exon = 0; exon<vardraw.getExons().size(); exon++) {
			if(!VariantHandler.allIsoforms.isSelected() && vardraw.getExons().get(exon).transcript.getGene().getCanonical() != null && !vardraw.getExons().get(exon).transcript.isCanonical()) {
				continue;
			}
			foundexon = -1;
			
			for(int v = 0; v<vardraw.vars.size(); v++) {
				 entry = vardraw.vars.get(v);
				if(Main.drawCanvas.hideNodeVar(vardraw, entry)) {
					continue;
				}			
			
				mutcount = 0;
				
				base = entry.getKey();
				
				for(int m = 0; m<entry.getValue().size(); m++) {
					
					if(!Main.drawCanvas.hideVar(entry.getValue().get(m))) {
						
						if(VariantHandler.onlyselected.isSelected()) {
							if(!entry.getValue().get(m).getSample().equals(Main.drawCanvas.selectedSample)) {
								continue;
							}
						}
						mutcount++;
					}
				}
				
				if(mutcount == 0) {
					continue;
				}
				
				
				if(!vardraw.coding) {
					
					if(VariantHandler.utr.isSelected()) {
						for(int i = 0; i<entry.getValue().size(); i++) {
							if(entry.getValue().get(i).alleles != null) {
								break;
							}
							if(Main.drawCanvas.hideVar(entry.getValue().get(i))) {
								continue;
							}
							if(VariantHandler.onlyselected.isSelected()) {
								if(!entry.getValue().get(i).getSample().equals(Main.drawCanvas.selectedSample)) {
									continue;
								}
							}
							if(!vardraw.getExons().get(exon).getTranscript().samples.contains(entry.getValue().get(i).getSample())) {
								vardraw.getExons().get(exon).getTranscript().samples.add(entry.getValue().get(i).getSample());
							}
						}
						
						if(!VariantHandler.table.transarray.contains(vardraw.getExons().get(exon).getTranscript()) && vardraw.getExons().get(exon).getTranscript().samples.size() >= VariantHandler.geneSlider.getValue()) {
							
							VariantHandler.table.addEntry(vardraw.getExons().get(exon).getTranscript());
						}							
						vardraw.inVarList = true;
						foundexon = exon;
						
						vardraw.getExons().get(exon).getTranscript().mutations+=mutcount;
						vardraw.getExons().get(exon).getTranscript().utr +=mutcount;
						VariantHandler.table.variants +=mutcount;
						continue;
					}
					else {
						continue;
					}						
				}				
				amino = Main.chromDraw.getAminoChange(vardraw,base,vardraw.getExons().get(exon));
				
				if(amino.contains("UTR")) {
					if(!VariantHandler.table.transarray.contains(vardraw.getExons().get(exon).getTranscript()) && vardraw.getExons().get(exon).getTranscript().samples.size() >= VariantHandler.geneSlider.getValue()) {
						
						VariantHandler.table.addEntry(vardraw.getExons().get(exon).getTranscript());
					}							
					vardraw.inVarList = true;
					foundexon = exon;
					
					vardraw.getExons().get(exon).getTranscript().mutations+=mutcount;
					vardraw.getExons().get(exon).getTranscript().utr +=mutcount;
					VariantHandler.table.variants +=mutcount;
					continue;
				}
				if(amino.equals("")) {
					continue;
				}
				if(VariantHandler.nonsense.isSelected()) {
					if(!MethodLibrary.aminoEffect(amino).equals("nonsense")) {							
						continue;
					}						
				}
				else if(VariantHandler.synonymous.isSelected()) {						
					if(MethodLibrary.aminoEffect(amino).equals("synonymous")) {							
						continue;
					}							
				}					
				
				foundexon = exon;
				if(amino.length() > 3) {
					for(int i = 0; i<entry.getValue().size(); i++) {
						if(entry.getValue().get(i).alleles != null) {
							break;
						}
						if(Main.drawCanvas.hideVar(entry.getValue().get(i))) {
							continue;
						}
						if(VariantHandler.onlyselected.isSelected()) {
							if(!entry.getValue().get(i).getSample().equals(Main.drawCanvas.selectedSample)) {
								continue;
							}
						}
						if(!vardraw.getExons().get(exon).getTranscript().samples.contains(entry.getValue().get(i).getSample())) {
							vardraw.getExons().get(exon).getTranscript().samples.add(entry.getValue().get(i).getSample());
							
						}
					}
					
					if(!VariantHandler.table.transarray.contains(vardraw.getExons().get(exon).getTranscript()) && vardraw.getExons().get(exon).getTranscript().samples.size() >= VariantHandler.geneSlider.getValue() ) {
						VariantHandler.table.addEntry(vardraw.getExons().get(exon).getTranscript());								
					}
				
					if(MethodLibrary.aminoEffect(amino).equals("nonsense")) {
						vardraw.getExons().get(exon).getTranscript().nonsense += mutcount;
					}
					else if(MethodLibrary.aminoEffect(amino).equals("synonymous")) {
						
						vardraw.getExons().get(exon).getTranscript().synonymous+= mutcount;
					}
					else if(MethodLibrary.aminoEffect(amino).equals("missense")) {
						vardraw.getExons().get(exon).getTranscript().missense += mutcount;
					}
					
				}					
				vardraw.inVarList = true;
				vardraw.getExons().get(exon).getTranscript().mutations+=mutcount;	
				VariantHandler.table.variants += mutcount;
			}
			
			if(foundexon > -1) {	
							
				vardraw.getExons().get(foundexon).getTranscript().varnodes.add(vardraw);
				VariantHandler.table.revalidate();
				VariantHandler.table.repaint();
			
			}
		}			
	}
}
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
	
	
	
	/*
	ResultSet rs = stmt.executeQuery(
	"SELECT chrom, txStart, txEnd, ensemblToGeneName.value, exonCount, strand, name2, ensGene.name, "
  + "ensemblSource.source, cdsStart, cdsEnd, exonStarts, exonEnds, exonFrames FROM ensGene, ensemblToGeneName, "
  + "ensemblSource WHERE exonCount > 100 AND ensGene.name = ensemblToGeneName.name AND ensemblSource.name = ensemblToGeneName.name"
  + ""
);
	*/
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
	
	
	
/*	
	ResultSet rs = stmt.executeQuery(
	"SELECT chrom, txStart, txEnd, ensemblToGeneName.value, exonCount, strand, name2, ensGene.name, "
  + "ensemblSource.source, cdsStart, cdsEnd, exonStarts, exonEnds, exonFrames FROM ensGene, ensemblToGeneName, "
  + "ensemblSource WHERE exonCount > 100 AND ensGene.name = ensemblToGeneName.name AND ensemblSource.name = ensemblToGeneName.name"
  + ""
);*/
	
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
	/*	ResultSet rs = conn.getMetaData().getTables(null, null, "%", null);
		
		
		while(rs.next()) {
			System.out.println(rs.getString(3));
		}
		*/
	//	ResultSet rs = stmt.executeQuery("select S.*,X.*,G.* from ensemblSource as S,ensGene as G,knownToEnsembl as KE, kgXref as X where S.name = G.name and G.name=KE.value and KE.name=X.kgID");
		System.out.println("Fetching...");
		ResultSet rs = stmt.executeQuery("select E.chrom, E.txStart, E.txEnd, N.value, E.exonCount, E.strand, E.name2, E.name, E.name, E.name, S.source, E.cdsStart, E.cdsEnd, E.exonStarts, E.exonEnds, E.exonFrames, E.name from "
				//+ "knownToEnsembl as K left outer join knownCanonical as C on K.name = C.transcript "
				+ "ensGene as E left outer join ensemblToGeneName as N on E.name = N.name "
				+ "left outer join ensemblSource as S on S.name = E.name");				
				
		System.out.println("Fecthed");
	/*	ResultSet rs = stmt.executeQuery(
		"SELECT chrom, txStart, txEnd, ensemblToGeneName.value, exonCount, strand, name2, ensGene.name, "
	  + "ensemblSource.source, cdsStart, cdsEnd, exonStarts, exonEnds, exonFrames FROM ensGene, ensemblToGeneName, "
	  + "ensemblSource WHERE exonCount > 100 AND ensGene.name = ensemblToGeneName.name AND ensemblSource.name = ensemblToGeneName.name"
	  + ""
	);
		*/
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
}
	/*
	public static SortedSet<Segment> getReads2(int chrom, int startpos, int endpos, Sample sample) {
				
		if(sample.samFile != null) {
			SortedSet<Segment> data = Collections.synchronizedSortedSet(new TreeSet<Segment>());
			boolean split = false;
			int level = 1;
			int count = 0;			
						
			String readString;
		//	if(Runtime.getRuntime().freeMemory() < 50000000) {
			//	System.out.println("Low memory. Free: " +(Runtime.getRuntime().freeMemory()/1000000) +"MB");
			/*	
				for(int i = 0; i< Draw.visibleStart; i++) {
					Main.drawVector.get(i).data.clear();
				}
				for(int i = Draw.drawVariables.visibleend; i< Main.drawVector.size(); i++) {
					Main.drawVector.get(i).data.clear();
				}
		//		System.gc();
		//	}	
			
		try {	
			    
				
					
				boolean chr = false;
				if(sample.samFileReader != null) {
					sample.samFileReader.close();
					
				}
			//	boolean first = false;
				
				
				sample.samFileReader = new SAMFileReader(sample.samFile);
				
				sample.samFileReader.setValidationStringency(ValidationStringency.SILENT);
				
				int chromTemp = chrom;
			//	System.out.println(chrom);
				if(chrom == 24) {
					chromTemp = sample.samFileReader.getFileHeader().getSequenceIndex("MT");
					if(chromTemp == -1) {
						chromTemp = sample.samFileReader.getFileHeader().getSequenceIndex("M");
					}				
				}
				else if(sample.samFileReader.getFileHeader().getSequenceIndex("" +Main.chromosomeDropdown.getItemAt(chrom)) == -1) {
					chromTemp = sample.samFileReader.getFileHeader().getSequenceIndex("chr" +Main.chromosomeDropdown.getItemAt(chrom));
					chr = true;				
				}
				else {
					chromTemp = sample.samFileReader.getFileHeader().getSequenceIndex("" +Main.chromosomeDropdown.getItemAt(chrom));
				}
				int calc = 0;
				
				
					
				try {
				
						sample.bamIterator = sample.samFileReader.iterator(sample.samFileReader.getIndex().getSpanOverlapping(chromTemp, startpos, endpos));					
						
				
				}
				catch(Exception ex) {
					
					ex.printStackTrace();
							
				}
			SortedSet<Level> levels = Collections.synchronizedSortedSet(new TreeSet<Level>());
			boolean dir = false;
			CigarElement cigar;
			int start = 0, end = 0;
			Level templevel = null;
			double pixel = pixel, startX; 	
		
			int endIndex = 0;
			String str;
			int length = 0, Nlen, j = 0;
			int distance = 0;
			int i = 0, s =0, startY;
			int len = 0,insertsize=0;
			Iterator<Level> levelIterator;
			Level lev,leveladd;
			
			
			if(sample.bamIterator != null && sample.bamIterator.hasNext()) {			
				
				length = 0;
				distance = 0;
				endIndex = 0;
				
				
				Segment seg= null;	
				
				if(pixel > 1) {
					for(i = 0; i<draw.covdata.length; i++) {
						draw.covdata[i][0] = Draw.start +i;
					}
				}
				else {
					draw.covdata[0][0] = Draw.start;
					for(i = 1; i<draw.covdata.length; i++) {
						draw.covdata[i][0] = draw.covdata[i-1][0] +(1/pixel);					
					}
					
				}
		
			while(sample.bamIterator.hasNext() && sample.bamIterator != null) {
				split = false;				
				sample.samRecord =sample.bamIterator.next(); 				
				
				if(sample.samRecord.getAlignmentStart() < start-10000) {
					continue;
				}
				if(sample.samRecord.getAlignmentEnd() > Main.drawCanvas.end+10000){
					break;
				}
				if(sample.samRecord.getReadUnmappedFlag()) {		
					
					continue;
				}			
			if(Main.clipped.isSelected())  {
				try{
					end = draw.samRecord.getUnclippedEnd();
				}
				catch(Exception e) {
					end = draw.samRecord.getAlignmentEnd();
				}
				try {
					start = draw.samRecord.getUnclippedStart();
				}
				catch(Exception e) {
					start = draw.samRecord.getAlignmentStart();
				}
			}
			else {
			
				end = sample.samRecord.getUnclippedEnd()-1;
				start = sample.samRecord.getUnclippedStart()-1;	
				length = (end-start)+1;
			//}
			insertsize = sample.samRecord.getInferredInsertSize();
			if(insertsize == 0) {
				insertsize = (end-start)+1;
			}
			
				if(Main.discordant.isSelected()) {
					
					if(draw.samRecord.getMateUnmappedFlag()) {
						continue;
					}
					if(!draw.complete) {
						if((draw.samRecord.getMateReferenceIndex() != draw.samRecord.getReferenceIndex()) || (Math.abs(insertsize) > Main.insertSlide.getValue()) || draw.samRecord.getReadNegativeStrandFlag() == draw.samRecord.getMateNegativeStrandFlag() || (draw.samRecord.getReadNegativeStrandFlag() && draw.samRecord.getAlignmentStart() < draw.samRecord.getMateAlignmentStart())  || (!draw.samRecord.getReadNegativeStrandFlag() && draw.samRecord.getAlignmentStart() > draw.samRecord.getMateAlignmentStart())) {
							
						}
						else {
							continue;
						}
					}
					else {
						if((draw.samRecord.getMateReferenceIndex() != draw.samRecord.getReferenceIndex()) || (Math.abs(insertsize) > Main.insertSlide.getValue()) || draw.samRecord.getReadNegativeStrandFlag() != draw.samRecord.getMateNegativeStrandFlag()) {
							
						}
						else {
							continue;
						}
					}				
						
				}
				*/
			/*	if(Draw.viewLength > 10000) {
					if(end < Draw.start) {
						continue;
					}
					if(start > Draw.end) {
						break;
					}		
				}
				else {
					if(end < Draw.start-10000) {
						continue;
					}
					if(start > Draw.end+10000) {
						break;
					}		
				}
				
*/
			/*	
					if(draw.samRecord != null && draw.samRecord.getCigarString() != null && draw.samRecord.getCigarString().contains("N")) {	
						
						len = 0;
						cigar = null;
						for(c = 0; c<draw.samRecord.getCigarLength(); c++) {
						
							try {
								cigar = draw.samRecord.getCigar().getCigarElement(c);
							}
							catch(Exception e) {
								break;
							}
							if(cigar.getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH)== 0) {
								len+=cigar.getLength();
							}
							if(cigar.getOperator().compareTo(CigarOperator.N)== 0 || c == draw.samRecord.getCigarLength()-1) {
								index = start-startpos; 
								length = len;
								
								if(pixel > 1) {
									
									for(i = index; i<index+length; i++) {
										if(i < 0) {
											continue;
										}
										if(i > draw.covdata.length-1) {
											break;
										}							
											draw.covdata[i][1]++;							
									}				
								}
								else {
									
									index = (int)(((start)-startpos)*pixel);
									endIndex = (int)(((start+length)-startpos)*pixel);
									
									if(endIndex != index) {
										
										for(i = index; i<=endIndex; i++) {
											if(i < 0) {
												continue;
											}
											if(endIndex > draw.covdata.length-1) {
												break;
											}
									//		draw.covdata[i][0] = start+i;
											if(i == index) {
												draw.covdata[i][1] += ((startpos+((1/pixel))*(i+1))-start)/(double)(1/pixel); //(double)((endIndex+1)-index);
											}
											else if(i != endIndex) {
												draw.covdata[i][1] += 1;
											}
											else {
												draw.covdata[i][1] += ((start+length)-(startpos+(1/pixel)*(endIndex)))/(double)(1/pixel); //(double)((endIndex+1)-index);
											}
										}
									}
									else {
										if(index < draw.covdata.length && index >= 0) {
											draw.covdata[index][1] += length/(double)(1/pixel);		
										}
									}							
							}
								start = start+len+cigar.getLength();
								len = 0;
						}
						}
					}
					else {

						
					index = start-startpos;
					length = (end-start)+1;
					
					if(pixel > 1) {
					
						for(i = index; i<index+length; i++) {
							if(i < 0) {
								continue;
							}
							if(i > draw.covdata.length-1) {
								break;
							}							
								draw.covdata[i][1]++;							
						}				
					}
					else {
						
						index = (int)(((start)-startpos)*pixel);
						endIndex = (int)(((start+length)-startpos)*pixel);
						
						if(endIndex != index) {
							
							for(i = index; i<=endIndex; i++) {
								
								if(i < 0) {
									continue;
								}
							/*	if(endIndex > draw.covdata.length-1) {
									break;
								}
								*/
							//	if(i > draw.covdata.length-1) {
							//		break;
							//	}
								
						//		draw.covdata[i][0] = start+i;
						/*		if(i == index) {
									draw.covdata[i][1] += ((startpos+((1/pixel))*(i+1))-start)/(double)(1/pixel); //(double)((endIndex+1)-index);
								}
								else if(i != endIndex) {
									draw.covdata[i][1] += 1;
								}
								else {
									draw.covdata[i][1] += ((start+length)-(startpos+(1/pixel)*(endIndex)))/(double)(1/pixel); //(double)((endIndex+1)-index);
								}
							}
						}
						else {
							if(index < draw.covdata.length && index >= 0) {
								draw.covdata[index][1] += length/(double)(1/pixel);		
							}
						}	
					}
					}
					
						if(Main.drawCanvas.viewLength < 30000 ) {//&& Draw.windowHeight/Draw.scale > 4) {
							
							dir = sample.samRecord.getReadNegativeStrandFlag();	
													
							
							distance = end;
						//	length = draw.samRecord.getReadLength();
							
							if(sample.samRecord.getCigar().numCigarElements() > 1) {
								
								count = 0;
								
								readString = sample.samRecord.getReadString(); //.substring((start - sample.samRecord.getUnclippedStart()),sample.samRecord.getReadLength()-(sample.samRecord.getUnclippedEnd()-end));
								
								
								for(i = 0; i<sample.samRecord.getCigar().numCigarElements(); i++) {
									try {
										cigar = sample.samRecord.getCigar().getCigarElement(i);
									}
									catch(Exception e) {
										e.printStackTrace();
										break;
									}
									if(cigar.getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH) == 0) {									
										count += cigar.getLength();									
									}	 
									
									else if(cigar.getOperator().compareTo(CigarOperator.DELETION) == 0) {	
										str = "";
							//			length += cigar.getLength();
										for(j =0 ; j<cigar.getLength();j++) {
											str += "D";
										}			
										
									//		readString = new StringBuffer(readString).delete(count, count+(cigar.getLength())).toString();
											readString = new StringBuffer(readString).insert(count, str).toString();		
											count+=cigar.getLength();
										
									}
									else if(cigar.getOperator().compareTo(CigarOperator.INSERTION) == 0) {		
								//		length-=cigar.getLength();
											
											readString = new StringBuffer(readString).delete(count, count+(cigar.getLength())).toString();
											readString = new StringBuffer(readString).replace(count, count+1, "I").toString();		
										
										
									}								
										
									else if(cigar.getOperator().compareTo(CigarOperator.N) == 0) {										
																	
											Nlen = 0;
											for(j = i-1; j>= 0; j--) {
												cigar = sample.samRecord.getCigar().getCigarElement(j);
												if(cigar.getOperator().compareTo(CigarOperator.MATCH_OR_MISMATCH) == 0) {
													Nlen += cigar.getLength();
												}	 
												if(cigar.getOperator().compareTo(CigarOperator.N) == 0) {											
													break;
												}										
											}
											cigar = sample.samRecord.getCigar().getCigarElement(i);
											if(split) {
												if(sample.samRecord.getReadPairedFlag() && sample.samRecord.getMateUnmappedFlag()) {
													
													seg.splits.add(new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, Nlen, 0, distance, readString.substring(0, Nlen),"", -1, false, sample.samRecord.getCigar(), sample.samRecord.getMappingQuality()));
									//				seg.splits.add(new Segment(Main.selectedChrom+1, start, draw.samRecord.getReadName(), dir, true, Nlen, 0, distance, readString.substring(0, Nlen),"", -1, false, draw.samRecord.getCigarString()));
													
													//	draw.data.add(seg);
												}
												else if (sample.samRecord.getReadPairedFlag() && !sample.samRecord.getMateUnmappedFlag()){
													
													if(!chr) {
														seg.splits.add(new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, Nlen, insertsize, distance, readString.substring(0, Nlen), sample.samRecord.getMateReferenceName(), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality()));
													}
													else {												
														seg.splits.add(new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, Nlen, insertsize, distance, readString.substring(0, Nlen), sample.samRecord.getMateReferenceName().substring(3), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality()));
													}
											//		draw.data.add(seg);
												}		
											}
											else {
												if(sample.samRecord.getReadPairedFlag() && sample.samRecord.getMateUnmappedFlag()) {
													
													seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, Nlen, 0, distance, readString.substring(0, Nlen),"", -1, false, sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
																									
													data.add(seg);
												}
												else if (sample.samRecord.getReadPairedFlag() && !sample.samRecord.getMateUnmappedFlag()){
													try {
														if(!chr) {
															seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, Nlen, insertsize, distance, readString.substring(0, Nlen), sample.samRecord.getMateReferenceName(), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
														}
														else {												
															seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, Nlen, insertsize, distance, readString.substring(0, Nlen), sample.samRecord.getMateReferenceName().substring(3), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
														}
														data.add(seg);
													}
													catch(Exception e) {
														
														e.printStackTrace();
													}
													
												}		
											}
											
											count = 0;
											start = start+Nlen+cigar.getLength();
											readString = readString.substring(Nlen);
											length -= Nlen;
											split = true;
										}	
									else if(cigar.getOperator().compareTo(CigarOperator.SOFT_CLIP) == 0) {	
									//	count+=cigar.getLength();
									}
								}
								
								if(!split) {
									
									if(sample.samRecord.getReadPairedFlag() && sample.samRecord.getMateUnmappedFlag()) {
											
											seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, 0, distance, readString,"", -1, false, sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
											
											data.add(seg);
									}
									else if (sample.samRecord.getReadPairedFlag() && !sample.samRecord.getMateUnmappedFlag()){
										try {
											if(sample.samRecord.getMateReferenceName().length() < 4) {
												seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, insertsize, distance, readString, sample.samRecord.getMateReferenceName(), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
											}
											else {
												seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, insertsize, distance, readString, sample.samRecord.getMateReferenceName().substring(3), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
											}
											data.add(seg);
										}
										catch( Exception e) {
											
											e.printStackTrace();
										//	break;
										}										
									}
									else {
										seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, 0, distance, readString,"", -1, false, sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
										
										data.add(seg);
									}
								}
								else {
									
									if(sample.samRecord.getReadPairedFlag() && sample.samRecord.getMateUnmappedFlag()) {
										
										seg.splits.add(new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, 0, distance, readString,"", -1, false, sample.samRecord.getCigar(), sample.samRecord.getMappingQuality()));
								//		sample.data.add(seg);
								}
								else if (sample.samRecord.getReadPairedFlag() && !sample.samRecord.getMateUnmappedFlag()){
									
										if(!chr) {
											seg.splits.add(new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, insertsize, distance, readString, sample.samRecord.getMateReferenceName(), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality()));
										}
										else {
											seg.splits.add(new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, insertsize, distance, readString, sample.samRecord.getMateReferenceName().substring(3), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality()));
										}
								//		sample.data.add(seg);
									}	
								}
								}
								else {
									
									if(sample.samRecord.getMateReferenceIndex() == -1) {								
									
										seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, 0, distance, sample.samRecord.getReadString(),"", -1, false, sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
										data.add(seg);
									}
									else if (sample.samRecord.getMateReferenceIndex() > -1 ){
										
											if(!chr) {
												seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, insertsize, distance, sample.samRecord.getReadString(), sample.samRecord.getMateReferenceName(), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
											}
											else {
												try {
													if(sample.samRecord.getMateReferenceName().length() < 4) {
														seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, insertsize, distance, sample.samRecord.getReadString(), sample.samRecord.getMateReferenceName(), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
														
													}
													else {
														seg = new Segment(Main.selectedChrom+1, start, sample.samRecord.getReadName(), dir, true, length, insertsize, distance, sample.samRecord.getReadString(), sample.samRecord.getMateReferenceName().substring(3), sample.samRecord.getMateAlignmentStart(), sample.samRecord.getMateNegativeStrandFlag(), sample.samRecord.getCigar(), sample.samRecord.getMappingQuality());
													}
												}
												catch(Exception e) {
													
													e.printStackTrace();
												}
												}
											
											data.add(seg);
											
										}
									}
							
							
			//			}
				if(seg != null) {
					templevel = null;
					if(seg.level != 0) {
						continue;
					}
					if(levels.size() != 0 && levels.last().pos < seg.pos) {
						levels.clear();
						level = 1;
					}
					
					if(levels.size() != 0 && levels.first().pos < seg.pos) {				
						levelIterator = levels.iterator();
						templevel = levels.first();
						level = templevel.level;
						
						while(levelIterator.hasNext()) {
							try {
							lev = levelIterator.next();
						
							if(lev.pos < seg.pos) {
								if(lev.level < level) {
									level = lev.level;
									templevel = lev;
								}
							}
							else {
								break;
							}
							}
							catch(Exception e) {
								e.printStackTrace();
							}
						}
						
						
						if(levels.size() > 0) {
							levels.remove(templevel);
							leveladd = new Level(level, seg.end+2);
							seg.level = level;	
							levels.add(leveladd);
							
						}
						else {
							
							leveladd = new Level(1, seg.end+2);
							seg.level = 1;
							levels.add(leveladd);
						}
						
						if(seg.splits.size() > 0) {
							for(s = 0; s<seg.splits.size(); s++) {
								seg.splits.get(s).level = seg.level;														
							}
						}
					//	found = true;	
						startY = (int)((sample.getIndex()+1)*Main.drawCanvas.sampleHeight-seg.level*Main.drawCanvas.readHeight);
						startX = (seg.pos - start)*pixel;			
						seg.x = (int)(startX);
						seg.y = startY;	
					
						continue;
					}				
		//		}
				
				if(levels.size() < 5000) {
					
					levels.add(new Level(levels.size()+1, seg.end +2));
					level = levels.size();
					seg.level = level;				
					
					if(seg.splits.size() > 0) {
						for(s = 0; s<seg.splits.size(); s++) {
							seg.splits.get(s).level = seg.level;														
						}
					}	
				}
				
					startY = (int)((sample.getIndex()+1)*Main.drawCanvas.sampleHeight-seg.level*Main.drawCanvas.readHeight);
					startX = (seg.pos - start)*pixel;			
					seg.x = (int)(startX);
					seg.y = startY;	
			
				}	
				}
			}
			
		//	index = draw.covdata.length-1;	
			
			}
			
			sample.samFileReader.close();	
		
			
		}
		
		catch(Exception e) {

			e.printStackTrace();
		
		}	
			return data;
		}
		return null;
	}
	
	void searchInsSites() {
		Main.drawCanvas.insertList.clear();
		HashMap<Integer, Character> bases = new HashMap<Integer,Character>();
		bases.put(0, 'A');
		bases.put(1, 'C');
		bases.put(2, 'G');
		bases.put(3, 'T');
		bases.put(4, 'N');
		
		HashMap<Character, Integer> bases2 = new HashMap<Character, Integer>();
		bases2.put('A', 0);
		bases2.put('C', 1);
		bases2.put('G', 2);
		bases2.put('T', 3);
		
		for(int i =Main.drawCanvas.visiblestart; i<=Main.drawCanvas.drawVariables.visibleend;i++ ) {
			if(i >= Main.samples) {
				break;
			}
			ArrayList<Integer> refReadQualities = new ArrayList<Integer>();
			Iterator<Segment> it = null; //Main.drawCanvas.sampleList.get(i).readData.iterator();					
			Segment cur;
			int readLengths = 0, clipLengths = 0, consensusStart = 0;
			String sequence = "", repeat = "";
			int maxbase = 0, allbases=0, disbases=0, disbasetemp = 0;
			List<Integer> maxindex = new ArrayList<Integer>();
			int del = 0;
			Double startX;
			String readString = "";
			HashMap<Integer, Integer> startVote = new HashMap<Integer, Integer>(), endVote = new HashMap<Integer, Integer>();
			int readPos =0,startY, miscounter=0, clusterStart = 0, startReadMisSum =0, endReadMisSum = 0, clusterEnd = 0, max = 0, sum = 0, refReads = 0, altReads5 = 0, altReads3 = 0, breakPointStart = 0, breakPointEnd =0, support5=0, support3=0;
			java.util.List<Integer> positionList = Collections.synchronizedList(new ArrayList<Integer>());
			java.util.List<Object[]> startReadList = Collections.synchronizedList(new ArrayList<Object[]>());
			java.util.List<Object[]> endReadList = Collections.synchronizedList(new ArrayList<Object[]>());
			String rate3 = "", rate5 = "";
			
			while (it.hasNext()) {
				
				try {
					cur = (Segment)it.next();
					
					if(cur.pos+cur.length < start ) {
						continue;
					}
					if(cur.pos > Main.drawCanvas.end ) {
						break;
					}
				
					miscounter=0;
					clusterStart = 0;
					clusterEnd = 0;
					startVote.clear();
					endVote.clear();
					refReads = 0;					
					refReadQualities.clear();
					startReadList.clear();
					endReadList.clear();
				
					if((cur.cigar.toString().contains("S"))) {
							
							breakPointStart = breakPointEnd = 0;
							clusterStart = cur.pos;
							clusterEnd = cur.pos + 150;
							
						if(cur.cigar.toString().endsWith("S")) {
							clipLengths = cur.cigar.getCigarElement(cur.cigar.getCigarElements().size()-1).getLength();
							readLengths = cur.length;
							Object[] read = {(cur.pos+cur.length)-cur.cigar.getCigarElement(cur.cigar.getCigarElements().size()-1).getLength()-1, cur.readString.substring(cur.readString.length()-cur.cigar.getCigarElement(cur.cigar.getCigarElements().size()-1).getLength(), cur.readString.length()) };
							endReadList.add(read);
						//	endVote.put((cur.pos+cur.length)-cur.cigar.getCigarElement(cur.cigar.getCigarElements().size()-1).getLength(), 1);
						}
						else {
							clipLengths = cur.cigar.getCigarElement(0).getLength();
							readLengths = cur.length;
							
							Object[] read = {cur.pos, cur.readString.substring(0, cur.cigar.getCigarElement(0).getLength()) };
							startReadList.add(read);
					//		startVote.put(cur.pos+cur.cigar.getCigarElement(0).getLength(), 1);
						}
						if(it.hasNext()) {
							cur = (Segment)it.next();
						}
						else {
							break;
						}
						
				//	while(cur.pos < clusterEnd) {
						
					while(clipLengths/(double)readLengths > 0.10 && cur.pos < clusterEnd) {	
						
						if(cur.direction && cur.cigar.getCigarElement(0).getOperator().compareTo(CigarOperator.SOFT_CLIP)== 0) {
						
							if(cur.quality >= 10) {
								
								clipLengths += cur.cigar.getCigarElement(0).getLength();
								readLengths += cur.length;
								Object[] read = {cur.pos, cur.readString.substring(0, cur.cigar.getCigarElement(0).getLength()) };
								startReadList.add(read);		
								
							}
						//	System.out.println("Start: " +(MethodLibrary.formatNumber(cur.pos+cur.cigar.getCigarElement(0).getLength())));
						}
						else if(!cur.direction && cur.cigar.toString().endsWith("S")) {
							
							if(cur.quality >= 10) {
								clipLengths += cur.cigar.getCigarElement(cur.cigar.getCigarElements().size()-1).getLength();
								readLengths += cur.length;
								Object[] read = {(cur.pos+cur.length)-cur.cigar.getCigarElement(cur.cigar.getCigarElements().size()-1).getLength()-1, cur.readString.substring(cur.readString.length()-cur.cigar.getCigarElement(cur.cigar.getCigarElements().size()-1).getLength(), cur.readString.length()) };
								endReadList.add(read);
								
							}
						//	System.out.println("End: " +(MethodLibrary.formatNumber((cur.pos+cur.length)-cur.cigar.getCigarElement(0).getLength())));
						}
						else {
							try {
								readLengths += cur.length;
								refReadQualities.add(cur.quality);
								refReads++;
							}
							catch(Exception e) {
								e.printStackTrace();
								System.out.println(cur.cigarString);
							}
						}
						if(it.hasNext()) {
							cur = (Segment)it.next();
						}
						else {
							break;
						}
				//		System.out.println(clipLengths/(double)readLengths);
						
					}
					
					if(startReadList.size() < 3 && endReadList.size() < 3) {
						
						continue;
					}
					
					clusterEnd = cur.pos + cur.length;
									
					Collections.sort(refReadQualities);
					if(refReadQualities.size() > 0 && refReadQualities.get(refReadQualities.size()/2) < 10) {
						
						continue;
					}
					startReadMisSum = 0;
					
					for(int r = 0; r< startReadList.size(); r++) {
						miscounter = 0;
						readPos = (Integer)startReadList.get(r)[0];
						readString = (String)startReadList.get(r)[1];
						
						for(int c = 0; c< readString.length(); c++) {												
							
							if(readPos+c > Main.chromDraw.getSeqStart() && (readPos+c)-Main.chromDraw.getSeqStart() <Main.chromDraw.getSequence().length()) {
								if(Main.chromDraw.getSequence().charAt((readPos+c)-Main.chromDraw.getSeqStart()) != readString.charAt(c)) {
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
							
							if(readPos+c > Main.chromDraw.getSeqStart() && (readPos+c)-Main.chromDraw.getSeqStart() <Main.chromDraw.getSequence().length()) {
								if(Main.chromDraw.getSequence().charAt((readPos+c)-Main.chromDraw.getSeqStart()+1) != readString.charAt(c)) {
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
					
					String seq;
					
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
					
					int[] adder = {breakPointStart, breakPointEnd};
					
					Main.drawCanvas.insertList.add(adder);
					Main.drawCanvas.repaint();
				
		//			sequence = Main.chromDraw.getSeq(clusterStart, clusterEnd, Main.referenceFile);
					
					
			
					
					int[][] consensus5;
					char[] result5 = {};
					int[][] consensus3;
					char[] result3 = {};
					
				//	System.out.println(breakPointStart +" " +breakPointEnd);
					
					if(breakPointStart > -1) {
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
									if(readPos+s-del-consensusStart < consensus5[0].length) {									
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
								}
							
						//	System.out.println();
						}
				//		System.out.print(new String(result5) +"\t" +disbases/(double)allbases);
				//		System.out.println();
						rate5 = ""+MethodLibrary.round(disbases/(double)allbases, 2);
					}
					
					if(breakPointEnd > -1) {
						Collections.sort(endReadList, new MethodLibrary.ReadListSorter());
						String readlen = (String)endReadList.get(endReadList.size()-1)[1];
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
						
						rate3 = ""+MethodLibrary.round(disbases/(double)allbases, 2);
					}
					
					if(breakPointStart > -1 && breakPointEnd > -1) {
						if(breakPointStart <= breakPointEnd) {
							repeat = sequence.substring(breakPointStart-clusterStart, breakPointEnd-clusterStart);
							System.out.println(cur.chrom +"\t" +breakPointStart +"\t" +breakPointEnd +"\t" +(breakPointEnd-breakPointStart) +"\t-\t+\t" +repeat +"\t" +(new String(result5)) +"\t" +rate5 +"\t" +(new String(result3)) +"\t" +rate3 +"\t" +(altReads5+altReads3) +"\t" +refReads +"\t" +support5 +"\t" +support3 );
						}
						else {
							repeat = sequence.substring(breakPointEnd-clusterStart, breakPointStart-clusterStart);
							System.out.println(cur.chrom +"\t" +breakPointEnd +"\t" +breakPointStart +"\t" +(breakPointStart-breakPointEnd) +"\t+\t-\t" +repeat +"\t" +(new String(result5)) +"\t" +rate5 +"\t" +(new String(result3)) +"\t" +rate3 +"\t" +(altReads5+altReads3) +"\t" +refReads +"\t" +support3 +"\t" +support5 );
							
						}
					}
					else {
						if(breakPointEnd == -1) {
							repeat = sequence.substring(breakPointStart-clusterStart, breakPointStart-clusterStart + 20);
							System.out.println(cur.chrom +"\t" +breakPointStart +"\t" +(breakPointStart+1) +"\t0\t-\t-\t" +repeat +"\t" +(new String(result5)) +"\t" +rate5 +"\t-\t-\t" +(altReads5+altReads3) +"\t" +refReads +"\t" +support5 +"\t" +support3 );
						}
						else {
							repeat = sequence.substring(breakPointEnd-clusterStart -20, breakPointEnd-clusterStart);
							System.out.println(cur.chrom +"\t" +breakPointEnd   +"\t" +(breakPointEnd+1)   +"\t0\t+\t+\t" +repeat +"\t-\t-\t" +(new String(result3)) +"\t" +rate3 +"\t" +(altReads5+altReads3) +"\t" +refReads +"\t" +support5 +"\t" +support3 );
							
						}
					}
					
				}
				}
				catch(Exception e) {
					e.printStackTrace();
				}
			}
			
		}
		System.out.println("--------");
		
	}
	*/
}
