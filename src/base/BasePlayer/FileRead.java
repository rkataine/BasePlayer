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
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.index.Block;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.samtools.cram.ref.ReferenceSource;
import java.awt.Dimension;
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
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.SwingWorker;

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
//	TabixReader tabixreader;
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
	private Short refcalls;
	private short altallele;
	private Short altcalls;
	private String altbase2;
	private String[] headersplit;
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
	private int seqstart;
	private int seqend;
	public boolean firstSample = false;
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
String setVCFFileStart(String chrom, int start, int end, Sample sample) {
	try {	Index index = null;	
			if(sample.getVCFReader() != null) {
				index = IndexFactory.loadIndex(sample.getTabixFile()+".idx");				
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
					
						sample.getVCFInput().seek(0);
						sample.getVCFInput().seek(blocks.get(0).getStartPosition());
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
	Main.opensamples.setText("Open samples");
	head.putNext(null);
	current = null;
	Main.drawCanvas.current = FileRead.head;
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
					if(vcf) {
						line = sample.getVCFReader().readLine();
						if(line != null && line.startsWith("#")) {
							continue;
						}
					}
					else {
						
						line = sample.getVCFInput().readLine();
					}
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
				
				if(sample.getVCFInput() != null) {
						if(sample.getVCFInput().getFilePointer() > sample.vcfEndPos ) {			
							
							break;
						}
				}
				else {
					
					if(Integer.parseInt(split[1]) > end || !split[0].equals(searchChrom)) {
						break;
					}
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
		if(sample.getVCFReader() != null) {
			sample.getVCFReader().close();
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
		Main.drawCanvas.loadbarAll = 0;
		Main.drawCanvas.loadBarSample= 0;
	}
	
	public ArrayList<Gene> getExons(String chrom) {		
		
		ArrayList<Gene> transcriptsTemp = new ArrayList<Gene>();
		
		try {			
			if(Main.genomehash.size() == 0 || Main.genomehash.get(Main.defaultGenome).size() == 0) {
				return new ArrayList<Gene>();
			}
			
			ChromDraw.exonReader = new TabixReader(Main.genomehash.get(Main.defaultGenome).get(Main.annotation).getCanonicalPath());
			if(chrom == null) {
				return null;
			}
			TabixReader.Iterator exonIterator = null;
				try {	
					if(!ChromDraw.exonReader.getChromosomes().contains(chrom)) {
						
						return new ArrayList<Gene>();
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
		      	 
		      
		      try {
		    	  if(currentSample.samFile.getName().endsWith(".cram")) {
		    		
		    		 CRAMReader = new CRAMFileReader(currentSample.samFile, new File(currentSample.samFile.getCanonicalFile() +".crai"),new ReferenceSource(Main.ref), ValidationStringency.LENIENT);
		    		 
		    		 if(CRAMReader.getFileHeader().getTextHeader().contains("SN:chr")) {
		    			 currentSample.chr = "chr";
		    		 }
		    		
		    			CRAMReader.close();
		    	  }
		    	  else {
		    		 
		    	//	  samFileReader = SamReaderFactory.makeDefault().open(currentSample.samFile);
		    		  samFileReader = new SAMFileReader(currentSample.samFile);
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
		      	added = true;
				  Main.samples++;
		         }   
	  		 
	  	 }
	  	if(!added) {
	  		return;
	  	}
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
		
		// Draw.updateReads = true;
		 //Main.drawCanvas.splits.get(0).updateReads = true;
		 Draw.updatevars = true;
		 Main.drawCanvas.repaint();
}	  	 
	
private void readVCF(File[] files) {
		try {	 
	  boolean cram = false;
  	  String line;
  	  Main.drawCanvas.loading("Loading variants...");
  	  File[] addDir;
  	  int sampletemp = Main.samples;
  	  Boolean added = false;
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
  		if(!files[fi].exists()) {
  			continue;
  		}
  		if(!files[fi].getName().endsWith(".vcf") && !files[fi].getName().endsWith(".vcf.gz")) {
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
  	    		if(cancelfileread) {
  	    			current = null;
  	    			FileRead.head.putNext(null);
  	    			return;
  	    		}
  	    		 addSample = new Sample(addDir[f].getName(), (Short)Main.samples, addDir[f].getCanonicalPath());      
  	    		 Main.drawCanvas.sampleList.add(addSample);
  	        	 Main.varsamples++;
  	        	 Main.samples++;
  	        	 added = true;
  	        	 VariantHandler.commonSlider.setMaximum(Main.varsamples);
  	        	 VariantHandler.commonSlider.setUpperValue(Main.varsamples);
  	        	 VariantHandler.geneSlider.setMaximum(Main.varsamples);
  	        	 sampleString.append(addSample.getName() +";");
  	    		
  	    		}
  	    		File[] bamfiles = files[fi].listFiles(new FilenameFilter() {
  	    		 
  	    	     public boolean accept(File dir, String name) {
  	    	        return name.toLowerCase().endsWith(".bam") || name.toLowerCase().endsWith(".cram") || name.toLowerCase().endsWith(".list");
  	    	     }
  	    	   	 });
  	    		/*File[] cramfiles = files[fi].listFiles(new FilenameFilter() {
  	  	    		 
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
  	    		  		// if(bamfiles[i].getName().indexOf(".bam") > 0) {
  	    		  			index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".")));
  	    		  		if(bamfiles[i].getName().endsWith(".cram") ) {
  	    		  			cram = true;
  	    		  		}
  	    		  	/*	 }
  	    		  		 else  if(bamfiles[i].getName().indexOf(".cram") > 0) {
  	    		  			 index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".cram")));
  	    		  			cram = true;
  	    		  		 }
  	    		  		 else {
  	    		  			index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".list")));
  	    		  		 }*/
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
  	    	  	 }
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
  	    	fileindex =fi;
  	    	
        	addSample = new Sample(files[fi].getName(), (Short)Main.samples, files[fi].getCanonicalPath());        	
        	
        	 Main.drawCanvas.sampleList.add(addSample);
          	 Main.varsamples++;
          	 Main.samples++;
          	 added = true;
          	 VariantHandler.commonSlider.setMaximum(Main.varsamples);
          	 VariantHandler.commonSlider.setUpperValue(Main.varsamples);
          	 VariantHandler.geneSlider.setMaximum(Main.varsamples);
          	 sampleString.append(addSample.getName() +";");
       	  }     	
  	 }
  	 if(!added) {
  		 return;
  	 }
  	 if(fileindex > -1) {
     File[] bamfiles = files[fileindex].getParentFile().listFiles(new FilenameFilter() {
     public boolean accept(File dir, String name) {
        return name.toLowerCase().endsWith(".bam") || name.toLowerCase().endsWith(".cram") || name.toLowerCase().endsWith(".list");
     }
   	 });
    /* File[] cramfiles = files[fileindex].getParentFile().listFiles(new FilenameFilter() {
         public boolean accept(File dir, String name) {
            return name.toLowerCase().endsWith(".cram");
         }
       	 });*/
  	 if(bamfiles.length > 0) {
  		
  		int index = -1, sampleindex;
  		 for(int i = 0; i<bamfiles.length; i++) {
  			 cram = false;
	  		 sampleindex = 0;
	  		 index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".")));
	  		if(bamfiles[i].getName().endsWith(".cram")){
	  			cram = true;
	  		}
	  		 /*if(bamfiles[i].getName().indexOf(".bam") > 0) {
	  			index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".bam")));
	  				  			
	  		 }
	  		 else if(bamfiles[i].getName().indexOf(".cram") > 0){
	  			 index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".cram")));
	  			cram = true;
	  		 }
	  		 else {
	  			index = sampleString.indexOf(bamfiles[i].getName().substring(0,bamfiles[i].getName().indexOf(".list")));
	  		 }*/
	  	 	 if (index < 0) continue; 
	  	 	 
	  	 	 for(char letter : sampleString.substring(0, index).toString().toCharArray()) {
	  	 		 if(letter == ';') sampleindex++;
	  	 	 }
	  	 	Main.drawCanvas.bam = true;
	  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).CRAM = cram;
	  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile = new File(bamfiles[i].getAbsolutePath());	 
	  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).resetreadHash();
	  	 	//samFileReader = SamReaderFactory.makeDefault().open(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
	  	 	samFileReader = new SAMFileReader(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
	  	    if(samFileReader.getFileHeader().getTextHeader().contains("SN:chr")) {
      	    	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).chr = "chr";
      	    }
	      	samFileReader.close();
	  	 }  	
  	 }
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
	  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile = new File(bamfiles[i].getAbsolutePath());	 
	  	 	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).resetreadHash();
	  	 	//samFileReader = SamReaderFactory.makeDefault().open(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
	  	 	samFileReader = new SAMFileReader(Main.drawCanvas.sampleList.get(sampleindex+sampletemp).samFile);
	  	    if(samFileReader.getFileHeader().getTextHeader().contains("SN:chr")) {
     	    	Main.drawCanvas.sampleList.get(sampleindex+sampletemp).chr = "chr";
     	    }
	      	samFileReader.close();
	  	 }  	
  	 }*/
  	/*	 else {
  	
  		 
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
  	 }*/
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
	    	/*  try {
	    		  tabixreader = new TabixReader(Main.drawCanvas.sampleList.get(i).getTabixFile());
	    	  }
	    	  catch(Exception ex) {
	    		  JOptionPane.showMessageDialog(Main.chromDraw, "Index file (tbi) not found for " +Main.drawCanvas.sampleList.get(i).getName(), "Error", JOptionPane.ERROR_MESSAGE);
	    		  ErrorLog.addError(ex.getStackTrace());
	    		  ex.printStackTrace();
	    	  }*/
		//      iterator=null;   
	    	  BufferedReader reader;
	    	  GZIPInputStream gzip = null;
	    	  if(Main.drawCanvas.sampleList.get(i).getTabixFile().endsWith(".gz")) {
			     gzip = new GZIPInputStream(new FileInputStream(Main.drawCanvas.sampleList.get(i).getTabixFile()));
			      reader = new BufferedReader(new InputStreamReader(gzip));		
	    	  }
	    	  else {
	    		  reader = new BufferedReader(new FileReader(Main.drawCanvas.sampleList.get(i).getTabixFile()));		
	    	  }
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
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
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
							}					
							current.getTranscripts().add(prevGene.getTranscripts().get(0));
							current.getTranscripts().add(gene.getTranscripts().get(0));
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
		if(linecounter  == 0 || pos -linecounter > (Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart)/100) {
		
				 Main.drawCanvas.loadBarSample = (int)(((pos-Main.drawCanvas.variantsStart)/(double)(Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart))*100);     
				 Draw.updatevars = true;
				 linecounter = pos; 
	
				
		}
		
		
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
				if(altbase.contains("*") || (altbase2!=null && altbase2.contains("*"))) {
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
								current.putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, split[2], currentSample, current, null));		
							}
							else {
								current.putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, null, currentSample, current, null));		
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
								current.getPrev().putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, split[2], currentSample, current.getPrev(), current));
							}
							else {
								current.getPrev().putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, null, currentSample, current.getPrev(), current));
								
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
		if(linecounter  == 0 || pos -linecounter > (Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart)/100) {
			
			 Main.drawCanvas.loadBarSample = (int)(((pos-Main.drawCanvas.variantsStart)/(double)(Main.drawCanvas.variantsEnd-Main.drawCanvas.variantsStart))*100);     
			 Draw.updatevars = true;
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
					//TODO jaettuna kolmella ei v�ltt�m�tt� pid� paikkaansa, jos format kent�ss� on esim: GT:GQ:FREQ:AD:RD....
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
		//		dpindex = split[7].indexOf("DP4");
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
		
			
		
		
		
		if(sample.maxCoverage < refcalls+altcalls) {
			sample.maxCoverage = (short)(refcalls+altcalls);
			
		}		
		
		while(current.getNext() != null && current != null && current.getPosition() < pos ){	
		
			current = current.getNext();		      				
		}		
		
		if(current.getNext() == null && current.getPosition() < pos) {	
			
			try {				
				if(!split[2].equals(".")) {
					current.putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, split[2], sample, current, null));		
				}
				else {
					current.putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, null, sample, current, null));		
					
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
					current.getPrev().putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, split[2], sample, current.getPrev(), current));
				}
				else {
					current.getPrev().putNext(new VarNode(pos, (byte)split[3].charAt(0), altbase, (short)(refcalls+altcalls), altcalls, genotype, quality, null, sample, current.getPrev(), current));
					
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
					//samFileReader = SamReaderFactory.make().open(file);
					samFileReader = new SAMFileReader(file);
				}
			}
			else {
				if(readClass.sample.CRAM) {					
					
					 CRAMReader = new CRAMFileReader(readClass.sample.samFile, new File(readClass.sample.samFile.getCanonicalFile() +".crai"),new ReferenceSource(Main.ref),Main.referenceFile, ValidationStringency.SILENT);
					 if(CRAMReader != null && !CRAMReader.hasIndex()) {			
							JOptionPane.showMessageDialog(null, "Index file is missing (.crai) for " +readClass.sample.samFile.getName(), "Error", JOptionPane.ERROR_MESSAGE);		    		
							return null;
						}
				}
				else {
					//samFileReader = SamReaderFactory.make().open(readClass.sample.samFile);
					samFileReader = new SAMFileReader(readClass.sample.samFile);
					if(samFileReader != null && !samFileReader.hasIndex()) {			
						JOptionPane.showMessageDialog(null, "Index file is missing (.bai) for " +readClass.sample.samFile.getName(), "Error", JOptionPane.ERROR_MESSAGE);		    		
						return null;
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
				//System.out.println(value.hasNext());
				return value;
			}
			else {
				return samFileReader.queryOverlapping(chrom, startpos, endpos);	
			}
		}
		catch(Exception e) {			
			try {
				if(readClass.sample.CRAM) {
					QueryInterval[] interval = { new QueryInterval(CRAMReader.getFileHeader().getSequence(readClass.sample.chr + "M").getSequenceIndex(), startpos, endpos) };
					
					Iterator<SAMRecord> value = CRAMReader.query(interval, false);
					//System.out.println(value.hasNext());
					return value;
				}
				else {
					return samFileReader.queryOverlapping(readClass.sample.chr + "M", startpos, endpos);	
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
						
						splitIndex.drawReference = new ReferenceSeq(splitIndex.chrom,startpos-searchwindow-1,  endpos+searchwindow +200, Main.referenceFile);
						
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
					if(splitIndex.drawReference == null) {
						splitIndex.drawReference = new ReferenceSeq(splitIndex.chrom,startpos-searchwindow-1,  endpos+searchwindow +200, Main.referenceFile);
						
					}
					else {
						splitIndex.drawReference.append( endpos+200);
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
					splitIndex.drawReference.appendToStart( startpos-1);
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
			
			while(bamIterator != null && bamIterator.hasNext()) {	
			
				try {	
					
					if(cancelreadread) {
			  			return null;
			  		}
					try {
						samRecord = bamIterator.next(); 	
						
					}
					catch(htsjdk.samtools.SAMFormatException ex) {
				//		ex.printStackTrace();
						
					}
					
				if(samRecord.getReadUnmappedFlag()) {					
					continue;
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
								
								ReadNode addNode = new ReadNode(samRecord, readClass.complete, chrom, readClass.sample, splitIndex,readClass, mismatches);
								readClass.setFirstRead(addNode);
								startY = (int)(((readClass.sample.getIndex()+1)*Main.drawCanvas.drawVariables.sampleHeight-Main.drawCanvas.bottom-(readClass.readHeight+2)));					
						//		addNode.setRectBounds((int)((addNode.getPosition()-start)*pixel), startY, (int)(samRecord.getReadLength()*pixel), Main.drawCanvas.readHeight);
														
								addNode.setPrev(null);
								addNode.setNext(readClass.getHeadAndTail().get(0)[headnode]);
								readClass.getHeadAndTail().get(0)[headnode].setPrev(addNode);
								lastAdded = addNode;
						
							}
							else {
								
								ReadNode addNode = new ReadNode(samRecord, readClass.complete, chrom, readClass.sample, splitIndex,readClass, mismatches);
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
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
	readClass.setCoverageEnd(readClass.getCoverageStart()+coverages.length);
	readClass.setCoverages(coverages);
	return coverages;
	
}
	
	
java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> getMismatches(SAMRecord samRecord, Reads readClass, double[][] coverageArray) {
	
	
	java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> mismatches = null;
	String MD = samRecord.getStringAttribute("MD");	
	String SA = samRecord.getStringAttribute("SA");
	
	try {
	
	if(MD == null) {
		
		if(readClass.sample.MD) {
			readClass.sample.MD = false;
		}
	}
	if(readClass.sample.MD || MD!=null && (samRecord.getCigarString().contains("S") && !Settings.softclips.isSelected())/* && ((!samRecord.getCigarString().contains("S") && SA == null) || SA !=null)*/) {		
		
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
					
						if(Settings.softclips.isSelected()) {
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
						
						readClass.getCoverages()[(readstart+readpos)- readClass.getCoverageStart()][Main.baseMap.get(""+bases.charAt(mispos))]++;
						
							if(samRecord.getBaseQualityString().length() != 1 && (int)qualities.charAt(mispos)-readClass.sample.phred < Settings.getBaseQuality() ) {
								
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
		while(splitIndex.drawReference == null) {
			Thread.sleep(100);		
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
						if(((readstart+r)-splitIndex.drawReference.getStartPos()-1) > splitIndex.drawReference.getSeq().length-1) {
							splitIndex.drawReference.append(samRecord.getUnclippedEnd()+1000);
						}
						if(((readstart+r)-splitIndex.drawReference.getStartPos()-1) < 0) {
							splitIndex.drawReference.appendToStart(samRecord.getUnclippedStart()-1000);
						}
						
						if(samRecord.getReadBases()[mispos] !=  splitIndex.drawReference.getSeq()[((readstart+r)-splitIndex.drawReference.getStartPos()-1)]) {	
							readClass.getCoverages()[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get(Main.getBase.get(samRecord.getReadBases()[mispos]))]++;
							
							if(samRecord.getBaseQualityString().length() != 1 && (int)samRecord.getBaseQualityString().charAt(mispos)-readClass.sample.phred < Settings.getBaseQuality() ) {
								
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
						readClass.getCoverages()[(readstart+readpos+d)- readClass.getCoverageStart()][Main.baseMap.get("D")]++;
						
					}
					readpos+=element.getLength();
				}
				else if(element.getOperator().compareTo(CigarOperator.INSERTION)== 0) {										
					readClass.getCoverages()[(readstart+readpos)- readClass.getCoverageStart()][Main.baseMap.get("I")]++;
					mispos+=element.getLength();					
				
				}
				else if(element.getOperator().compareTo(CigarOperator.SOFT_CLIP)== 0 ) {
						if(Settings.softclips.isSelected()) {
							for(int r = readpos; r<readpos+element.getLength(); r++) {
								if(((readstart+r)-splitIndex.drawReference.getStartPos()-1) > splitIndex.drawReference.getSeq().length-1) {
									splitIndex.drawReference.append(samRecord.getUnclippedEnd()+1000);
								}
								if(((readstart+r)-splitIndex.drawReference.getStartPos()-1) < 0) {
									splitIndex.drawReference.appendToStart(samRecord.getUnclippedStart()-1000);
								}
								/*	if(((readstart+r)-readClass.reference.getStartPos()-1) < 0) {
									continue;
								}*/
								if(samRecord.getReadBases()[mispos] != splitIndex.drawReference.getSeq()[((readstart+r)-splitIndex.drawReference.getStartPos()-1)]) {													
									if(SA == null) {
									
										readClass.getCoverages()[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get(Main.getBase.get(samRecord.getReadBases()[mispos]))]++;
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
					if(((readstart+r)-splitIndex.drawReference.getStartPos()-1) > splitIndex.drawReference.getSeq().length-1) {
						splitIndex.drawReference.append(samRecord.getUnclippedEnd()+1000);
					}
					if(((readstart+r)-splitIndex.drawReference.getStartPos()-1) < 0) {
						splitIndex.drawReference.appendToStart(samRecord.getUnclippedStart()-1000);
					}
				/*	if(((readstart+r)-readClass.reference.getStartPos()-1) > readClass.reference.getSeq().length-1) {
						readClass.reference.append(samRecord.getUnclippedEnd()+1000);
					}
					if(((readstart+r)-readClass.reference.getStartPos()-1) < 0) {
						continue;
					}*/
					if(samRecord.getReadBases()[r] != splitIndex.drawReference.getSeq()[((readstart+r)-splitIndex.drawReference.getStartPos()-1)]) {									
						
						if((readstart+r)-readClass.getCoverageStart() < readClass.getCoverages().length-1) {
							readClass.getCoverages()[(readstart+r)- readClass.getCoverageStart()][Main.baseMap.get(Main.getBase.get(samRecord.getReadBases()[r]))]++;
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
		  add.setName(track.file.getName().substring(0, 10) +"...");
	  }
	  else {
		  VariantHandler.tabs.add(track.file.getName(), addpane);
		  add.setName(track.file.getName());
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

void readSam(String chrom, Reads readClass, SAMRecord samRecordBuffer, java.util.ArrayList<java.util.Map.Entry<Integer,Byte>> mismatches) {
	try {
		
		if(readClass.getReads().isEmpty()) {
			
			addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, readClass.sample, splitIndex,readClass, mismatches);			
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
											addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, readClass.sample, splitIndex,readClass, mismatches);	
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
								addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, readClass.sample, splitIndex,readClass, mismatches);	
								
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
								addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, readClass.sample, splitIndex,readClass, mismatches);	
								
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
						addNode = new ReadNode(samRecordBuffer, readClass.complete, chrom, readClass.sample, splitIndex,readClass, mismatches);	
						
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
	VariantHandler.table.genearray.clear();
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
		split.getGenes().get(g).mutations = 0;
		split.getGenes().get(g).missense = 0;
		split.getGenes().get(g).nonsense = 0;
		split.getGenes().get(g).synonymous = 0;
		split.getGenes().get(g).intronic = 0;
		split.getGenes().get(g).utr = 0;
		split.getGenes().get(g).samples.clear();
		split.getGenes().get(g).varnodes.clear();
		split.getGenes().get(g).transcriptString = new StringBuffer();
		/*for(int i = 0; i<split.getGenes().get(g).getTranscripts().size(); i++) {
			split.getGenes().get(g).getTranscripts().get(i).mutations = 0;
			split.getGenes().get(g).getTranscripts().get(i).missense = 0;
			split.getGenes().get(g).getTranscripts().get(i).nonsense = 0;
			split.getGenes().get(g).getTranscripts().get(i).synonymous = 0;
			split.getGenes().get(g).getTranscripts().get(i).intronic = 0;
			split.getGenes().get(g).getTranscripts().get(i).utr = 0;
			split.getGenes().get(g).getTranscripts().get(i).samples.clear();
			split.getGenes().get(g).getTranscripts().get(i).varnodes.clear();
		}*/
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
	Main.drawCanvas.loadBarSample= 0;
	while(true) {
		
		if(cancelvarcount) {
			cancelFileRead();
			vardraw = null;
  			break;
  		}
		
		if(vardraw == null) {
			
			if(VariantHandler.writetofile.isSelected()) {
			/*	if(VariantHandler.table.transarray.size() > 0) {
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
				*/
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
						split.getGenes().get(g).mutations = 0;
						split.getGenes().get(g).missense = 0;
						split.getGenes().get(g).nonsense = 0;
						split.getGenes().get(g).synonymous = 0;
						split.getGenes().get(g).intronic = 0;
						split.getGenes().get(g).utr = 0;
						split.getGenes().get(g).samples.clear();
						split.getGenes().get(g).varnodes.clear();
						split.getGenes().get(g).transcriptString = new StringBuffer();
					/*	for(int i = 0; i<split.getGenes().get(g).getTranscripts().size(); i++) {
							split.getGenes().get(g).getTranscripts().get(i).mutations = 0;
							split.getGenes().get(g).getTranscripts().get(i).missense = 0;
							split.getGenes().get(g).getTranscripts().get(i).nonsense = 0;
							split.getGenes().get(g).getTranscripts().get(i).synonymous = 0;
							split.getGenes().get(g).getTranscripts().get(i).intronic = 0;
							split.getGenes().get(g).getTranscripts().get(i).utr = 0;
							split.getGenes().get(g).getTranscripts().get(i).samples.clear();
							split.getGenes().get(g).getTranscripts().get(i).varnodes.clear();
						}*/
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
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
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
	VariantHandler.table.genearray.clear();
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
		split.getGenes().get(g).mutations = 0;
		split.getGenes().get(g).missense = 0;
		split.getGenes().get(g).nonsense = 0;
		split.getGenes().get(g).synonymous = 0;
		split.getGenes().get(g).intronic = 0;
		split.getGenes().get(g).utr = 0;
		split.getGenes().get(g).samples.clear();
		split.getGenes().get(g).varnodes.clear();
		split.getGenes().get(g).transcriptString = new StringBuffer();
	/*	for(int i = 0; i<split.getGenes().get(g).getTranscripts().size(); i++) {
			split.getGenes().get(g).getTranscripts().get(i).mutations = 0;
			split.getGenes().get(g).getTranscripts().get(i).missense = 0;
			split.getGenes().get(g).getTranscripts().get(i).nonsense = 0;
			split.getGenes().get(g).getTranscripts().get(i).synonymous = 0;
			split.getGenes().get(g).getTranscripts().get(i).intronic = 0;
			split.getGenes().get(g).getTranscripts().get(i).utr = 0;
			split.getGenes().get(g).getTranscripts().get(i).samples.clear();
			split.getGenes().get(g).getTranscripts().get(i).varnodes.clear();
		}*/
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
	Main.drawCanvas.loadBarSample= 0;
	while(true) {
		
		if(cancelvarcount || cancelfileread) {
			cancelFileRead();
			vardraw = null;
  			break;
  		}
		if(vardraw == null) {
			if(VariantHandler.writetofile.isSelected()) {
		/*		if(VariantHandler.table.transarray.size() > 0) {
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
				}				*/
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
						split = Main.drawCanvas.splits.get(0);
						split.getGenes().get(g).mutations = 0;
						split.getGenes().get(g).missense = 0;
						split.getGenes().get(g).nonsense = 0;
						split.getGenes().get(g).synonymous = 0;
						split.getGenes().get(g).intronic = 0;
						split.getGenes().get(g).utr = 0;
						split.getGenes().get(g).samples.clear();
						split.getGenes().get(g).varnodes.clear();
						split.getGenes().get(g).transcriptString = new StringBuffer();
					/*	for(int i = 0; i<Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().size(); i++) {
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).mutations = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).missense = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).nonsense = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).synonymous = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).intronic = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).utr = 0;
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).samples.clear();
							Main.drawCanvas.splits.get(0).getGenes().get(g).getTranscripts().get(i).varnodes.clear();
						}*/
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
	Main.drawCanvas.loadbarAll = 0;
	Main.drawCanvas.loadBarSample= 0;
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
		/*if(VariantHandler.table.transarray.size() > 0) {
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
		}	*/	
	}
	Map.Entry<String, ArrayList<SampleNode>> entry = null;
	String base = null, amino = null;
	int pretrack = -1;
	
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
		//	}
		}
	if(mutcount == 0) {
		continue;
	}
	if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {
		if(Main.drawCanvas.clusterNodes.size() == 0) {
			ClusterNode cluster = new ClusterNode();
			cluster.nodecount=mutcount;
			cluster.ID = vardraw.clusterId;
			cluster.varnodes.add(vardraw);			
			cluster.width = 1;
			Main.drawCanvas.clusterNodes.add(cluster);
		}
		else if(vardraw.clusterId-1 >= Main.drawCanvas.clusterNodes.size()) {
			Main.drawCanvas.clusterNodes.get(vardraw.clusterId-2).width = Main.drawCanvas.clusterNodes.get(vardraw.clusterId-2).varnodes.get(Main.drawCanvas.clusterNodes.get(vardraw.clusterId-2).varnodes.size()-1).getPosition() -Main.drawCanvas.clusterNodes.get(vardraw.clusterId-2).varnodes.get(0).getPosition() +1;
			
			ClusterNode cluster = new ClusterNode();
			cluster.nodecount+= mutcount;
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
	
			pretrack = -1;
			for(int i = 0; i<vardraw.getBedHits().size(); i++) {
				if(pretrack != vardraw.getBedHits().get(i).getTrack().trackIndex) {
					vardraw.getBedHits().get(i).getTrack().getTable().variants += mutcount;
					pretrack = vardraw.getBedHits().get(i).getTrack().trackIndex;
				}
				if(!vardraw.getBedHits().get(i).getTrack().getTable().bedarray.contains(vardraw.getBedHits().get(i))) {
					vardraw.getBedHits().get(i).mutations+=mutcount;
					vardraw.getBedHits().get(i).getTrack().getTable().bedarray.add(vardraw.getBedHits().get(i));
					vardraw.getBedHits().get(i).getTrack().getTable().setPreferredSize(new Dimension(vardraw.getBedHits().get(i).getTrack().getTable().tablescroll.getViewport().getWidth(), (vardraw.getBedHits().get(i).getTrack().getTable().getTableSize()+2+samplecount)*vardraw.getBedHits().get(i).getTrack().getTable().rowHeight));										
					vardraw.getBedHits().get(i).getTrack().getTable().revalidate();		
					vardraw.inVarList = true;
				}
				else {
					vardraw.getBedHits().get(i).mutations+=mutcount;
					vardraw.inVarList = true;
				}
			}
		}		
	//}
	// INTRONIC
	
	if(VariantHandler.intronic.isSelected() && vardraw.isInGene() && vardraw.getTranscripts() != null && vardraw.getExons() == null) {
		
		for(int t = 0; t<vardraw.getTranscripts().size(); t++) {				
			
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
				if(!vardraw.getTranscripts().get(t).getGene().samples.contains(entry.getValue().get(i).getSample())) {
					vardraw.getTranscripts().get(t).getGene().samples.add(entry.getValue().get(i).getSample());
				}
				
			}
			if(vardraw.getTranscripts().get(t).getGene().mutations == 0 && vardraw.getTranscripts().get(t).getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
				VariantHandler.table.addEntry(vardraw.getTranscripts().get(t).getGene());	
				vardraw.getTranscripts().get(t).getGene().mutations+=mutcount;
									
				VariantHandler.table.variants += mutcount;
				VariantHandler.table.revalidate();
				VariantHandler.table.repaint();	
			}
			vardraw.getTranscripts().get(t).getGene().intronic+=mutcount;	
			
			vardraw.inVarList = true;
			
								
			
			if(vardraw.inVarList) {
				if(!vardraw.getTranscripts().get(t).getGene().varnodes.contains(vardraw)) {
					vardraw.getTranscripts().get(t).getGene().varnodes.add(vardraw);
				}
			}
		}	
	}
		
	if(vardraw != null && vardraw.getExons() != null) {	
	
		for(int exon = 0; exon<vardraw.getExons().size(); exon++) {	
					
				amino = Main.chromDraw.getAminoChange(vardraw,base,vardraw.getExons().get(exon));
				
				if(amino.contains("UTR")) {
					if(VariantHandler.utr.isSelected()) {
						if(!VariantHandler.table.genearray.contains(vardraw.getExons().get(exon).getTranscript().getGene()) && vardraw.getExons().get(exon).getTranscript().getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
							
							VariantHandler.table.addEntry(vardraw.getExons().get(exon).getTranscript().getGene());
						}							
						vardraw.inVarList = true;
					
					//	vardraw.getExons().get(exon).getTranscript().getGene().mutations+=mutcount;
						vardraw.getExons().get(exon).getTranscript().getGene().utr +=mutcount;
						//VariantHandler.table.variants +=mutcount;
					}
					else {
						continue;
					}
				}
				/*if(!vardraw.coding) {
					
					if(VariantHandler.utr.isSelected()) {
						if(!VariantHandler.table.genearray.contains(vardraw.getExons().get(exon).getTranscript().getGene()) && vardraw.getExons().get(exon).getTranscript().getGene().samples.size() >= VariantHandler.geneSlider.getValue()) {
						
							VariantHandler.table.addEntry(vardraw.getExons().get(exon).getTranscript().getGene());
						}	
						vardraw.inVarList = true;
						
						vardraw.getExons().get(exon).getTranscript().getGene().mutations+=mutcount;
						vardraw.getExons().get(exon).getTranscript().getGene().utr+=mutcount;
						VariantHandler.table.variants+=mutcount;
					
						
						continue;
					}
					else {
						continue;
					}						
				}			
				/*if(amino.equals("")) {
					continue;
				}*/
				if(VariantHandler.nonsense.isSelected()) {
					if(!MethodLibrary.aminoEffect(amino).contains("nonsense")) {							
						continue;
					}						
				}
				if(VariantHandler.synonymous.isSelected()) {
					
					if(MethodLibrary.aminoEffect(amino).contains("synonymous")) {	
						
						continue;
					}							
				}					
				
				
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
						if(!vardraw.getExons().get(exon).getTranscript().getGene().samples.contains(entry.getValue().get(i).getSample())) {
							vardraw.getExons().get(exon).getTranscript().getGene().samples.add(entry.getValue().get(i).getSample());
							
						}
					}
					
					if(!VariantHandler.table.genearray.contains(vardraw.getExons().get(exon).getTranscript().getGene()) && vardraw.getExons().get(exon).getTranscript().getGene().samples.size() >= VariantHandler.geneSlider.getValue() ) {
						VariantHandler.table.addEntry(vardraw.getExons().get(exon).getTranscript().getGene());								
					}
				
					if(MethodLibrary.aminoEffect(amino).contains("nonsense")) {
						vardraw.getExons().get(exon).getTranscript().getGene().nonsense +=mutcount;
					}
					else if(MethodLibrary.aminoEffect(amino).contains("missense")) {
						
						vardraw.getExons().get(exon).getTranscript().getGene().missense  +=mutcount;
					}
					else if(MethodLibrary.aminoEffect(amino).contains("synonymous")) {
						
						vardraw.getExons().get(exon).getTranscript().getGene().synonymous +=mutcount;
					}
					
					
				}					
				vardraw.inVarList = true;
				if(!vardraw.getExons().get(exon).getTranscript().getGene().varnodes.contains(vardraw)) {			
					vardraw.getExons().get(exon).getTranscript().getGene().varnodes.add(vardraw);
					vardraw.getExons().get(exon).getTranscript().getGene().mutations+=mutcount;	
					VariantHandler.table.variants +=mutcount;
				}
				
				
				
				VariantHandler.table.revalidate();
				VariantHandler.table.repaint();
			}
			
		
	//	}			
		}
		
	}
	if(!vardraw.isInGene() && VariantHandler.intergenic.isSelected()) {
		
		VariantHandler.writeVariantToTSVFile(vardraw, output);
				
	}
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
	public static void main(String args[]) {
		/*String cram = "X:/cg7/Heikki/CRAMtest/cram3/Y_crc47_1_13-0818_cram3_test.cram";
		String crai = cram + ".crai";
		File ref = new File("C:/HY-Data/RKATAINE/BasePlayer/BasePlayer/genomes/Homo_sapiens_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa");
		*/
		
		try {
			String file = "X:/cg8/Riku/PlatinumGenomes/Ohutsuoli9_1_15-1064TP_T_160129_K00110_0062_AH55TWBBXX_1.snps_indels.ug.vcf.idx";
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
