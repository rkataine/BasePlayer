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
import java.awt.Color;
import java.awt.Composite;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.RandomAccessFile;
import java.util.ArrayList;

import htsjdk.tribble.readers.TabixReader;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map.Entry;

import javax.swing.JPanel;



public class ChromDraw extends JPanel implements MouseMotionListener, MouseListener {

private static final long serialVersionUID = 1L;
	int queryStart=0, queryEnd=0, memoryUsage = 0;
	Color[] colorPalette = new Color[11];
	VarNode varnode = null, vardraw;
	String s, startPhaseCodon, endPhaseCodon;

	int cytoHeight = 15, exonInfoWidth, aminoNro, codonStartPos, lastAmino;
	Runtime instance = Runtime.getRuntime();
	int mouseY=0, pressY=0, maxGeneLevel = 0;
	SplitClass splitClass = null;
	int middle = 0, baselevel, toMegabytes = 1048576;
	int exonFetchStart =0, exonFetchEnd = 0;
	static final Font middleFont = new Font("SansSerif", Font.BOLD, 20), seqFont=new Font("SansSerif", Font.BOLD, 12);

	int transIndex = 0, isoformIndex = 0;

	byte[] seqresult;
	int startPixel = 0, areaWidth = 0;
	static ArrayList<String[]> bandVector = new ArrayList<String[]>();
	ArrayList<String[]> chromBands = new ArrayList<String[]>();
	String resString, zoomtext;
	Transcript.Exon selectedExon = null;

	static HashMap<String, Integer> chromPos;
	ArrayList<Double> gerpList = new ArrayList<Double>();
	
	BufferedImage chromImage, splitCytoImage;
	Graphics2D chromImageBuffer;
	BufferedImage exonImage;
	Graphics2D exonImageBuffer;
	BufferedImage selectImage;
	Graphics2D selectImageBuffer;
	
	boolean zoomDrag = false, foundlevel = false, nameDraw = false;
	final static Color exonBarColor = new Color(20,100,20), forwardExon = new Color(50, 200, 100, 30), reverseExon = new Color(255, 50, 50, 30), seqpaint = new Color(100,200,200,150), highlight = new Color(255, 50, 50, 15);
	ArrayList<Integer> geneLevelMatrix= new ArrayList<Integer>();
	StringBuffer seqBuffer;
	RandomAccessFile chromo;
	public int zoompostemp;
	BufferedImage cytoImage;
	Graphics2D cytoImageBuffer;
	Rectangle mouseRect = new Rectangle();
	public int mouseX;//, tempDrag;
	private int bandwidth;
	private int Xpos;
	private String[] color;
	private String[] exonString = new String[8];
	static Color lightGray = new Color(230,230,230,100);
	
	private int level;

	private int screenPos, geneStartPos, geneEndPos, levelEndPos;
	private String gene;
	public String splitChrom = "";
	static Hashtable<String, String> aminoacids = new Hashtable<String, String>();
	static Hashtable<String, String> codons = new Hashtable<String, String>();
	static float[] dash = {5f};
	static TabixReader exonReader;
	
	ArrayList<String> exonRemove = new ArrayList<String>();
	
	private int exonwidth, genewidth;
	private FontMetrics fm;
	private Rectangle2D textWidth;	
	
	private int exonDrawY = 30;
	
	Composite backupe;
//	private int transcript.ypos;
	boolean lineZoomer;
	private int phase;
	
	private String base;
//	public boolean split = false;
	private BufferedImage tempImage;
	private Transcript transcript;
	private Transcript.Exon exon;
	private boolean found;
	private int aminopos;
	
	boolean updateExons = false;
	private int mutScreenPos;
	static Color backTransparent = new Color(255,255,230,255);
	private int baselength;
	private Transcript transcriptSelect;
	private Color varColor = new Color(100,100,100,30), gray = new Color(100,100,100,30);
	private boolean seqDrag;
	long timer;
	String message;
	private int mutcount;
	private Entry<String, ArrayList<SampleNode>> entry;
	private int prevExonEnd;
	private char[] array;
	private int nextExonStart;
	private boolean foundexon;
	private int foundexonindex;
	private Transcript.Exon clickedExon;
	private short foundcursor;	
	private boolean resizeSidebar;
	
	private SplitClass clickedSplit = null;
	private int seqstart;
	private int seqend;
	
	ChromDraw(int width, int height) {		
	//	super();		
		
		chromImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
		chromImageBuffer = (Graphics2D)chromImage.getGraphics();
		exonImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
		exonImageBuffer = (Graphics2D)exonImage.getGraphics();
		selectImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));
		selectImageBuffer = (Graphics2D)selectImage.getGraphics();		
		backupe = exonImageBuffer.getComposite();		
		addMouseMotionListener(this);
		addMouseListener(this);
		colorPalette[10] = new Color(255,60,60);
		colorPalette[9] = new Color(255,100,100);
		colorPalette[8] = new Color(255,140,140);
		colorPalette[7] = new Color(255,180,180);
		colorPalette[6] = new Color(255,220,220);
		colorPalette[5] = new Color(255,255,255);
		colorPalette[4] = new Color(220,220,255);
		colorPalette[3] = new Color(180,180,255);
		colorPalette[2] = new Color(140,140,255);
		colorPalette[1] = new Color(100,100,255);
		colorPalette[0] = new Color(80,80,255);
		aminoacids.put("ATT", "Ile");
		aminoacids.put("ATC", "Ile");
		aminoacids.put("ATA", "Ile");
	    
		aminoacids.put("CTT", "Leu");
		aminoacids.put("CTC", "Leu");
		aminoacids.put("CTA", "Leu");
		aminoacids.put("CTG", "Leu");
		aminoacids.put("TTA", "Leu");
		aminoacids.put("TTG", "Leu");
		
		aminoacids.put("GTT", "Val");
		aminoacids.put("GTC", "Val");
		aminoacids.put("GTA", "Val");
		aminoacids.put("GTG", "Val");
	    
		aminoacids.put("TTT", "Phe");
		aminoacids.put("TTC", "Phe");
		
		aminoacids.put("ATG", "Met");
	    
		aminoacids.put("TGT", "Cys");
		aminoacids.put("TGC", "Cys");
			 
		aminoacids.put("GCT", "Ala");
		aminoacids.put("GCC", "Ala");
		aminoacids.put("GCA", "Ala");
		aminoacids.put("GCG", "Ala");
		
		aminoacids.put("GGT", "Gly");
		aminoacids.put("GGC", "Gly");
		aminoacids.put("GGA", "Gly");
		aminoacids.put("GGG", "Gly");
		
		aminoacids.put("CCT", "Pro");
		aminoacids.put("CCC", "Pro");
		aminoacids.put("CCA", "Pro");
		aminoacids.put("CCG", "Pro");
		
		aminoacids.put("ACT", "Thr");
		aminoacids.put("ACC", "Thr");
		aminoacids.put("ACA", "Thr");
		aminoacids.put("ACG", "Thr");
		
		aminoacids.put("TCT", "Ser");
		aminoacids.put("TCC", "Ser");
		aminoacids.put("TCA", "Ser");
		aminoacids.put("TCG", "Ser");
		aminoacids.put("AGT", "Ser");
		aminoacids.put("AGC", "Ser");
		
		aminoacids.put("TAT", "Tyr");
		aminoacids.put("TAC", "Tyr");
		
		aminoacids.put("TGG", "Trp");
		
		aminoacids.put("CAA", "Gln");
		aminoacids.put("CAG", "Gln");
		
		aminoacids.put("AAT", "Asn");
		aminoacids.put("AAC", "Asn");
			 
		aminoacids.put("CAT", "His");
		aminoacids.put("CAC", "His");
			 
		aminoacids.put("GAA", "Glu");
		aminoacids.put("GAG", "Glu");
		
		aminoacids.put("GAT", "Asp");
		aminoacids.put("GAC", "Asp");
		
		aminoacids.put("AAA", "Lys");
		aminoacids.put("AAG", "Lys");
		
		aminoacids.put("CGT", "Arg");
		aminoacids.put("CGC", "Arg");
		aminoacids.put("CGA", "Arg");
		aminoacids.put("CGG", "Arg");
		aminoacids.put("AGA", "Arg");
		aminoacids.put("AGG", "Arg");
		
		aminoacids.put("TAA", "Stop");	 
		aminoacids.put("TAG", "Stop");	
		aminoacids.put("TGA", "Stop");
		
		codons.put("ATT", "ATT");
		codons.put("ATC", "ATC");
		codons.put("ATA", "ATA");
	    
		codons.put("CTT", "CTT");
		codons.put("CTC", "CTC");
		codons.put("CTA", "CTA");
		codons.put("CTG", "CTG");
		codons.put("TTA", "TTA");
		codons.put("TTG", "TTG");
		
		codons.put("GTT", "GTT");
		codons.put("GTC", "GTC");
		codons.put("GTA", "GTA");
		codons.put("GTG", "GTG");
	    
		codons.put("TTT", "TTT");
		codons.put("TTC", "TTC");
		
		codons.put("ATG", "ATG");
	    
		codons.put("TGT", "TGT");
		codons.put("TGC", "TGC");
			 
		codons.put("GCT", "GCT");
		codons.put("GCC", "GCC");
		codons.put("GCA", "GCA");
		codons.put("GCG", "GCG");
		
		codons.put("GGT", "GGT");
		codons.put("GGC", "GGC");
		codons.put("GGA", "GGA");
		codons.put("GGG", "GGG");
		
		codons.put("CCT", "CCT");
		codons.put("CCC", "CCC");
		codons.put("CCA", "CCA");
		codons.put("CCG", "CCG");
		
		codons.put("ACT", "ACT");
		codons.put("ACC", "ACC");
		codons.put("ACA", "ACA");
		codons.put("ACG", "ACG");
		
		codons.put("TCT", "TCT");
		codons.put("TCC", "TCC");
		codons.put("TCA", "TCA");
		codons.put("TCG", "TCG");
		codons.put("AGT", "AGT");
		codons.put("AGC", "AGC");
		
		codons.put("TAT", "TAT");
		codons.put("TAC", "TAC");
		
		codons.put("TGG", "TGG");
		
		codons.put("CAA", "CAA");
		codons.put("CAG", "CAG");
		
		codons.put("AAT", "AAT");
		codons.put("AAC", "AAC");
			 
		codons.put("CAT", "CAT");
		codons.put("CAC", "CAC");
			 
		codons.put("GAA", "GAA");
		codons.put("GAG", "GAG");
		
		codons.put("GAT", "GAT");
		codons.put("GAC", "GAC");
		
		codons.put("AAA", "AAA");
		codons.put("AAG", "AAG");
		
		codons.put("CGT", "CGT");
		codons.put("CGC", "CGC");
		codons.put("CGA", "CGA");
		codons.put("CGG", "CGG");
		codons.put("AGA", "AGA");
		codons.put("AGG", "AGG");
		
		codons.put("TAA", "TAA");	 
		codons.put("TAG", "TAG");	
		codons.put("TGA", "TGA");
		
	}
	void getDrawSeq(SplitClass split) {
		if(!ReferenceSeq.wait) {
			if(split.getReference() == null) {
				if((int)split.start-(int)split.viewLength < 0) seqstart = 0;				
				else seqstart = (int)split.start-(int)split.viewLength;
				if((int)split.end+(int)split.viewLength > split.chromEnd) seqend = split.chromEnd; 
				else seqend = (int)split.end+(int)split.viewLength;				
			
				split.setReference(new ReferenceSeq(split.chrom,seqstart, seqend, Main.referenceFile));
				
			}
			
			else if(split.getReference().getStartPos() > 1 && split.getReference().getStartPos() > split.start-200) {
				if(split.getReference().getStartPos()-(int)split.viewLength < 1) seqstart = 1;
				else seqstart = split.getReference().getStartPos()-(int)split.viewLength;
				split.getReference().appendToStart(seqstart);								
			}
			else if(split.getReference().getEndPos() < split.end+200) {
				if(split.getReference().getEndPos()+(int)split.viewLength > split.chromEnd) seqend = split.chromEnd;
				else seqend = split.getReference().getEndPos()+(int)split.viewLength+200;
				split.getReference().append(seqend);											
				
			}
		
		}
	}
	void drawSeq(SplitClass split) {
		
		if(split.viewLength <= Settings.readDrawDistance && split.viewLength > 10){
			
			split.getExonImageBuffer().setColor(Color.black);
			
			
			if(split.getReference() == null || split.getReference().getSeq() == null) {
				return;
			}
			if(Main.noreadthread) {
				
				return;
			}
			if(split.pixel >= 1) {
				split.getExonImageBuffer().fillRect(0, Main.chromScroll.getViewport().getHeight()-32, this.getWidth(), 15);
				
				for(int i = Math.max(0,(int)(split.start-split.getReference().getStartPos()-1)); i< split.getReference().getSeq().length; i++) {
					try {
				
					if(split.getReference().getStartPos()+i < split.start-1 ) {						
						continue;
					}
					else if(split.getReference().getStartPos()+i > split.end-1) {
						break;
					}
					else if((int)((split.getReference().getStartPos()+i -split.start)*split.pixel)+(int)(split.pixel/3) > (Main.drawCanvas.getDrawWidth())/2 && (int)((split.getReference().getStartPos()+i+1 -split.start)*split.pixel)+(int)(split.pixel/3+chromImageBuffer.getFont().getSize()) < (Main.drawCanvas.getDrawWidth())/2+split.pixel) {
						split.getExonImageBuffer().setColor(Color.black);
						split.getExonImageBuffer().fillRect((int)((split.getReference().getStartPos()+i -split.start)*split.pixel), Main.chromScroll.getViewport().getHeight()-40, (int)split.pixel, 20);
					
						split.getExonImageBuffer().setColor(Color.white);
						split.getExonImageBuffer().setFont(middleFont);
						
						split.getExonImageBuffer().drawString(Main.getBase.get(split.getReference().getSeq()[i]), (int)((split.getReference().getStartPos()+i -split.start)*split.pixel)+(int)(split.pixel/3), Main.chromScroll.getViewport().getHeight()-20);
					}
					else {	
						if(split.getReference().getSeq()[i] == (byte)'A') {
							split.getExonImageBuffer().setColor(Color.green);
						}
						else if(split.getReference().getSeq()[i] == (byte)'C') {
							split.getExonImageBuffer().setColor(Color.cyan);
						}
						else if(split.getReference().getSeq()[i] == (byte)'G') {
							split.getExonImageBuffer().setColor(Color.orange);
						}
						else if(split.getReference().getSeq()[i] == (byte)'T') {
							split.getExonImageBuffer().setColor(Color.red);
						}
						else {
							split.getExonImageBuffer().setColor(Color.gray);
						}
						split.getExonImageBuffer().setFont(seqFont);
						if(split.viewLength < 300) {
							split.getExonImageBuffer().drawString(Main.getBase.get(split.getReference().getSeq()[i]), (int)((split.getReference().getStartPos()+i+1 -split.start)*split.pixel)+(int)(split.pixel/3), Main.chromScroll.getViewport().getHeight()-20);
							
						}
						else {
							split.getExonImageBuffer().fillRect((int)((split.getReference().getStartPos()+i -split.start)*split.pixel), Main.chromScroll.getViewport().getHeight()-30, (int)split.pixel, 10);
						}
					}
					
					}
					catch(Exception e) {
						ErrorLog.addError(e.getStackTrace());
						e.printStackTrace();
					}
				}			
				// Middle line		
				split.getExonImageBuffer().setColor(Color.black);
				split.getExonImageBuffer().setStroke(Draw.dashed);			
				split.getExonImageBuffer().drawLine((Main.drawCanvas.getDrawWidth())/2, 0, ((Main.drawCanvas.getDrawWidth()))/2, Main.chromScroll.getViewport().getHeight());
				split.getExonImageBuffer().drawLine((int)((Main.drawCanvas.getDrawWidth())/2+split.pixel), 0, (int)((Main.drawCanvas.getDrawWidth())/2+split.pixel), Main.chromScroll.getViewport().getHeight());
				split.getExonImageBuffer().setStroke(Draw.basicStroke);
				split.getExonImageBuffer().drawString(MethodLibrary.formatNumber(getPosition((int)(Main.drawCanvas.getDrawWidth()/2.0 +split.pixel/2), split)), (int)((Main.drawCanvas.getDrawWidth())/2+split.pixel)+4, Main.chromScroll.getViewport().getHeight()-40);
			
			}
		}
		else {
			split.nullRef();
		
		}
	}

public int getPosition(int pixel, SplitClass split) {
	//return (int)(split.start+(pixel+split.pixel/2)/split.pixel);
	return (int)(split.start+(pixel/split.pixel));
}
	
void drawSideBar() {
	
	chromImageBuffer.setColor(Draw.sidecolor);

	chromImageBuffer.fillRect(0, 0, Main.sidebarWidth, Main.drawCanvas.splits.get(0).getExonImage().getHeight());
	chromImageBuffer.setColor(Color.black);
	
		
	
	if(Main.ref != null) {
		chromImageBuffer.drawString("Chromosome " +Main.chromosomeDropdown.getSelectedItem(), 10,15);
		if(Main.ref.getName().contains(".fa")) {
			chromImageBuffer.drawString("Ref: " +Main.ref.getName().substring(0, Main.ref.getName().indexOf(".fa")), 10,35);
		}
		else {
			chromImageBuffer.drawString("Ref: " +Main.ref.getName(), 10,35);
		}
	}
	else {
		chromImageBuffer.drawString("Add reference fasta file", 10,15);
	}
	chromImageBuffer.drawString("Genes: " +Main.annotationfile, 10,55);
	if(memoryUsage > (int)(((instance.totalMemory()-instance.freeMemory())/toMegabytes))) {
		if(instance.maxMemory()/toMegabytes  < 1000) {
			
			if((instance.maxMemory()-(instance.totalMemory()-instance.freeMemory()))/toMegabytes < 100) {		
				System.gc();
				if((instance.maxMemory()-(instance.totalMemory()-instance.freeMemory()))/toMegabytes < 100) {
					Main.drawCanvas.clearReads();
					System.gc();
					
					if((instance.maxMemory()-(instance.totalMemory()-instance.freeMemory()))/toMegabytes < 100) {
						Loader.memory = true;
					}
				}
			}
		}
		else if((instance.maxMemory()-(instance.totalMemory()-instance.freeMemory()))/toMegabytes < 300) {		
			System.gc();
			if((instance.maxMemory()-(instance.totalMemory()-instance.freeMemory()))/toMegabytes < 300) {		
				Main.drawCanvas.clearReads();
				System.gc();
				if((instance.maxMemory()-(instance.totalMemory()-instance.freeMemory()))/toMegabytes < 300) {
					Loader.memory = true;
				}
			}
			
		}
	}
	memoryUsage = (int)(((instance.totalMemory()-instance.freeMemory())/toMegabytes));
	if((instance.totalMemory()-instance.freeMemory())/(double)instance.maxMemory() > 0.8) {chromImageBuffer.setFont(seqFont);chromImageBuffer.setColor(Color.red);} else {chromImageBuffer.setFont(seqFont);chromImageBuffer.setColor(Color.black);}
	chromImageBuffer.drawString(""+memoryUsage +" / " +(instance.maxMemory()/toMegabytes) +"MB (" +(int)(MethodLibrary.round((instance.totalMemory()-instance.freeMemory())/(double)instance.maxMemory(),2)*100) +"%) " /*+FileRead.readcount*/, 10,Main.splitPane.getDividerLocation()-10);
	
}
void drawCyto(SplitClass split) {
	
	if(split.getCytoImage() == null) {
		split.setCytoImage(createBands(split.chrom));			
	}
	chromImageBuffer.setFont(Draw.defaultFont);
	chromImageBuffer.drawImage(split.getCytoImage(),split.offset,0,Main.drawCanvas.getDrawWidth(),cytoHeight,null);
	chromImageBuffer.setStroke(Draw.doubleStroke);
	chromImageBuffer.setColor(Color.black);
	chromImageBuffer.drawRect(split.offset, 1, Main.drawCanvas.getDrawWidth(), cytoHeight-2);
	chromImageBuffer.setColor(Color.red);	
	chromImageBuffer.setStroke(Draw.strongStroke);
	
	startPixel = (int)((Main.drawCanvas.getDrawWidth())*(split.start/(double)split.chromEnd))+split.offset;
	areaWidth = (int)((Main.drawCanvas.getDrawWidth())*split.viewLength/(double)split.chromEnd);
	if(areaWidth < 3) {
		chromImageBuffer.drawRect(startPixel+2, 2, areaWidth, cytoHeight-4);
	}
	else {
		chromImageBuffer.drawRect(startPixel+2, 2, areaWidth-2, cytoHeight-4);
	}
	chromImageBuffer.setColor(Color.black);	
	chromImageBuffer.setStroke(Draw.doubleStroke);
	
	for(int i =chromBands.size()-1 ; i>=0; i--) {
		try {
		if((Integer.parseInt(chromBands.get(i)[2])- Integer.parseInt(chromBands.get(i)[1]))*(Main.drawCanvas.getDrawWidth()/(double)chromPos.get(Main.refchrom +Main.drawCanvas.selectedSplit.chrom)) > 30) {
			if(chromBands.get(i)[4].contains(",0,") || chromBands.get(i)[4].equals("gpos100")) {
				if(chromImageBuffer.getColor() != Color.white) {
					chromImageBuffer.setColor(Color.white);
				}
				chromImageBuffer.drawString(chromBands.get(i)[3], (int)(Integer.parseInt(chromBands.get(i)[1])*((Main.drawCanvas.getWidth()-split.offset)/(double)chromPos.get(split.chrom)))+split.offset+4, 10);
				chromImageBuffer.setColor(Color.black);
			}
			else {
				chromImageBuffer.drawString(chromBands.get(i)[3], (int)(Integer.parseInt(chromBands.get(i)[1])*((Main.drawCanvas.getWidth()-split.offset)/(double)chromPos.get(split.chrom)))+split.offset+4, 10);
				
			}
			
		}
		if(!Main.drawCanvas.loading && mouseY < cytoHeight) {
			if(mouseX-split.offset > (int)(Integer.parseInt(chromBands.get(i)[1])*((Main.drawCanvas.getWidth()-split.offset)/(double)chromPos.get(split.chrom))) && mouseX-split.offset <= (int)(Integer.parseInt(chromBands.get(i)[2])*((Main.drawCanvas.getWidth()-split.offset)/(double)chromPos.get(split.chrom)))) {
				chromImageBuffer.setColor(Color.black);
				chromImageBuffer.fillRect(mouseX+18, 0, 30, cytoHeight);
				chromImageBuffer.setColor(Color.white);	
				chromImageBuffer.drawString(chromBands.get(i)[3], mouseX+20, 10);
				chromImageBuffer.setColor(Color.black);
			}
		}
		}
		catch(Exception e) {
			
		}
	}
	if(mouseY < cytoHeight) {
		
		if(getCursor().getType() != Cursor.TEXT_CURSOR) {
			
			setCursor(Cursor.getPredefinedCursor(Cursor.TEXT_CURSOR));
		}
		
	}
			
			
	chromImageBuffer.setFont(seqFont);
	
	
}
/*
void drawCytoSplit() {
	
	if(this.splitClass != null) {
	
		chromImageBuffer.setColor(Color.white);
		
		chromImageBuffer.fillRect((Main.drawCanvas.getDrawWidth())/2,0,(Main.drawCanvas.getDrawWidth()),this.getHeight());
		chromImageBuffer.setColor(Color.black);
		startPixel = (int)((Main.drawCanvas.getDrawWidth())/2*(this.splitClass.start/this.splitClass.chromEnd))+(Main.drawCanvas.getWidth())/2;
		areaWidth = (int)((Main.drawCanvas.getDrawWidth())/2*this.splitClass.viewLength/this.splitClass.chromEnd);
		chromImageBuffer.drawImage(this.splitClass.cytoImage, (Main.drawCanvas.getDrawWidth())/2,0,null);
		
		chromImageBuffer.setStroke(Draw.doubleStroke);
		
		chromImageBuffer.drawRect((Main.drawCanvas.getDrawWidth())/2, 1, (Main.drawCanvas.getDrawWidth())/2, cytoHeight);
		chromImageBuffer.setColor(Color.red);	
		chromImageBuffer.setStroke(Draw.strongStroke);
		chromImageBuffer.drawRect(startPixel, 2, areaWidth, cytoHeight-2);
		chromImageBuffer.setColor(Color.black);	
	}
}*/
void drawZoom() {	
	
	if(lineZoomer) {
		
		chromImageBuffer.setColor(Color.black);
		chromImageBuffer.setStroke(Draw.dashed);
		
		chromImageBuffer.drawLine(Main.drawCanvas.pressX, pressY, mouseX, mouseY);
		chromImageBuffer.setStroke(Draw.basicStroke);
		
	}
	else if(zoomDrag) {
		
		if(mouseY <= cytoHeight) {
			chromImageBuffer.setColor(Draw.zoomColor);
			chromImageBuffer.fillRect(Main.drawCanvas.pressX, 0, this.mouseX-Main.drawCanvas.pressX, cytoHeight);	
			
			
			return;
		}
		chromImageBuffer.setStroke(Draw.dashed);
		chromImageBuffer.setColor(Color.black);
		chromImageBuffer.setFont(seqFont);
		if(this.mouseX-Main.drawCanvas.pressX >= 0) {
			chromImageBuffer.drawRect(Main.drawCanvas.pressX, 0, this.mouseX-Main.drawCanvas.pressX, this.getHeight());	
			if(Main.drawCanvas.getDrawWidth()-this.mouseX > 200) {
				chromImageBuffer.drawString("" +MethodLibrary.formatNumber((int)((this.mouseX-Main.drawCanvas.pressX)/Main.drawCanvas.selectedSplit.pixel)) +"bp", Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX)+4, this.mouseY-35);
				chromImageBuffer.drawString("Right click to cancel zoom" ,Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX)+4, this.mouseY-6);
			}
			else {
				fm = exonImageBuffer.getFontMetrics();
				zoomtext = ""+MethodLibrary.formatNumber((int)((this.mouseX-Main.drawCanvas.pressX)/Main.drawCanvas.selectedSplit.pixel)) +"bp";
				textWidth = fm.getStringBounds(zoomtext, exonImageBuffer);
				this.zoompostemp = (int)((Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX))-textWidth.getWidth());
				chromImageBuffer.drawString(zoomtext, this.zoompostemp, this.mouseY-35);
				
				textWidth = fm.getStringBounds("Right click to cancel zoom", exonImageBuffer);
				this.zoompostemp = (int)((Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX))-textWidth.getWidth())-6;
				chromImageBuffer.drawString("Right click to cancel zoom" , this.zoompostemp, this.mouseY-6);
	
			}
			chromImageBuffer.setColor(Draw.zoomColor);
			chromImageBuffer.fillRect(Main.drawCanvas.pressX, 0, this.mouseX-Main.drawCanvas.pressX, this.getHeight());
			
		}
		else {
			
			lineZoomer = true;
			Main.drawCanvas.lineZoomer = true;
			zoomDrag = false;
		//	updateExons = true;
		
		//	Main.drawCanvas.repaint();
		
		}
		chromImageBuffer.setStroke(Draw.basicStroke);
	}
	
}
void drawScreen(Graphics g) {
	
//	chromImageBuffer.drawImage(Draw.image, Main.drawCanvas.selectedSplit.offset,0,Main.drawCanvas.getDrawWidth(),Main.chromScroll.getViewport().getHeight(), this);
//	chromImageBuffer.drawImage(Draw.image, 0,this.cytoHeight,Main.drawCanvas.getWidth(),Main.chromScroll.getViewport().getHeight(), this);
	
//	chromImageBuffer.setColor(new Color(255,255,230,240));
	
//	chromImageBuffer.fillRect(Main.drawCanvas.selectedSplit.offset,0,Main.drawCanvas.getDrawWidth(),this.getHeight());	
//	chromImageBuffer.fillRect(0,this.cytoHeight,Main.drawCanvas.getWidth(),this.getHeight());	
//	System.out.println(Main.drawCanvas.selectedSplit +" " +updateExons);
	
	
	drawSideBar();	
	drawCyto(Main.drawCanvas.selectedSplit);
	
	if(clickedExon != null) {
		drawClickedExon(clickedSplit);
	}
	
	drawExons(Main.drawCanvas.selectedSplit);
	for(int s = 0; s<Main.drawCanvas.splits.size(); s++) {
		
		chromImageBuffer.drawImage(Main.drawCanvas.splits.get(s).getExonImage(), Main.drawCanvas.splits.get(s).offset, this.cytoHeight,null);	
		if(Main.ref!=null) {
			chromImageBuffer.drawString("chr" +Main.drawCanvas.splits.get(s).chrom, Main.drawCanvas.splits.get(s).offset+4, this.cytoHeight+17);
		}
		else {
			chromImageBuffer.drawString("Add reference directory and fasta file to genomes -folder", Main.drawCanvas.splits.get(s).offset+4, this.cytoHeight+17);
			
		}
		chromImageBuffer.setStroke(Draw.strongStroke);
		chromImageBuffer.drawLine(Main.drawCanvas.splits.get(s).offset, 0, Main.drawCanvas.splits.get(s).offset, this.getHeight());
	}	
	
	chromImageBuffer.drawImage(selectImage, Main.drawCanvas.selectedSplit.offset, this.cytoHeight, null);

	if(VariantHandler.table != null && VariantHandler.table.hoverVar != null && VariantHandler.table.hoverVar.getExons() != null && VariantHandler.table.hoverVar.getExons().get(0).getTranscript().getChrom().equals(Main.chromosomeDropdown.getSelectedItem().toString())) {
		chromImageBuffer.setColor(Color.white);	
		chromImageBuffer.fillRect((int)((VariantHandler.table.hoverVar.getPosition()+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth-2, 17, 4, this.getHeight());
		chromImageBuffer.setColor(Color.black);
		chromImageBuffer.drawLine((int)((VariantHandler.table.hoverVar.getPosition()+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, 17, (int)((VariantHandler.table.hoverVar.getPosition()+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, this.getHeight());
		
	}
	if(VariantHandler.table != null && VariantHandler.table.hoverNode != null && VariantHandler.table.hoverNode.getChrom().equals(Main.chromosomeDropdown.getSelectedItem().toString())) {
		if(VariantHandler.table.hoverNode.getStart() > Main.drawCanvas.splits.get(0).start-1 && VariantHandler.table.hoverNode.getStart() < Main.drawCanvas.splits.get(0).end) {
			chromImageBuffer.setColor(Color.white);				
			chromImageBuffer.fillRect((int)((VariantHandler.table.hoverNode.getStart()+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth-2,17, 4, this.getHeight());
			chromImageBuffer.setColor(Color.black);
			chromImageBuffer.drawLine((int)((VariantHandler.table.hoverNode.getStart()+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth,17, (int)((VariantHandler.table.hoverNode.getStart()+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, this.getHeight());
		}
	}
	
	if(timer > 0) {
		chromImageBuffer.setColor(Color.black);
		if(System.currentTimeMillis()-timer < 4000) {
			chromImageBuffer.drawString(message, Main.sidebarWidth+10, Main.chromScroll.getViewport().getHeight()-40);			
		} 
		else {
			timer = 0;
		}		
	}
	if(mouseX > Main.sidebarWidth-3 && mouseX < Main.sidebarWidth+3) {
		
		if(getCursor().getType() != Cursor.W_RESIZE_CURSOR) {
			
			Main.drawCanvas.sampleZoomer = false;
			setCursor(Cursor.getPredefinedCursor(Cursor.W_RESIZE_CURSOR));
		}
	}
	else {
		if(getCursor().getType() == Cursor.W_RESIZE_CURSOR) {
			
			setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
		}
	}
	
	if(Main.searchChrom.equals(Main.drawCanvas.splits.get(0).chrom)) {
	
		drawSearchInterval();
	}
	if(this.seqDrag && mouseX-Main.drawCanvas.pressX > 0) {		
		chromImageBuffer.setColor(seqpaint);		
		chromImageBuffer.fillRect(Main.drawCanvas.pressX, Main.chromScroll.getViewport().getHeight()-17, mouseX-Main.drawCanvas.pressX, 15);
		chromImageBuffer.setColor(Color.black);
	}
	
	
	drawZoom();	
	g.drawImage(chromImage, 0, 0, null);
	
}	
void drawSearchInterval() {

	if(Main.searchStart > 0 && Main.searchEnd < 1) {
		
		if(Main.searchStart > Main.drawCanvas.splits.get(0).start && Main.searchStart < Main.drawCanvas.splits.get(0).end) {
			chromImageBuffer.setColor(Color.black);		
			chromImageBuffer.drawLine((int)((Main.searchStart-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, cytoHeight, (int)((Main.searchStart-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, this.getHeight());
			chromImageBuffer.drawLine((int)((Main.searchStart+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, cytoHeight, (int)((Main.searchStart+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, this.getHeight());
			chromImageBuffer.setColor(highlight);
			chromImageBuffer.fillRect((int)((Main.searchStart-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, cytoHeight, (int)Main.drawCanvas.splits.get(0).pixel, this.getHeight());
		}
	}
	else if(Main.searchStart > 0 && Main.searchEnd > 0) {
		chromImageBuffer.setColor(Color.black);	
		if(Main.searchStart > Main.drawCanvas.splits.get(0).start && Main.searchEnd < Main.drawCanvas.splits.get(0).end) {			
			chromImageBuffer.drawLine((int)((Main.searchStart-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, cytoHeight, (int)((Main.searchStart-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, this.getHeight());
			chromImageBuffer.drawLine((int)((Main.searchEnd+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, cytoHeight, (int)((Main.searchEnd+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, this.getHeight());
			chromImageBuffer.setColor(highlight);
			chromImageBuffer.fillRect((int)((Main.searchStart-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, cytoHeight, (int)((Main.searchEnd+1-Main.searchStart)*Main.drawCanvas.splits.get(0).pixel), this.getHeight());
			
		}
		else if(Main.searchStart > Main.drawCanvas.splits.get(0).start && Main.searchStart < Main.drawCanvas.splits.get(0).end &&  Main.searchEnd > Main.drawCanvas.splits.get(0).end) {
			chromImageBuffer.drawLine((int)((Main.searchStart-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, cytoHeight, (int)((Main.searchStart-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, this.getHeight());
			chromImageBuffer.setColor(highlight);
			chromImageBuffer.fillRect((int)((Main.searchStart-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, cytoHeight, (int)((Main.drawCanvas.splits.get(0).end-Main.searchStart)*Main.drawCanvas.splits.get(0).pixel), this.getHeight());
			
		}
		else if(Main.searchStart < Main.drawCanvas.splits.get(0).start && Main.searchEnd < Main.drawCanvas.splits.get(0).end && Main.searchEnd > Main.drawCanvas.splits.get(0).start){
			chromImageBuffer.drawLine((int)((Main.searchEnd+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, cytoHeight, (int)((Main.searchEnd+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel)+Main.sidebarWidth, this.getHeight());
			chromImageBuffer.setColor(highlight);
			chromImageBuffer.fillRect(Main.sidebarWidth, cytoHeight, (int)((Main.searchEnd-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel), this.getHeight());
			
		}
		else if(Main.searchStart < Main.drawCanvas.splits.get(0).start && Main.searchEnd > Main.drawCanvas.splits.get(0).end) {
			chromImageBuffer.setColor(highlight);
			chromImageBuffer.fillRect(Main.sidebarWidth, cytoHeight, (int)(Main.drawCanvas.splits.get(0).viewLength*Main.drawCanvas.splits.get(0).pixel), this.getHeight());
			
		}
		
		
		
	}
}
public void paint(Graphics g) {
	try {	
		
		drawScreen(g);
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
}

public StringBuffer getSeq(String chrom, int start, int end, RandomAccessFile seqchrom) {
	
	try {
		if(chrom == null) {
			return null;
		}
		if(start < 0) {
			start = 0;			
			if(end > Main.chromIndex.get(Main.refchrom + chrom)[1]) {
				end = Main.chromIndex.get(Main.refchrom +chrom)[1].intValue();
			}			
		}
		
		seqresult = new byte[(end-start+1)+((end-start)/(Main.chromIndex.get(Main.refchrom +chrom)[2].intValue()-1))];
		
		if(seqresult.length == 0 || seqresult.length > 200000) {
			return new StringBuffer();
		}
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
		return new StringBuffer();
	}
	
	seqBuffer = new StringBuffer();
	//chromo = seqchrom;
	
	try {
	
		try {
			seqchrom.seek((Main.chromIndex.get(Main.refchrom +chrom)[0]+(start)+((start)/Main.chromIndex.get(Main.refchrom +chrom)[2].intValue())));
		}
		catch(Exception e) {
			
			e.printStackTrace();
		}
		seqchrom.readFully(seqresult);
		
		if(seqresult[0] == 10) {			
			seqchrom.seek((Main.chromIndex.get(Main.refchrom +chrom)[0]+(start-1)+((start)/Main.chromIndex.get(Main.refchrom +chrom)[2].intValue())));
			seqchrom.readFully(seqresult);			
		}			
				
		for(int i= 0; i< seqresult.length; i++){
			
			if(seqresult[i] != 10) {
				
				seqBuffer.append((char)seqresult[i]);
				
			}
			
		
			if(seqBuffer.length() >= end-start) {
				break;
			}			
		}
		
	}
	catch (Exception e) {
		e.printStackTrace();
	
	}	
	
	
	return seqBuffer;			
}
static int[] makeLine(String[] split) {
	
	if(split.length > 16) {	   	
			
			int[] temp = { Integer.parseInt(split[8]), Integer.parseInt(split[9]), Integer.parseInt(split[12]), Integer.parseInt(split[13]), Integer.parseInt(split[14]), Integer.parseInt(split[15]), Integer.parseInt(split[16]), Integer.parseInt(split[17].trim()), Integer.parseInt(split[6]), Integer.parseInt(split[7])};
			return temp;
			}
	else {	      
			
			int[] temp = { Integer.parseInt(split[8]), Integer.parseInt(split[9]), Integer.parseInt(split[12]), Integer.parseInt(split[13]), Integer.parseInt(split[14]), Integer.parseInt(split[15].trim()), Integer.parseInt(split[6]), Integer.parseInt(split[7])};
			return temp;
	}	      
	
	
}
void drawExons(SplitClass split) {
	try {
	
	if(updateExons) {
		
		if(split.getGenes() == null) {
				split.getExonImageBuffer().setColor(backTransparent );			
				split.getExonImageBuffer().fillRect(0,0,Main.drawCanvas.getDrawWidth(),this.getHeight());  //(int)split.getExonImage().getHeight());
			//	return;
			}
	
		split.getExonImageBuffer().setFont(Draw.defaultFont);
		try {
	//		split.getExonImageBuffer().drawImage(Draw.image, 0,0,Main.drawCanvas.getDrawWidth(),(int)split.getExonImage().getHeight(), this);			
			split.getExonImageBuffer().setColor(backTransparent );			
			split.getExonImageBuffer().fillRect(0,0,Main.drawCanvas.getDrawWidth(),this.getHeight());  //(int)split.getExonImage().getHeight());	
						
			level = 1;
	
	if(split.getGenes() != null && split.getGenes().size() > 0) {
		
	
	
	if(split.equals(Main.drawCanvas.splits.get(0))) {
		drawMutations((exonDrawY*level)+2);
	}
	if(split.transStart > split.getGenes().size()-1) {
		split.transStart = 0;
	}
	if(split.getGenes().get(split.transStart) != null && split.getGenes().get(split.transStart).getEnd()+100000 < (int)split.start && split.transStart <split.getGenes().size()-2 && split.getGenes().get(split.transStart+1).getEnd()+100000 < (int)split.start) {
		while(split.getGenes().get(split.transStart) != null && split.getGenes().get(split.transStart).getEnd()+100000 < (int)split.start) {
			split.transStart++;	
			if(split.transStart >split.getGenes().size()-1 ) {
				break;
			}
		}		
	}
	else if(split.getGenes().get(split.transStart) != null && split.getGenes().get(split.transStart).getEnd()+100000 > (int)split.start) {
		while(split.transStart>0 && split.getGenes().get(split.transStart).getEnd()+100000 > (int)split.start) {
			split.transStart--;					
		}		
	}
	if(split.viewLength <= Settings.readDrawDistance && split.viewLength > 10){
		getDrawSeq(split);
	}
	for(int g = split.transStart; g < split.getGenes().size(); g++) {
		for(int i=0; i< split.getGenes().get(g).getTranscripts().size(); i++) {		
			
			transcript = split.getGenes().get(g).getTranscripts().get(i);
		
			if(!transcript.getGene().showIsoforms() && !transcript.equals(split.getGenes().get(g).getLongest())) {
				
				continue;
			}
			if(transcript.getStart() > (int)split.end) {
				break;
			}	
			geneStartPos = (int)((transcript.getStart()+1-split.start)*split.pixel);
			geneEndPos = (int)((transcript.getEnd()-split.start)*split.pixel);
			genewidth = (int)((transcript.getEnd() - (transcript.getStart()))*split.pixel);
			
		if(genewidth < 2 && split.viewLength > 10000000) {
				if(transcript.getStrand()) {
					split.getExonImageBuffer().setColor(exonBarColor);			
				}
				else {					
					split.getExonImageBuffer().setColor(Color.RED);
				}
				split.getExonImageBuffer().drawLine(geneStartPos, 5, geneStartPos, 2);
				continue;
			}
			
			gene = transcript.getGenename();
			fm = split.getExonImageBuffer().getFontMetrics();		
			textWidth = fm.getStringBounds(gene, split.getExonImageBuffer());
			
			if(geneEndPos < geneStartPos+textWidth.getWidth()) {
				levelEndPos = (int)(geneStartPos+textWidth.getWidth())+5;
			}
			else {
				levelEndPos = geneEndPos+5;
			}
						
				if(geneLevelMatrix.isEmpty()) {
				
					geneLevelMatrix.add(levelEndPos);
					level = 1;
				}
				else {
					
					foundlevel = false;
					for(int c=0;c<geneLevelMatrix.size(); c++) {
						
						if(geneLevelMatrix.get(c) < geneStartPos) {
							
							level = c+1;
							foundlevel = true;				
							geneLevelMatrix.set(c, levelEndPos);
							break;
						}
					
					}
				if(!foundlevel) {			
						geneLevelMatrix.add(levelEndPos);
						level = geneLevelMatrix.size();
						
					}
				}
				transcript.ypos =  (exonDrawY*level)+2;
				if(genewidth > 2) {
					split.getExonImageBuffer().setColor(Color.black);
					if(geneEndPos-geneStartPos >0) {
						split.getExonImageBuffer().drawLine(geneStartPos, transcript.ypos+5, geneEndPos, transcript.ypos+5);
					}
				}
				nameDraw = false;
			
			
			for (int e = 0; e<transcript.getExons().length; e++) { 
					try {
						if(e > transcript.getExons().length-1) {
							break;
						}
						exon = transcript.getExons()[e];
					}
					catch(Exception ex) {
						
						ex.printStackTrace();
					}
					if (exon.getEnd() < split.start) {
						continue;
					}
					else if ( exon.getStart() > split.end) {
						break;
					}
					
					screenPos = (int)((exon.getStart()+1-split.start)*split.pixel);						
					exonwidth = (int)((exon.getEnd() - (exon.getStart()))*split.pixel)+1;				
					exon.getRectangle().setBounds(screenPos, transcript.ypos, exonwidth, exonDrawY/2);				
				
					if(exonwidth < 1) {
						exonwidth = 1;
					}				
				
					split.getExonImageBuffer().setColor(Color.black);	
					
					if(!nameDraw && screenPos >= 0) {
						split.getExonImageBuffer().drawString(gene, screenPos, transcript.ypos-1);
						nameDraw = true;
					}
					
					if(transcript.getStrand()) {
						split.getExonImageBuffer().setColor(exonBarColor);			
					}
					else {					
						split.getExonImageBuffer().setColor(Color.RED);
					}
					if(exon.getStartPhase() == -1) {
						split.getExonImageBuffer().setColor(Color.lightGray);
					}				
					
					split.getExonImageBuffer().fillRect(screenPos, transcript.ypos, exonwidth, exonDrawY/2);
					//UTR draw
					
					exonwidth = (int)((transcript.getCodingStart()-exon.getStart())*split.pixel)+1;
					
					//Start & End exon
					if(exon.getEnd() > transcript.getCodingEnd() && exon.getStart() < transcript.getCodingEnd() && exon.getEnd() > transcript.getCodingStart() && exon.getStart() < transcript.getCodingStart()) {
						split.getExonImageBuffer().setColor(Color.lightGray);
						
						exonwidth = (int)((transcript.getCodingStart()-exon.getStart())*split.pixel)+1;					
						split.getExonImageBuffer().fillRect(screenPos, transcript.ypos, exonwidth, exonDrawY/2);			
						
						screenPos = (int)((transcript.getCodingEnd()+1-split.start)*split.pixel)+1;
						exonwidth = (int)((exon.getEnd()-transcript.getCodingEnd())*split.pixel)+1;
						
						split.getExonImageBuffer().fillRect(screenPos, transcript.ypos, exonwidth, exonDrawY/2);
						
					}
					//Start exon
					else if(exon.getStart() < transcript.getCodingStart() && exon.getEnd() > transcript.getCodingStart()) {
						split.getExonImageBuffer().setColor(Color.lightGray);					
						
						exonwidth = (int)((transcript.getCodingStart()-exon.getStart())*split.pixel)+1;
						
						split.getExonImageBuffer().fillRect(screenPos, transcript.ypos, exonwidth, exonDrawY/2);
						
					}
					//End exon
					else if (exon.getStart() < transcript.getCodingEnd() && exon.getEnd() > transcript.getCodingEnd()) {
						split.getExonImageBuffer().setColor(Color.lightGray);
						
						exonwidth = (int)((exon.getEnd()-transcript.getCodingEnd())*split.pixel)+1;
						screenPos = (int)((transcript.getCodingEnd()+1-split.start)*split.pixel)+1;					
						split.getExonImageBuffer().fillRect(screenPos, transcript.ypos,exonwidth, exonDrawY/2);
				
					}
					
					if(split.viewLength < 200) {
						if(split.getReference() == null) {
							
							break;
						}
						// FORWARD STRAND aminoacid sequence draw TODO split
						
						if(transcript.getStrand()) {
							
							if(split.end > transcript.getCodingStart() && split.start < transcript.getCodingEnd()) {
								if(exon.getStartPhase() == -1){
									continue;
								}
								// Start exon
								else if(exon.getEnd() > transcript.getCodingStart() && exon.getStart() < transcript.getCodingStart()) {
									aminoNro = exon.getFirstAmino();
									
									for(int codon = transcript.getCodingStart(); codon<exon.getEnd()-exon.getEndPhase(); codon+=3) {
										if(codon >= transcript.getCodingEnd() || (codon-split.getReference().getStartPos())+3 > split.getReference().getSeq().length-1 || codon-split.getReference().getStartPos() < 0) {
											break;
										}
																	
										codonStartPos = (int)((codon+1-split.start)*split.pixel);
										split.getExonImageBuffer().setColor(Color.white);								
										split.getExonImageBuffer().drawLine(codonStartPos, transcript.ypos, codonStartPos, transcript.ypos+this.exonDrawY/2-2);
										//MethodLibrary.getAminoAcid(MethodLibrary.getAminoAcid(new String(Arrays.copyOfRange(split.getReference().getSeq(), codon-split.getReference().getStartPos(),(codon-split.getReference().getStartPos())+3)))); 
										
										//if(split.sequence != null) {
											split.getExonImageBuffer().drawString(MethodLibrary.getAminoAcid(new String(Arrays.copyOfRange(split.getReference().getSeq(), codon-split.getReference().getStartPos(),(codon-split.getReference().getStartPos())+3))), codonStartPos, transcript.ypos+this.exonDrawY/2-2);
											split.getExonImageBuffer().drawString("" +aminoNro, codonStartPos+(int)(3*split.pixel)/2, transcript.ypos+this.exonDrawY/2-2);								
										//}
										aminoNro++;
									
									}
									split.getExonImageBuffer().drawLine((int)((exon.getEnd()-exon.getEndPhase()+1-split.start)*split.pixel), transcript.ypos, (int)((exon.getEnd()-exon.getEndPhase()+1-split.start)*split.pixel), transcript.ypos+this.exonDrawY/2-2);
									
								}
								else {
									// Other exons
									aminoNro = exon.getFirstAmino();
									for(int codon = exon.getStart()+exon.getStartPhase(); codon<exon.getEnd()-exon.getEndPhase(); codon+=3) {
										if(codon >= transcript.getCodingEnd() || codon > split.end) {
											break;
										}
										
										codonStartPos = (int)((codon+1-split.start)*split.pixel);
										split.getExonImageBuffer().setColor(Color.white);
										split.getExonImageBuffer().drawLine(codonStartPos, transcript.ypos, codonStartPos, transcript.ypos+(this.exonDrawY/2-2));
										try {
										if(codon-split.getReference().getStartPos() > 0) {
											split.getExonImageBuffer().drawString(MethodLibrary.getAminoAcid(new String(Arrays.copyOfRange(split.getReference().getSeq(), codon-split.getReference().getStartPos(),(codon-split.getReference().getStartPos())+3))), codonStartPos, transcript.ypos+this.exonDrawY/2-2);
											}
										}
										catch(Exception ex) {
											ex.printStackTrace();
										}
										split.getExonImageBuffer().drawString("" +aminoNro, codonStartPos+(int)(3*split.pixel)/2, transcript.ypos+this.exonDrawY/2-2);
										aminoNro++;
									}
									split.getExonImageBuffer().drawLine((int)((exon.getEnd()-exon.getEndPhase()+1-split.start)*split.pixel), transcript.ypos, (int)((exon.getEnd()-exon.getEndPhase()+1-split.start)*split.pixel), transcript.ypos+this.exonDrawY/2-2);
									
								}
							}					
							
						}
						else {
							//REVERSE STRAND
							
							if(split.end > transcript.getCodingStart() && split.start < transcript.getCodingEnd()) {
								if(exon.getStartPhase() == -1) {
									continue;
								}
								
								// Start exon
								else if(exon.getEnd() > transcript.getCodingEnd() && exon.getStart() < transcript.getCodingEnd()) {
									aminoNro = 1;
									
									for(int codon = transcript.getCodingEnd(); codon-3>=exon.getStart(); codon-=3) {
									//for(int codon = transcript.getCodingEnd(); codon-3>=exon.getStart(); codon-=3) {	
										
										if(codon < transcript.getCodingStart()+1 ) {
											break;
										}
										
										if(codon-3-split.getReference().getStartPos() < 0|| codon-split.getReference().getStartPos()>split.getReference().getSeq().length-1) {
											continue;
										}
											
										codonStartPos = (int)((codon-2-split.start)*split.pixel);
										split.getExonImageBuffer().setColor(Color.white);								
										split.getExonImageBuffer().drawLine(codonStartPos, transcript.ypos, codonStartPos, transcript.ypos+(this.exonDrawY/2-2));
										
										split.getExonImageBuffer().drawString(MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(Arrays.copyOfRange(split.getReference().getSeq(), codon-3-split.getReference().getStartPos(),codon-split.getReference().getStartPos())))), codonStartPos, transcript.ypos+this.exonDrawY/2-2);
										
										split.getExonImageBuffer().drawString("" +aminoNro, codonStartPos+(int)(3*split.pixel)/2, transcript.ypos+this.exonDrawY/2-2);								
										aminoNro++;								
									}
									split.getExonImageBuffer().drawLine((int)((exon.getEnd()+1-exon.getStartPhase()-split.start)*split.pixel), transcript.ypos, (int)((exon.getEnd()+1-exon.getStartPhase()-split.start)*split.pixel), transcript.ypos+this.exonDrawY);
									
								}
								// Other exons
								else{
									
									
									aminoNro = exon.getFirstAmino();
									if(split.getReference() == null) {
										continue;
									}
									for(int codon = exon.getEnd()-exon.getStartPhase(); codon-3>=exon.getStart(); codon-=3) {
										if(codon-3 < transcript.getCodingStart()) {
											break;
										}
										codonStartPos = (int)((codon-2-split.start)*split.pixel);
										split.getExonImageBuffer().setColor(Color.white);
										split.getExonImageBuffer().drawLine(codonStartPos, transcript.ypos, codonStartPos, transcript.ypos+(this.exonDrawY/2));
										try {
											if(codon > split.start && codon < split.end) {
												
												split.getExonImageBuffer().drawString(MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(Arrays.copyOfRange(split.getReference().getSeq(), codon-3-split.getReference().getStartPos(),codon-split.getReference().getStartPos())))), codonStartPos, transcript.ypos+this.exonDrawY/2-2);
												split.getExonImageBuffer().drawString("" +aminoNro, codonStartPos+(int)(3*split.pixel)/2, transcript.ypos+this.exonDrawY/2-2);
											}
											
										}
										catch(Exception ex) {
											ex.printStackTrace();
										}
										aminoNro++;
									}
									split.getExonImageBuffer().drawLine((int)((exon.getEnd()+1-exon.getStartPhase()-split.start)*split.pixel), transcript.ypos, (int)((exon.getEnd()+1-exon.getStartPhase()-split.start)*split.pixel), transcript.ypos+this.exonDrawY/2);
									
								}
							}
						}
					}
					
			}	
			
		}
	}
	//exon for end
	
	this.setPreferredSize(new Dimension(Main.drawScroll.getViewport().getWidth(), ((geneLevelMatrix.size()+1) * exonDrawY)));
	this.revalidate();		
	
	geneLevelMatrix.clear();	
	
	updateExons = false;
	}
	if(Main.referenceFile != null) {	
		drawSeq(split);
	}
	}
	catch(Exception e) {
		e.printStackTrace();
		
	}
	}
	else {
		
		selectImageBuffer.setComposite( Main.drawCanvas.composite);		
		selectImageBuffer.fillRect(0,0,Main.drawCanvas.getWidth(),this.getHeight());					
		selectImageBuffer.setComposite(backupe);	
		if(split.transStart > split.getGenes().size()-1) {
			split.transStart = 0;
		}
		if(Main.drawCanvas.loading) {
			return;
		}
		for(int g = split.transStart; g<split.getGenes().size(); g++ ) {
			for(int i=0; i< split.getGenes().get(g).getTranscripts().size(); i++) {
				transcriptSelect = split.getGenes().get(g).getTranscripts().get(i);
				
				
				if(!transcriptSelect.getGene().showIsoforms() && !transcriptSelect.equals(split.getGenes().get(g).getLongest())) {
					
					continue;
				}
				
				if(transcriptSelect.getStart() > (int)split.end) {
					break;
				}	
				if(transcriptSelect != null && transcriptSelect.getEnd()+100000 < (int)split.start) {
					continue;
				}
				if(transcriptSelect.getStart() > (int)split.end) {
					break;
				}
				if((transcriptSelect.getEnd()-transcriptSelect.getStart())*split.pixel < 2 && split.viewLength > 10000000) {				
					continue;
				}
				
			for (int e = 0; e<transcriptSelect.getExons().length; e++) { 	
				exon = transcriptSelect.getExons()[e];
				if (exon.getEnd() < split.start) {
					continue;
				}
				else if ( exon.getStart() > split.end) {
					break;
				}		
		
			if(selectedExon == null && exon.getRectangle().intersects(mouseRect)){					
					
					selectedExon = exon;					
					transIndex = i;
					//TODO Splice Amino
					/*if(viewLength < 200 && (exon.getStartPhase() > 0 || exon.getEndPhase() > 0)) {
						if(exon.getStartPhase() > 0) {
							this.startPhaseCodon = MethodLibrary.getSpliceAmino(Main.chromosomeDropdown.getItemAt(Main.selectedChrom), exon, true);
						}
						if(exon.getEndPhase() > 0) {
							this.endPhaseCodon = MethodLibrary.getSpliceAmino(Main.chromosomeDropdown.getItemAt(Main.selectedChrom),exon, false);
						}						
					}*/				
				}
				if(selectedExon != null && exon == selectedExon) {
					
					screenPos = (int)((exon.getStart()+1-split.start)*split.pixel);									
					exonwidth = (int)((exon.getEnd() - (exon.getStart()))*split.pixel)+1;		
					transcript.ypos = (int)selectedExon.getRectangle().getY();
					selectImageBuffer.setColor(Color.red);					
					selectImageBuffer.drawRect(screenPos-1, transcript.ypos-1, exonwidth+1, exonDrawY/2+1);	
					
					//Splice amino draw
				
					if(split.viewLength < 200 && (exon.getStartPhase() > 0 || exon.getEndPhase() > 0)) {
						
						selectImageBuffer.setColor(Color.black);	
						//Forward strand
						
						if(transcriptSelect.getStrand()) {
							if(exon.getStartPhase() > 0) {
								selectImageBuffer.drawString(startPhaseCodon +" " +(exon.getFirstAmino()-1), (int)(screenPos-((3-exon.getStartPhase())*split.pixel)), transcript.ypos+exonDrawY);
							}
							if(exon.getEndPhase() > 0) {
								
								if(exon.getFirstAmino() == 1) {
									lastAmino = ((exon.getEnd()-transcriptSelect.getCodingStart())/3)+1;
								}
								else {
									lastAmino = ((exon.getFirstAmino()+(exon.getEnd()-(exon.getStart()+exon.getStartPhase()))/3));
								}
						//		selectImageBuffer.drawString(endPhaseCodon +" " +lastAmino, (int)((exon.getEnd()+1-split.start)*split.pixel)+4, transcript.ypos+exonDrawY/2+4);
							}
						}
						//Reverse strand
						else {
							
							if(exon.getStartPhase() > 0) {
								selectImageBuffer.drawString(startPhaseCodon +" " +(exon.getFirstAmino()-1), (int)((exon.getEnd()+1-split.start)*split.pixel)+4, transcript.ypos+exonDrawY/2+4);
								
							}
							if(exon.getEndPhase() > 0) {
								if(exon.getFirstAmino() == 1) {
									lastAmino = ((transcriptSelect.getCodingEnd() -exon.getStart())/3)+1;
								}
								else {
									lastAmino = ((exon.getFirstAmino()+(exon.getEnd()-(exon.getStart()+exon.getStartPhase()))/3));
								}
						//		selectImageBuffer.drawString(endPhaseCodon +" " +lastAmino, (int)(screenPos-((3-exon.getEndPhase())*split.pixel)), transcript.ypos+exonDrawY);
								
								}
						}					
					}		
				}
				}
			}
		}
		//select exon end
		
		if(selectedExon != null && !selectedExon.equals(clickedExon)) {
			if(getCursor().getType() != Cursor.HAND_CURSOR) {
				setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
			}
			screenPos = this.mouseRect.x;
			transcript.ypos = (int)selectedExon.getRectangle().getY();
			exonString[0] = "Exon "+selectedExon.getNro();
			exonString[1] = selectedExon.transcript.getGenename();
			exonString[2] = selectedExon.transcript.getENSG();
			
			
			exonString[3] = selectedExon.transcript.getENST();
					
			
			exonString[4] = selectedExon.transcript.getBiotype();
			exonString[5] = selectedExon.transcript.getUniprot();
			exonString[6] = selectedExon.transcript.getDescription();	
			if(selectedExon.transcript.getGene().getTranscripts().size() > 1) {
				if(selectedExon != null && selectedExon.transcript.getGene().showIsoforms()) {
					exonString[7] = "Right click to collapse isoforms.";
				}
				else {
					exonString[7] = "Right click to expand all isoforms.";
				}
			}
			else {
				exonString[7] = "";
			}
			fm = selectImageBuffer.getFontMetrics(); 
			exonInfoWidth = (int)fm.getStringBounds(exonString[maxWidth(exonString)], exonImageBuffer).getWidth()+2;
			selectImageBuffer.setColor(Color.white);
			
			selectImageBuffer.fillRect(screenPos, transcript.ypos+exonDrawY, exonInfoWidth, (int)(textWidth.getHeight()*(exonString.length))+2);
			selectImageBuffer.setColor(Color.black);
			selectImageBuffer.drawRect(screenPos-1, transcript.ypos+exonDrawY-1, exonInfoWidth+2, (int)(textWidth.getHeight()*(exonString.length))+4);
			
			for(int i=0;i<exonString.length;i++) {
				selectImageBuffer.drawString(exonString[i], screenPos+2, (int)(transcript.ypos+exonDrawY+(textWidth.getHeight())*(i+1)));
			}
			
		}
		if(selectedExon != null && !selectedExon.getRectangle().intersects(mouseRect) && mouseY > cytoHeight) {
			
			selectedExon = null;
			
			if(getCursor().getType() != Cursor.DEFAULT_CURSOR) {
				setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
			}
		}	
	}
	
	
	}
	catch(Exception ex) {
		
		ex.printStackTrace();
		return;
	}
	
}
void drawClickedExon(SplitClass split) {
	if(clickedExon != null) {
		
		screenPos = (int)clickedExon.getRectangle().getX();
		if(screenPos < 12) {
			screenPos = 12;
		}
		
		transcript.ypos = (int)clickedExon.getRectangle().getY();
		exonString[0] = "Exon "+clickedExon.getNro();
		exonString[1] = clickedExon.transcript.getGenename();
		if(clickedExon.transcript.getENSG().contains("ENS")) {
			exonString[2] = clickedExon.transcript.getENSG() +" (view in Ensembl)";
		}
		else if(clickedExon.transcript.getENSG().startsWith("GeneID")) {
			exonString[2] = clickedExon.transcript.getENSG() +" (view in RefSeq)";
		}
		else {
			exonString[2] = clickedExon.transcript.getENSG();
		}
		if(clickedExon.transcript.getENST().contains("ENS")) {
			exonString[3] = clickedExon.transcript.getENST() +" (view in Ensembl)";
				
		}
		exonString[4] = clickedExon.transcript.getBiotype();
		exonString[5] = "View in GeneCards";
		exonString[6] = clickedExon.transcript.getDescription();	
		if(clickedExon.transcript.getGene().getTranscripts().size() > 1) {
			if(clickedExon != null && clickedExon.transcript.getGene().showIsoforms()) {
				exonString[7] = "Right click to collapse isoforms.";
			}
			else {
				exonString[7] = "Right click to expand all isoforms.";
			}
		}
		else {
			exonString[7] = "";
		}
		fm = split.getExonImageBuffer().getFontMetrics(); 
		exonInfoWidth = (int)fm.getStringBounds(exonString[maxWidth(exonString)], exonImageBuffer).getWidth()+2;
		split.getExonImageBuffer().setColor(Color.white);
		
		split.getExonImageBuffer().fillRect(screenPos, transcript.ypos+exonDrawY, exonInfoWidth, (int)(textWidth.getHeight()*(exonString.length))+2);
		split.getExonImageBuffer().setColor(Color.black);
		split.getExonImageBuffer().drawRect(screenPos-1, transcript.ypos+exonDrawY-1, exonInfoWidth+2, (int)(textWidth.getHeight()*(exonString.length))+4);
		foundcursor = -1;
		for(int i=0;i<exonString.length-1;i++) {
			if(mouseRect.getX() > screenPos && mouseRect.getX() < screenPos +exonInfoWidth) {
				if(i<6 && mouseRect.getY() > (int)(transcript.ypos+exonDrawY+(textWidth.getHeight())*(i)) && mouseRect.getY() < (int)(transcript.ypos+exonDrawY+(textWidth.getHeight())*(i+1))) {
					split.getExonImageBuffer().setFont(seqFont);
					foundcursor = (short)i;
					if(getCursor().getType() != Cursor.HAND_CURSOR ) {
						setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
						
					}
				}
				else {
					split.getExonImageBuffer().setFont(Draw.defaultFont);
				}
			}
			else {
				split.getExonImageBuffer().setFont(Draw.defaultFont);
				
			}
			
			if(split.getExonImageBuffer().getFont().equals(seqFont)) {
				if(i == 0) {
					split.getExonImageBuffer().drawString(exonString[i] +"  (zoom to exon)", screenPos+2, (int)(transcript.ypos+exonDrawY+(textWidth.getHeight())*(i+1)));
				}
				else if(i==1) {
					split.getExonImageBuffer().drawString(exonString[i] +"  (zoom to gene)", screenPos+2, (int)(transcript.ypos+exonDrawY+(textWidth.getHeight())*(i+1)));
				}
				else if(i==2 || i == 3) {
					split.getExonImageBuffer().drawString(exonString[i], screenPos+2, (int)(transcript.ypos+exonDrawY+(textWidth.getHeight())*(i+1)));
				}
				else if(i>=4) {
					split.getExonImageBuffer().drawString(exonString[i] , screenPos+2, (int)(transcript.ypos+exonDrawY+(textWidth.getHeight())*(i+1)));
				}
			}
			else {
				split.getExonImageBuffer().drawString(exonString[i], screenPos+2, (int)(transcript.ypos+exonDrawY+(textWidth.getHeight())*(i+1)));
			}
		}
		if(foundcursor < 0 && selectedExon == null && mouseY > cytoHeight && mouseY < Main.chromScroll.getViewport().getHeight() -15) {
			
			setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
		}
		
	}
}

void drawMutations(int ylevel, Transcript trans) {
	try {
	
	if(varnode != null && Main.drawCanvas.selectedSplit.viewLength < 1000000) {
		
		vardraw = varnode;
		
		Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(varColor);
		
		while(vardraw != null && vardraw.getPosition() < Main.drawCanvas.splits.get(0).end+1) {		
			
			if(Main.drawCanvas.hideNode(vardraw)) {
				vardraw = vardraw.getNext();
				continue;
			}
		
			found = false;
			if(vardraw.getPrev().getPosition()-vardraw.getPosition() != -1) {
				baselevel = 1;
			}
			for(int j = 0; j<vardraw.vars.size(); j++) { //Map.Entry<String, ArrayList<SampleNode>> entry : vardraw.vars) {
				entry = vardraw.vars.get(j);
				if(Main.drawCanvas.hideNodeVar(vardraw, entry)) {
					continue;
				}
				
				for(int i = 0; i<entry.getValue().size(); i++) {					
					if(!Main.drawCanvas.hideVar(entry.getValue().get(i))) {
						found = true;						
						break;
					}
				}
				if(!found) {
				
					continue;
				}
				 if(vardraw.indel && entry.getKey().length() > 1) {
					 mutScreenPos = (int)((vardraw.getPosition()+2-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel);
				 }
				 else {
					 mutScreenPos = (int)((vardraw.getPosition()+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel);
				 }
				if(vardraw.getExons().size() > 0) {
					foundexon = false;
					
					for(int i = 0; i< vardraw.getExons().size(); i++) {
						if(!vardraw.getExons().get(i).getTranscript().equals(trans)) {
							return;
						}
							
						if(vardraw.indel && entry.getKey().length() > 1) {
							if(vardraw.getPosition()+1 > vardraw.getExons().get(i).getTranscript().getCodingStart() && vardraw.getPosition()+1 < vardraw.getExons().get(i).getTranscript().getCodingEnd()) {
								if(vardraw.getPosition()+1 >= vardraw.getExons().get(i).getStart()-2 && vardraw.getPosition()+1 < vardraw.getExons().get(i).getEnd()+2) {
									foundexon = true;
									foundexonindex = i;
									
									break;
								}
							}							
						}
						else {
							if(vardraw.getPosition() >= vardraw.getExons().get(i).getTranscript().getCodingStart() && vardraw.getPosition() < vardraw.getExons().get(i).getTranscript().getCodingEnd()) {
								foundexon = true;
								foundexonindex = i;
								break;
							}
						}
					}
					
					if(foundexon) {
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(reverseExon);
						
					}
					else {
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(varColor);
					}		
					
					Main.drawCanvas.splits.get(0).getExonImageBuffer().fillRect(mutScreenPos, 0, (int)Main.drawCanvas.splits.get(0).pixel+1, this.getHeight());
					Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(varColor);
				}
				else {					
					Main.drawCanvas.splits.get(0).getExonImageBuffer().fillRect(mutScreenPos, 0, (int)Main.drawCanvas.splits.get(0).pixel+1, this.getHeight());					
				}
				
				if( Main.drawCanvas.splits.get(0).viewLength < 200 && vardraw.getExons().size() > 0 && trans != null && vardraw.getExons().get(0).transcript.equals(trans)) {
					 //TODO draw aminochanges	 				
				
					mutcount = 0;
					
					for(int m = 0; m<entry.getValue().size(); m++) {
						if(!Main.drawCanvas.hideVar(entry.getValue().get(m))) {
							mutcount++;
						}
					}
					if(mutcount == 0) {
						break;
					}					
										
						base = entry.getKey();						
							
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setFont(seqFont);
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(Color.black);		
						
						Main.drawCanvas.splits.get(0).getExonImageBuffer().drawString(mutcount +" " +getAminoChange(vardraw,base,vardraw.getExons().get(foundexonindex)), mutScreenPos, ylevel+(30*baselevel));
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(varColor);
						baselevel++;
						
					}
					Main.drawCanvas.splits.get(0).getExonImageBuffer().setFont(Draw.defaultFont);					
				
				}
				if(vardraw != null) {
					vardraw = vardraw.getNext();
				
				}		
		}		
		
	}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	Main.drawCanvas.splits.get(0).getExonImageBuffer().setFont(Draw.defaultFont);
}
void drawMutations(int ylevel) {
	try {
	
	if(varnode != null && Main.drawCanvas.selectedSplit.viewLength < 1000000) {
		
		VarNode vardraw = varnode;

		while(vardraw != null && vardraw.getPosition() < Main.drawCanvas.splits.get(0).end+1) {		
			
			if(Main.drawCanvas.hideNode(vardraw)) {
				vardraw = vardraw.getNext();
				continue;
			}
		
			found = false;
			if(vardraw.getPrev() != null && vardraw.getPrev().getPosition()-vardraw.getPosition() != -1) {
				baselevel = 0;
			}
			if(vardraw.vars == null) {
				System.out.println("vardraw vars == null (drawMutations)" +vardraw.getPosition());
				continue;
			}
			for(int j = 0; j<vardraw.vars.size(); j++) { //Map.Entry<String, ArrayList<SampleNode>> entry : vardraw.vars) {
				entry = vardraw.vars.get(j);
				if(Main.drawCanvas.hideNodeVar(vardraw, entry)) {
					continue;
				}
				mutcount = 0;
				for(int i = 0; i<entry.getValue().size(); i++) {					
					if(!Main.drawCanvas.hideVar(entry.getValue().get(i))) {
						mutcount++;						
					}
				}
				if(mutcount == 0) {
				
					continue;
				}
				baselevel++;
				base = entry.getKey();	
				 
				 if(vardraw.indel && entry.getKey().length() > 1) {
					 mutScreenPos = (int)((vardraw.getPosition()+2-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel);
					 
				 }
				 else {
					
					 mutScreenPos = (int)((vardraw.getPosition()+1-Main.drawCanvas.splits.get(0).start)*Main.drawCanvas.splits.get(0).pixel);
				 }
				if(vardraw.coding) {
					
					/*foundexon = false;
					
					for(int i = 0; i< vardraw.getExons().size(); i++) {
						
							
						if(vardraw.indel && entry.getKey().length() > 1) {
							if(vardraw.getPosition()+1 > vardraw.getExons().get(i).getTranscript().getCodingStart() && vardraw.getPosition()+1 < vardraw.getExons().get(i).getTranscript().getCodingEnd()) {
								if(vardraw.getPosition()+1 >= vardraw.getExons().get(i).getStart()-2 && vardraw.getPosition()+1 < vardraw.getExons().get(i).getEnd()+2) {
									foundexon = true;
									foundexonindex = i;
									
									break;
								}
							}							
						}
						else {
							if(vardraw.getPosition() >= vardraw.getExons().get(i).getTranscript().getCodingStart() && vardraw.getPosition() < vardraw.getExons().get(i).getTranscript().getCodingEnd()) {
								foundexon = true;
								foundexonindex = i;
								break;
							}
						}
					}
					
					if(foundexon) {
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(this.reverseExon);
						
					}
					else {
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(varColor);
					}		
					*/
					if(base.length() > 1) {
						if(!varColor.equals(forwardExon)) {
							varColor = forwardExon;
						}
					}
					else {
					
						if(!varColor.equals(reverseExon)) {
							varColor = reverseExon;
						}						
					}
					if(!Main.drawCanvas.splits.get(0).getExonImageBuffer().getColor().equals(varColor)) {
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(varColor);
					}
					Main.drawCanvas.splits.get(0).getExonImageBuffer().fillRect(mutScreenPos, 0, (int)Main.drawCanvas.splits.get(0).pixel+1, this.getHeight());
					
				}
				else {		
					if(!varColor.equals(gray)) {
						varColor = gray;
					}
					if(!Main.drawCanvas.splits.get(0).getExonImageBuffer().getColor().equals(varColor)) {
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(varColor);
					}
					Main.drawCanvas.splits.get(0).getExonImageBuffer().fillRect(mutScreenPos, 0, (int)Main.drawCanvas.splits.get(0).pixel+1, this.getHeight());	
					
					if(Main.drawCanvas.splits.get(0).viewLength < 200) {
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setFont(seqFont);
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(Color.black);		
						if(base.length() > 1) {
							Main.drawCanvas.splits.get(0).getExonImageBuffer().drawString(mutcount +" " +base, mutScreenPos, 27*baselevel);
						}
						else {
							Main.drawCanvas.splits.get(0).getExonImageBuffer().drawString(mutcount +" " +Main.getBase.get(vardraw.getRefBase()) +"->" +base, mutScreenPos, 27*baselevel);
							
						}
					//	Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(varColor);		
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setFont(Draw.defaultFont);
					}
				
					
				}
				
				if( Main.drawCanvas.splits.get(0).viewLength < 200 && vardraw.getExons() != null ) {
					 //TODO draw aminochanges	 				
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setFont(seqFont);
						Main.drawCanvas.splits.get(0).getExonImageBuffer().setColor(Color.black);		
						
						for(int i = 0; i< vardraw.getExons().size(); i++) {
							if(vardraw.getExons().get(i).getTranscript().getGene().showIsoforms()) {
								
									Main.drawCanvas.splits.get(0).getExonImageBuffer().drawString(mutcount +" " +getChange(vardraw,base,vardraw.getExons().get(i)), mutScreenPos, vardraw.getExons().get(i).getTranscript().ypos+27*baselevel);
									
								
							}
							else {
								if(vardraw.getExons().get(i).getTranscript().equals(vardraw.getExons().get(i).getTranscript().getGene().getLongest())) {
									
										Main.drawCanvas.splits.get(0).getExonImageBuffer().drawString(mutcount +" " +getChange(vardraw,base,vardraw.getExons().get(i)), mutScreenPos, vardraw.getExons().get(i).getTranscript().ypos+27*baselevel);
									
								}
							}
						}					
						
					}
					Main.drawCanvas.splits.get(0).getExonImageBuffer().setFont(Draw.defaultFont);					
				
				}
				if(vardraw != null) {
					vardraw = vardraw.getNext();
				
				}		
		}		
		vardraw = null;
	}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
	Main.drawCanvas.splits.get(0).getExonImageBuffer().setFont(Draw.defaultFont);
}
/*
String getAminoChange(VarNode node, String base) {
	
	
	if(node.getPosition() >= Main.drawCanvas.selectedSplit.start && node.getPosition() < Main.drawCanvas.selectedSplit.end && node.getExons().size() > 0 && Main.drawCanvas.selectedSplit.viewLength < 200) {
		
		if(base.length() > 1) {
			if(base.substring(3).matches("\\d+")) {
				baselength = Integer.parseInt(base.substring(3));
			}
			else {
				baselength = base.substring(3).length();
			}
									
			if(baselength % 3 == 0) {
				return "if-" +base;
			}
			else {
				return "fs-" +base;
			}
		}
		Transcript.Exon exon = node.getExons().get(0);
		if(exon.getTranscript().getStrand()) {
			if(exon.getFirstAmino() == 1) {
				phase = (node.getPosition()-exon.getTranscript().getCodingStart())%3;
			}
			else {
				phase = (node.getPosition()-(exon.getStart()+exon.getStartPhase()))%3;
			}
		
			if(phase < 0 || (phase == 0 && (exon.getEnd()-node.getPosition()<3))) {
				
				return "";
			}
			else {
				if(node.getCodon() == null) {
					node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getItemAt(Main.selectedChrom),node.getPosition()-phase, node.getPosition()-phase+3, Main.referenceFile)));
				}
				char[] array = node.getCodon().toCharArray();
				
				array[phase] = base.charAt(0);
				return MethodLibrary.getAminoAcid(new String(array));
			}
		}
		else {
			if(exon.getFirstAmino() == 1) {
				phase = (exon.getTranscript().getCodingEnd()-node.getPosition())%3;
			}
			else {
				phase = ((exon.getEnd()-exon.getStartPhase())-node.getPosition())%3;
			
			}
			if(phase > 0) {
				phase = 3-phase;
			}
			
			if(phase < 0 || (phase == 0 && (exon.getEnd()-node.getPosition()<3))) {
				
				return "";
			}
			else {
				
				if(node.getCodon() == null) {
					node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getItemAt(Main.selectedChrom),node.getPosition()-phase, node.getPosition()-phase+3, Main.referenceFile)));
				}
				
				char[] array = node.getCodon().toCharArray();
				array[phase] = base.charAt(0);
				return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
			}
		}
	}	
	
	return "";
}
*/
String getChange(VarNode node, String base, Transcript.Exon exon) {
	try {
	
		aminopos = 0;	
		if(exon == null) {
			return "";
		}
		
		if(base.length() > 1) {
			if(exon.getTranscript().getStrand()) {
				if(exon.getFirstAmino() == 1) {					
					aminopos = (node.getPosition()+1-exon.getTranscript().getCodingStart())/3+1;
				}
				else {					
					aminopos = (node.getPosition()+1-(exon.getStart()+exon.getStartPhase()))/3+exon.getFirstAmino();
					
				}
			}
			else {
				if(exon.getFirstAmino() == 1) {
					
					aminopos = (exon.getTranscript().getCodingEnd()-(node.getPosition()+1))/3+1;
				}
				else {
					
					aminopos = ((exon.getEnd()-exon.getStartPhase()-1)-(node.getPosition()+1))/3+exon.getFirstAmino();
					
				}
			}
			
			if(base.substring(3).matches("\\d+")) {
				baselength = Integer.parseInt(base.substring(3));
			}
			else {
				baselength = base.substring(3).length();
			}
			
			if(node.getPosition()+1 < exon.getTranscript().getCodingStart() || node.getPosition()+1 > exon.getTranscript().getCodingEnd()) {
				
				return base +" (UTR)";
			}
			
			if(node.getPosition()+1 > exon.getEnd()) {			
				if(node.getPosition()+2 - exon.getEnd() > 2) {
					return "";
				}
				return "spl" +(node.getPosition()+2 - exon.getEnd()) +"-" +base;	
			}
			else if (node.getPosition()+1 < exon.getStart()) {	
				
				if(base.contains("ins")) {
					return "spl" +(exon.getStart()-(node.getPosition()+1)) +"-" +base;	
				}
				else {					
					if(((node.getPosition()+1+baselength)-exon.getStart()) % 3 == 0 ) {
						return "spl-" +base +"-if";
					}
					else {
						return "spl-" +base +"-fs";
					}
				}
			}
			else if(base.contains("del") && node.getPosition()+1+baselength >= exon.getEnd()) {
				
				if((baselength-(node.getPosition()+1+baselength - exon.getEnd())) % 3 == 0 ) {
					
					return "if-spl-" +base;
				}
				else {
				
					return "fs-spl-" +base;
				}
				
			}
			else if(baselength % 3 == 0) {
				
				return "if-" +base;				
			}
			else {
				return "fs-" +base;
			}			
		}
		
		
		if(node.getPosition() < exon.getTranscript().getCodingStart() || node.getPosition() >= exon.getTranscript().getCodingEnd()) {
			if(node.getPosition() >= exon.getEnd()) {			
				return "spl" +(node.getPosition()+1 - exon.getEnd()) +" (UTR)";	
			}
			else if (node.getPosition() < exon.getStart()) {			
				return "spl" +(exon.getStart()-(node.getPosition()))+" (UTR)";	
			}
			else {
				return base +" (UTR)";
			}
		}
		if(node.getPosition() >= exon.getEnd()) {			
			return "spl" +(node.getPosition()+1 - exon.getEnd());	
		}
		else if (node.getPosition() < exon.getStart()) {			
			return "spl" +(exon.getStart()-(node.getPosition()));	
		}
		if(exon.getTranscript().getStrand()) {
			
			if(exon.getFirstAmino() == 1) {
				phase = (node.getPosition()-exon.getTranscript().getCodingStart())%3;
				aminopos = (node.getPosition()-exon.getTranscript().getCodingStart())/3+1;
			}
			else {
				phase = (node.getPosition()-(exon.getStart()+exon.getStartPhase()))%3;
				aminopos = (node.getPosition()-(exon.getStart()+exon.getStartPhase()))/3+exon.getFirstAmino();
				
			}
			// SPLICE CODON
			if(phase < 0 || (exon.getEnd() - node.getPosition() <= exon.getEndPhase())) {
				
				if(exon.getEnd() - node.getPosition() <= exon.getEndPhase()) {
					if(phase == 1) {
						if(node.getCodon() == null) {
							nextExonStart = exon.getTranscript().getExons()[exon.getNro()].getStart();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-2, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), nextExonStart, nextExonStart+1, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(new String(array));	
					}
					else {
						if(exon.getEndPhase() == 2) {
							if(node.getCodon() == null) {
								nextExonStart = exon.getTranscript().getExons()[exon.getNro()].getStart();
								node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-2, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), nextExonStart, nextExonStart+1, Main.referenceFile))));
							}
							array = node.getCodon().toCharArray();
							array[phase] = base.charAt(0);
							return MethodLibrary.getAminoAcid(new String(array));	
						}
						else {
							if(node.getCodon() == null) {
								nextExonStart = exon.getTranscript().getExons()[exon.getNro()].getStart();
								node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-1, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), nextExonStart, nextExonStart+2, Main.referenceFile))));
							}
							array = node.getCodon().toCharArray();
							array[phase] = base.charAt(0);
							return MethodLibrary.getAminoAcid(new String(array));	
						}
						
					}
			//		System.out.println(Main.chromosomeDropdown.getSelectedItem().toString() +":" +node.getPosition() +" " +phase);
					
				}				
				else if(phase == -2) {
					phase = 1;
					if(node.getCodon() == null) {
						prevExonEnd = exon.getTranscript().getExons()[exon.getNro()-2].getEnd();
						node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd-1, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getStart(), exon.getStart()+2, Main.referenceFile))));
					}
					array = node.getCodon().toCharArray();
					array[phase] = base.charAt(0);
					return MethodLibrary.getAminoAcid(new String(array));				
					
				}
				else {
					phase = 2;
					if(node.getPosition() == exon.getStart()) {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getNro()-2].getEnd();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd-2, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getStart(), exon.getStart()+1, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(new String(array));
					}
					else {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getNro()-2].getEnd();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd-1, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getStart(), exon.getStart()+2, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(new String(array));
					}					
				}			
			}
			else {
				
				if(node.getCodon() == null) {
					node.setCodon(new String(this.getSeq(exon.getTranscript().getChrom(),node.getPosition()-phase, node.getPosition()-phase+3, Main.referenceFile)));
				}
				
				array = node.getCodon().toCharArray();
				array[phase] = base.charAt(0);				
				return MethodLibrary.getAminoAcid(new String(array));
			}
		}
		else {			
			
			if(exon.getFirstAmino() == 1) {
				phase = (exon.getTranscript().getCodingEnd()-node.getPosition())%3;
				aminopos = (exon.getTranscript().getCodingEnd()-(node.getPosition()+1))/3+1;
				
			}
			else {
				phase = ((exon.getEnd()-exon.getStartPhase())-node.getPosition())%3;
				aminopos = ((exon.getEnd()-exon.getStartPhase()-1)-node.getPosition())/3+exon.getFirstAmino();
				
			}
			if(phase > 0) {
				phase = 3-phase;
			}
			if(phase < 0 || (phase == 0 && (exon.getEnd()-node.getPosition()<3)) || (node.getPosition()-exon.getStart() < exon.getEndPhase())) {
				
				if(node.getPosition()-exon.getStart() < exon.getEndPhase()) {
					
					if(exon.getEndPhase() == 1) {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()-1].getEnd();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(),prevExonEnd-2, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getStart(),  exon.getStart()+1, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
					}
					else {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()-1].getEnd();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(),prevExonEnd-1, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(),  exon.getStart(),  exon.getStart()+2, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
									
					}
					
				}
				else if(phase == -1) {
					
					phase  = 1;
					if(node.getCodon() == null) {
						prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()+1].getStart();
						node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-2, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd, prevExonEnd+1, Main.referenceFile))));
					}
					array = node.getCodon().toCharArray();
					array[phase] = base.charAt(0);
					return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
				}
				else {
					if(exon.getEnd()-node.getPosition() == 2) {
						
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()+1].getStart();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-2, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd, prevExonEnd+1, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
			
					}
					else {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()+1].getStart();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-1, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd, prevExonEnd+2, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
			
					}
				}				
				
			}
			else {
				if(node.getCodon() == null) {
					node.setCodon(new String(this.getSeq(exon.getTranscript().getChrom(),node.getPosition()-phase, node.getPosition()-phase+3, Main.referenceFile)));
				}
				
				array = node.getCodon().toCharArray();
				array[phase] = base.charAt(0);
				return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
			}
		}
	}
	catch(Exception e) {
		//System.out.println(node.getExons().get(0).getTranscript().getGenename() +" " +node.getPosition());
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
		return "";
	}

}
String getAminoChange(VarNode node, String base, Transcript.Exon exon) {
	try {
		aminopos = 0;	
		
		if(exon == null) {
			return "";
		}
		if(base.length() > 1) {
			if(exon.getTranscript().getStrand()) {
				if(exon.getFirstAmino() == 1) {					
					aminopos = (node.getPosition()+1-exon.getTranscript().getCodingStart())/3+1;
				}
				else {					
					aminopos = (node.getPosition()+1-(exon.getStart()+exon.getStartPhase()))/3+exon.getFirstAmino();
					
				}
			}
			else {
				if(exon.getFirstAmino() == 1) {
					
					aminopos = (exon.getTranscript().getCodingEnd()-(node.getPosition()+1))/3+1;
				}
				else {
					
					aminopos = ((exon.getEnd()-exon.getStartPhase()-1)-(node.getPosition()+1))/3+exon.getFirstAmino();
					
				}
			}
			if(base.substring(3).matches("\\d+")) {
				baselength = Integer.parseInt(base.substring(3));
			}
			else {
				baselength = base.substring(3).length();
			}
			
			if(node.getPosition()+1 < exon.getTranscript().getCodingStart() || node.getPosition()+1 > exon.getTranscript().getCodingEnd()) {
				
				return base +" (UTR)";
			}
			if(node.getPosition()+1 > exon.getEnd()) {			
				if(node.getPosition()+2 - exon.getEnd() > 2) {
					return "";
				}
				return "spl" +(node.getPosition()+2 - exon.getEnd()) +"-" +base;	
			}
			else if (node.getPosition()+1 < exon.getStart()) {	
				
				if(base.contains("ins")) {
					return "spl" +(exon.getStart()-(node.getPosition()+1)) +"-" +base;	
				}
				else {					
					if(((node.getPosition()+1+baselength)-exon.getStart()) % 3 == 0 ) {
						return "spl-" +base +"-if";
					}
					else {
						return "spl-" +base +"-fs";
					}
				}
			}
			else if(base.contains("del") && node.getPosition()+1+baselength >= exon.getEnd()) {
				
				if((baselength-(node.getPosition()+1+baselength - exon.getEnd())) % 3 == 0 ) {
					
					return aminopos +"-if-spl-" +base;
				}
				else {
				
					return aminopos +"-fs-spl-" +base;
				}
				
			}
			else if(baselength % 3 == 0) {
				
				return aminopos +"-if-" +base;				
			}
			else {
				return aminopos +"-fs-" +base;
			}			
		}
		if(node.getPosition() < exon.getTranscript().getCodingStart() || node.getPosition() >= exon.getTranscript().getCodingEnd()) {
			
			if(node.getPosition() >= exon.getEnd()) {			
				return "spl" +(node.getPosition()+1 - exon.getEnd()) +" (UTR)";	
			}
			else if (node.getPosition() < exon.getStart()) {			
				return "spl" +(exon.getStart()-(node.getPosition()))+" (UTR)";	
			}
			else {
				return base +" (UTR)";
			}
		}
		if(node.getPosition() >= exon.getEnd()) {			
			return "spl" +(node.getPosition()+1 - exon.getEnd());	
		}
		else if (node.getPosition() < exon.getStart()) {			
			return "spl" +(exon.getStart()-(node.getPosition()));	
		}
		
		if(exon.getTranscript().getStrand()) {
			
			if(exon.getFirstAmino() == 1) {
				phase = (node.getPosition()-exon.getTranscript().getCodingStart())%3;
				aminopos = (node.getPosition()-exon.getTranscript().getCodingStart())/3+1;
			}
			else {
				phase = (node.getPosition()-(exon.getStart()+exon.getStartPhase()))%3;
				aminopos = (node.getPosition()-(exon.getStart()+exon.getStartPhase()))/3+exon.getFirstAmino();
				
			}
			// SPLICE CODON
			if(phase < 0 || (exon.getEnd() - node.getPosition() <= exon.getEndPhase())) {
				if(exon.getNro() >= exon.getTranscript().getExons().length ) {
					return "";
				}
				if(exon.getEnd() - node.getPosition() <= exon.getEndPhase()) {
					if(phase == 1) {
						if(node.getCodon() == null) {
							nextExonStart = exon.getTranscript().getExons()[exon.getNro()].getStart();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-2, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), nextExonStart, nextExonStart+1, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(node.getCodon()) +(aminopos) +MethodLibrary.getAminoAcid(new String(array));	
					}
					else {
						if(exon.getEndPhase() == 2) {
							if(node.getCodon() == null) {
								
								nextExonStart = exon.getTranscript().getExons()[exon.getNro()].getStart();
								node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-2, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), nextExonStart, nextExonStart+1, Main.referenceFile))));
								
							}
							array = node.getCodon().toCharArray();
							array[phase] = base.charAt(0);
							return MethodLibrary.getAminoAcid(node.getCodon()) +(aminopos) +MethodLibrary.getAminoAcid(new String(array));	
						}
						else {
							if(node.getCodon() == null) {
								nextExonStart = exon.getTranscript().getExons()[exon.getNro()].getStart();
								node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-1, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), nextExonStart, nextExonStart+2, Main.referenceFile))));
							}
							array = node.getCodon().toCharArray();
							array[phase] = base.charAt(0);
							return MethodLibrary.getAminoAcid(node.getCodon()) +(aminopos) +MethodLibrary.getAminoAcid(new String(array));	
						}
						
					}
					
				}				
				else if(phase == -2) {
					phase = 1;
					if(node.getCodon() == null) {
						prevExonEnd = exon.getTranscript().getExons()[exon.getNro()-2].getEnd();
						node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd-1, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getStart(), exon.getStart()+2, Main.referenceFile))));
					}
					array = node.getCodon().toCharArray();
					array[phase] = base.charAt(0);
					return MethodLibrary.getAminoAcid(node.getCodon()) +(aminopos-1) +MethodLibrary.getAminoAcid(new String(array));				
					
				}
				else {
					phase = 2;
					if(node.getPosition() == exon.getStart()) {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getNro()-2].getEnd();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd-2, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getStart(), exon.getStart()+1, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(node.getCodon()) +(aminopos-1) +MethodLibrary.getAminoAcid(new String(array));
					}
					else {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getNro()-2].getEnd();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd-1, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getStart(), exon.getStart()+2, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(node.getCodon()) +(aminopos-1) +MethodLibrary.getAminoAcid(new String(array));
					}					
				}			
			}
			else {
				
				if(node.getCodon() == null) {
					node.setCodon(new String(this.getSeq(exon.getTranscript().getChrom(),node.getPosition()-phase, node.getPosition()-phase+3, Main.referenceFile)));
				}
				
				array = node.getCodon().toCharArray();
				array[phase] = base.charAt(0);				
				return MethodLibrary.getAminoAcid(node.getCodon()) +aminopos +MethodLibrary.getAminoAcid(new String(array));
			}
		}
		else {			
			
			if(exon.getFirstAmino() == 1) {
				phase = (exon.getTranscript().getCodingEnd()-node.getPosition())%3;
				aminopos = (exon.getTranscript().getCodingEnd()-(node.getPosition()+1))/3+1;
				
			}
			else {
				phase = ((exon.getEnd()-exon.getStartPhase())-node.getPosition())%3;
				aminopos = ((exon.getEnd()-exon.getStartPhase()-1)-node.getPosition())/3+exon.getFirstAmino();
				
			}
			if(phase > 0) {
				phase = 3-phase;
			}
			if(phase < 0 || (phase == 0 && (exon.getEnd()-node.getPosition()<3)) || (node.getPosition()-exon.getStart() < exon.getEndPhase())) {
				
				if(node.getPosition()-exon.getStart() < exon.getEndPhase()) {
					
					if(exon.getEndPhase() == 1) {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()-1].getEnd();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(),prevExonEnd-2, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getStart(),  exon.getStart()+1, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(node.getCodon())) +(aminopos) +MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
					}
					else {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()-1].getEnd();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(),prevExonEnd-1, prevExonEnd, Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(),  exon.getStart(),  exon.getStart()+2, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(node.getCodon())) +(aminopos) +MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
									
					}
					
				}
				else if(phase == -1) {
					
					phase  = 1;
					if(node.getCodon() == null) {
						prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()+1].getStart();
						node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-2, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd, prevExonEnd+1, Main.referenceFile))));
					}
					array = node.getCodon().toCharArray();
					array[phase] = base.charAt(0);
					return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(node.getCodon())) +(aminopos+1) +MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
				}
				else {
					if(exon.getEnd()-node.getPosition() == 2) {
						
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()+1].getStart();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-2, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd, prevExonEnd+1, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(node.getCodon())) +(aminopos+1) +MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
			
					}
					else {
						if(node.getCodon() == null) {
							prevExonEnd = exon.getTranscript().getExons()[exon.getTranscript().getExons().length -exon.getNro()+1].getStart();
							node.setCodon(new String(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), exon.getEnd()-1, exon.getEnd(), Main.referenceFile).append(this.getSeq(Main.chromosomeDropdown.getSelectedItem().toString(), prevExonEnd, prevExonEnd+2, Main.referenceFile))));
						}
						array = node.getCodon().toCharArray();
						array[phase] = base.charAt(0);
						return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(node.getCodon())) +(aminopos+1) +MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
			
					}
				}				
				
			}
			else {
				if(node.getCodon() == null) {
					node.setCodon(new String(this.getSeq(exon.getTranscript().getChrom(),node.getPosition()-phase, node.getPosition()-phase+3, Main.referenceFile)));
				}
				
				array = node.getCodon().toCharArray();
				array[phase] = base.charAt(0);
				return MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(node.getCodon())) +aminopos +MethodLibrary.getAminoAcid(MethodLibrary.reverseComplement(new String(array)));
			}
		}
	}
	catch(Exception e) {
		//System.out.println(node.getExons().get(0).getTranscript().getGenename() +" " +node.getPosition());
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
		return "";
	}

}

int maxWidth(String[] list) {
	int max = 0, index = 0;
	for(int i = 0; i<list.length; i++) {
		if(list[i].length() > max) {
			max = list[i].length();
			index = i;
		}
	}
	return index;	
	
}
int getMousePos(int mousex) {
	return (int)(Main.drawCanvas.selectedSplit.start+((mousex-Main.drawCanvas.selectedSplit.offset)/Main.drawCanvas.selectedSplit.pixel));
}
	@Override
public void mouseDragged(MouseEvent event) {
	
	switch(event.getModifiers()) {	
	
		case InputEvent.BUTTON1_MASK: {	
			this.mouseX = event.getX();
			this.mouseY = event.getY();
			if(resizeSidebar) {
				Main.drawCanvas.resizeSidebar(mouseX);
		//		updateExons = true;
		//		repaint();
				break;
			}
			else if(!seqDrag && (mouseY < Main.chromScroll.getViewport().getHeight()-15 || Main.drawCanvas.selectedSplit.viewLength > 300)) {
				if(Main.drawCanvas.lineZoomer) {	
					
					 if (Main.drawCanvas.selectedSplit.start > 1 || Main.drawCanvas.selectedSplit.end < Main.drawCanvas.selectedSplit.chromEnd ) {						
						 Main.drawCanvas.gotoPos(Main.drawCanvas.selectedSplit.start-(Main.drawCanvas.tempDrag-mouseX)/Main.drawCanvas.selectedSplit.pixel*2, Main.drawCanvas.selectedSplit.end+(Main.drawCanvas.tempDrag-mouseX)/Main.drawCanvas.selectedSplit.pixel*2);
					 }
				
					 Main.drawCanvas.tempDrag = mouseX;
					updateExons = true;
					repaint();
					Draw.updatevars = true;
					Main.drawCanvas.repaint();
				}
				else {							
					zoomDrag = true;				
					
					repaint();
					return;
				}
			}
			else {
				seqDrag = true;
				
					repaint();
				return;
			}
			
			break;
		}
		case InputEvent.BUTTON3_MASK: {	
				if((int)Main.drawCanvas.selectedSplit.start == 1 && (int)Main.drawCanvas.selectedSplit.end == Main.drawCanvas.selectedSplit.chromEnd) {
					break;
				}			
				Main.drawCanvas.mouseDrag = true;
				Main.drawCanvas.moveX = event.getX();
				Main.drawCanvas.drag(Main.drawCanvas.moveX);
		
			
			break;
		}
		case 17: {
			if((int)Main.drawCanvas.selectedSplit.start == 1 && (int)Main.drawCanvas.selectedSplit.end == Main.drawCanvas.selectedSplit.chromEnd) {
				break;
			}			
			Main.drawCanvas.mouseDrag = true;
			Main.drawCanvas.moveX = event.getX();
			Main.drawCanvas.drag(Main.drawCanvas.moveX);		
			break;
		}
	}
	
}

	@Override
public void mouseMoved(MouseEvent event) {
		this.mouseX = event.getX();
		this.mouseY = event.getY();
		if((this.mouseX-Main.sidebarWidth)/(Main.drawCanvas.getDrawWidth()) > -1 && (this.mouseX-Main.sidebarWidth)/(Main.drawCanvas.getDrawWidth()) < Main.drawCanvas.splits.size()) {
			if(Main.drawCanvas.selectedSplit != Main.drawCanvas.splits.get((this.mouseX-Main.sidebarWidth)/(Main.drawCanvas.getDrawWidth()))) {
				Main.drawCanvas.selectedSplit = Main.drawCanvas.splits.get((this.mouseX-Main.sidebarWidth)/(Main.drawCanvas.getDrawWidth()));
				if(Main.drawCanvas.selectedSplit == null) {
					Main.drawCanvas.selectedSplit = Main.drawCanvas.splits.get(0);
				}
				/*
				if(Main.drawCanvas.selectedSplit > Main.drawCanvas.splits.size()-1) {
					Main.drawCanvas.selectedSplit = Main.drawCanvas.splits.size()-1;
				}
				else if(Main.drawCanvas.selectedSplit < 0) {
					Main.drawCanvas.selectedSplit = 0;
				}*/
			}
		}
	if(event.getY() > Main.chromScroll.getViewport().getHeight() -15 && Main.drawCanvas.selectedSplit.viewLength < 300) {
		
		if(getCursor().getType() != Cursor.TEXT_CURSOR) {
			setCursor(Cursor.getPredefinedCursor(Cursor.TEXT_CURSOR));
		}
	}
	else {
		if(getCursor().getType() == Cursor.TEXT_CURSOR && getCursor().getType() != Cursor.DEFAULT_CURSOR && mouseY > cytoHeight) {
			
			setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
		}
		
	}
	mouseRect.setBounds(event.getX()-Main.drawCanvas.selectedSplit.offset, event.getY()-this.cytoHeight, 2,2);
	
	repaint();

}

@Override
public void mouseClicked(MouseEvent event) {
	switch(event.getModifiers()) {	
	
		case InputEvent.BUTTON1_MASK: {
			if(clickedExon != null && foundcursor > -1) {
				if(foundcursor == 0) {
					Main.drawCanvas.clearReads();
					Main.drawCanvas.gotoPos(clickedExon.getStart()+1, clickedExon.getEnd()+1);
				}
				else if(foundcursor == 1) {
					Main.drawCanvas.clearReads();
					Main.drawCanvas.gotoPos(clickedExon.transcript.getStart()+1, clickedExon.transcript.getEnd()+1);
				}
				else if(foundcursor == 2) {
					if(exonString[2].contains("ENS")) {
						Main.gotoURL("http://ensembl.org/Multi/Search/Results?q=" +exonString[2].substring(exonString[2].indexOf("ENS")).split(" ")[0]);
					}
					else {
						Main.gotoURL("http://www.ncbi.nlm.nih.gov/gene/?term=" +exonString[2].split(":")[1].split(",")[0].split(" ")[0]);
					}
				}
				else if(foundcursor == 3) {
					if(exonString[3].contains("ENS")) {
						Main.gotoURL("http://ensembl.org/Multi/Search/Results?q=" +exonString[3].substring(exonString[3].indexOf("ENS")).split(" ")[0]);
					}
				}
				else if(foundcursor == 5) {
					Main.gotoURL("http://www.genecards.org/cgi-bin/carddisp.pl?gene=" +exonString[1]);
				}
				updateExons = true;
						repaint();
				Draw.updateReads = true;
				Draw.updatevars = true;
				Main.drawCanvas.repaint();
				break;
			}
			else {
		//	System.out.println(MethodLibrary.formatNumber((int)(Main.drawCanvas.selectedSplit.start +((mouseX-Main.drawCanvas.selectedSplit.offset)/(double)Main.drawCanvas.selectedSplit.pixel))));
	//			System.out.println(MethodLibrary.getRegion((int)(Main.drawCanvas.selectedSplit.start +((mouseX-Main.drawCanvas.selectedSplit.offset)/(double)Main.drawCanvas.selectedSplit.pixel)), Main.drawCanvas.selectedSplit,0));
			}
			
			break;
		}	
		case InputEvent.BUTTON3_MASK: {	
			
			if(selectedExon != null && selectedExon.transcript.getGene().getTranscripts().size() > 1 && !selectedExon.transcript.getGene().showIsoforms()) {
				
				selectedExon.transcript.getGene().setShowIsoforms(true);
				
		
				updateExons = true;
				repaint();
			}
			else if(selectedExon != null && selectedExon.transcript.getGene().getTranscripts().size() > 1 && selectedExon.transcript.getGene().showIsoforms()) {
				selectedExon.transcript.getGene().setShowIsoforms(false);				
		
				updateExons = true;
				repaint();
			}
			else {
				String copy = Main.drawCanvas.selectedSplit.chrom +":" +getPosition(mouseX-Main.drawCanvas.selectedSplit.offset, Main.drawCanvas.selectedSplit);
				StringSelection stringSelection= new StringSelection(copy);
				Main.putMessage("Position " +copy +" copied to clipboard.");
				Clipboard clpbrd = Toolkit.getDefaultToolkit().getSystemClipboard();
				clpbrd.setContents(stringSelection, null);
			}
			
			break;
		}	
	}
}

@Override
public void mouseEntered(MouseEvent arg0) {
	
	
}

@Override
public void mouseExited(MouseEvent arg0) {
	
	
}

@Override
public void mousePressed(MouseEvent event) {
	Logo.frame.setVisible(false);
	
	Main.drawCanvas.pressX = event.getX();
	this.pressY = event.getY();
	switch(event.getModifiers()) {	
	
		case InputEvent.BUTTON1_MASK: {	
	
			if(selectedExon != null) {		
				clickedSplit = Main.drawCanvas.selectedSplit;
				clickedExon = selectedExon;
				repaint();
			}
			if(!lineZoomer && mouseX > Main.sidebarWidth-3 && mouseX < Main.sidebarWidth+3) {
				resizeSidebar = true;
			}
			else {
				resizeSidebar = false;
			}
			
			break;
		}
		case InputEvent.BUTTON3_MASK: {
			
			clickedExon = null;			
			this.zoomDrag = false;
			updateExons = true;
			repaint();
		}
	}
	Main.drawCanvas.tempDrag = Main.drawCanvas.pressX;	
}

	@Override
public void mouseReleased(MouseEvent arg0) {
		
		if(seqDrag) {
		//	System.out.println(Main.drawCanvas.selectedSplit.sequence.substring(getMousePos(pressX)-1-Main.drawCanvas.selectedSplit.seqStartPos, getMousePos(mouseX)-1-Main.drawCanvas.selectedSplit.seqStartPos+1));
			seqDrag = false;
			if(mouseX-Main.drawCanvas.pressX > 0)  {
				String myString = new String( Arrays.copyOfRange(Main.drawCanvas.selectedSplit.getReference().getSeq(),(getMousePos(Main.drawCanvas.pressX)-1-Main.drawCanvas.selectedSplit.getReference().getStartPos()), getMousePos(mouseX)-1-Main.drawCanvas.selectedSplit.getReference().getStartPos()+1));
				StringSelection stringSelection = new StringSelection(myString);
				Clipboard clpbrd = Toolkit.getDefaultToolkit().getSystemClipboard();
				clpbrd.setContents(stringSelection, null);
				message = "Sequence copied to clipboard.";
				timer = System.currentTimeMillis();
			}
				repaint();
			
		}
		if(resizeSidebar) {
			Main.drawCanvas.resizeCanvas(Main.drawCanvas.getWidth(), Main.drawCanvas.getHeight());
			Main.controlDraw.repaint();
			Main.bedCanvas.repaint();
			resizeSidebar = false;
		}
		Main.drawCanvas.lineZoomer = false;
		lineZoomer = false;
		if(Main.drawCanvas.mouseDrag) {
			
			Main.drawCanvas.mouseDrag = false;
		}
		
	if(zoomDrag) {
		if(mouseY <= cytoHeight) {
			if(mouseX-Main.drawCanvas.pressX > 0) {
				
				Main.drawCanvas.gotoPos((Main.drawCanvas.pressX-Main.drawCanvas.selectedSplit.offset)/(Main.drawCanvas.getDrawWidth()/(double)chromPos.get(Main.refchrom +Main.drawCanvas.selectedSplit.chrom)), (mouseX-Main.drawCanvas.selectedSplit.offset)/(Main.drawCanvas.getDrawWidth()/(double)chromPos.get(Main.refchrom +Main.drawCanvas.selectedSplit.chrom)));
			}
			
		}
		else if(mouseX-Main.drawCanvas.pressX > 0) {
			/*if(split && pressX > (Main.drawCanvas.getDrawWidth())/2) {
				splitClass.gotoPos(splitClass.start+(pressX-(Main.sidebarWidth+(Main.drawCanvas.getDrawWidth())/2))/splitClass.pixel,splitClass.start+(mouseX-(Main.sidebarWidth+(Main.drawCanvas.getDrawWidth())/2))/splitClass.pixel);
				//TODO
			}
			else {
				*/
				Main.drawCanvas.gotoPos(Main.drawCanvas.selectedSplit.start+(Main.drawCanvas.pressX-Main.drawCanvas.selectedSplit.offset)/Main.drawCanvas.selectedSplit.pixel, Main.drawCanvas.selectedSplit.start+(mouseX-Main.drawCanvas.selectedSplit.offset)/Main.drawCanvas.selectedSplit.pixel);
		//	}
		}
	//	System.out.println(pressX/Main.drawCanvas.selectedSplit.pixel +" " + mouseX/Main.drawCanvas.selectedSplit.pixel);
		zoomDrag = false;	
		updateExons = true;
		
		Draw.updatevars = true;
		
	}
	Main.bedCanvas.repaint();
		Main.drawCanvas.repaint();
		repaint();
	
}
String[] getGposColor(String gpos) {
	String[] color = new String[3];
	if(gpos.equals("gneg")) {
		color[0] = "255";
		color[1] = "255";
		color[2] = "255";
		
		return color;		
	}
	else if(gpos.equals("gpos50")) {
		color[0] = "200";
		color[1] = "200";
		color[2] = "200";
		return color;
	}
	else if(gpos.equals("gpos100")) {
		color[0] = "0";
		color[1] = "0";
		color[2] = "0";
		return color;
	}
	else if(gpos.equals("gpos75")) {
		color[0] = "130";
		color[1] = "130";
		color[2] = "130";
		return color;
	}
	else if(gpos.equals("gpos66")) {
		color[0] = "160";
		color[1] = "160";
		color[2] = "160";
		return color;
	}
	else if(gpos.equals("gpos33")) {
		color[0] = "210";
		color[1] = "210";
		color[2] = "210";
		return color;
	}
	else if(gpos.equals("gpos25")) {
		color[0] = "200";
		color[1] = "200";
		color[2] = "200";
		return color;
	}
	else if(gpos.equals("gvar")) {
		color[0] = "220";
		color[1] = "220";
		color[2] = "220";
		return color;
	}
	else if(gpos.equals("acen")) {
		color[0] = "217";
		color[1] = "47";
		color[2] = "39";
		return color;
	}
	else if(gpos.equals("stalk")) {
		color[0] = "100";
		color[1] = "127";
		color[2] = "164";
		return color;
	}
	else if(gpos.equals("gpos")) {
		color[0] = "0";
		color[1] = "0";
		color[2] = "0";
		return color;
	}
	
	return color;
}
public BufferedImage createBands(String chrom) {
	try {
		tempImage = MethodLibrary.toCompatibleImage(new BufferedImage((Main.drawCanvas.getWidth()), cytoHeight, BufferedImage.TYPE_INT_RGB));
		cytoImageBuffer = (Graphics2D)tempImage.getGraphics();
		chromBands.clear();
		
		for(int i = 0; i<bandVector.size(); i++) {
			
				if(!bandVector.get(i)[0].equals(chrom) && !bandVector.get(i)[0].equals("chr" +chrom)) {
					continue;
				}
				
					chromBands.add(bandVector.get(i));
			//	if(!split) {
					bandwidth = (int)((Integer.parseInt(bandVector.get(i)[2]) - Integer.parseInt(bandVector.get(i)[1]))*((Main.drawCanvas.getWidth())/(double)chromPos.get(Main.refchrom +chrom)));
					Xpos = (int)(Integer.parseInt(bandVector.get(i)[1])*((Main.drawCanvas.getWidth())/(double)chromPos.get(chrom)));
			/*	}
				else {
					bandwidth = (int)((Integer.parseInt(bandVector.get(i)[2]) - Integer.parseInt(bandVector.get(i)[1]))*((Main.drawCanvas.getWidth())/2/(double)chromPos.get(Main.refchrom +chrom)));
					Xpos = (int)(Integer.parseInt(bandVector.get(i)[1])*((Main.drawCanvas.getWidth())/2/(double)chromPos.get(chrom)));
				}*/
					if( bandVector.get(i)[4].contains(",")) {
						color = bandVector.get(i)[4].split(",");
					}
					else {
						color = getGposColor(bandVector.get(i)[4]);
					}
				cytoImageBuffer.setColor(new Color(Integer.parseInt(color[0]), Integer.parseInt(color[1]), Integer.parseInt(color[2])));
				
				cytoImageBuffer.fillRect(Xpos, 0, bandwidth, 15);
				
				if(color[0].equals("0")) {
					
					cytoImageBuffer.setColor(Color.white);
				}
				else {
					
					cytoImageBuffer.setColor(Color.black);
				}
			/*	if(bandwidth > 20) {
					
					cytoImageBuffer.drawString(bandVector.get(i)[3], Xpos, cytoHeight-2);
				}		
				*/		
		}
	}
	catch(Exception e) {
		e.printStackTrace();
	}
		return tempImage;
}


/*
public StringBuffer getSequence() {
	return Main.drawCanvas.selectedSplit.sequence;
}
public int getSeqStart() {
	return Main.drawCanvas.selectedSplit.seqStartPos;
}*/
/*
void resizeDraw(int Main.drawCanvas.getDrawWidth(), int this.getHeight()) {
	this.Main.drawCanvas.getDrawWidth() = Main.drawCanvas.getDrawWidth();
	this.this.getHeight() = this.getHeight();
	chromImage = new BufferedImage(Main.drawCanvas.getDrawWidth(), this.getHeight(), BufferedImage.TYPE_INT_RGB);	
	chromImageBuffer = (Graphics2D)chromImage.getGraphics();
	exonImage = new BufferedImage(Main.drawCanvas.getDrawWidth(), this.getHeight(), BufferedImage.TYPE_INT_ARGB);	
	exonImageBuffer = (Graphics2D)exonImage.getGraphics();
	}*/
}
