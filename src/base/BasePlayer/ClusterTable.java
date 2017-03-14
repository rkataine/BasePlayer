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
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Map;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import org.apache.commons.io.FilenameUtils;


public class ClusterTable  extends JPanel implements MouseMotionListener, MouseListener, MouseWheelListener{

	
	private static final long serialVersionUID = 1L;
	BufferedImage bufImage;
	ListSorter sorter = new ListSorter();
	int rowHeight = 15, geneheaderlength = 0;
	int variants = 0;
	Graphics2D buf;
	String chrom = "";
	StringBuffer bedbuffer = new StringBuffer("");
	int width, height;
	ArrayList<AminoEntry> aminoarray = new ArrayList<AminoEntry>();
//	ArrayList<Object[]> sampleArray = new ArrayList<Object[]>();
	ArrayList<Sample> samples = new ArrayList<Sample>();
	int mouseY=0, mouseX=0, pressX=0, pressY=0;
//	VarNode hoverNode;
	final JScrollPane tablescroll;
	ClusterNode hoverNode, selectedNode;
	String[] posSplit, hoverString;
	int samplecount = 0;
	Enumeration<String> e;
	String base;
	
	ArrayList<Object[]> geneheader =  new ArrayList<Object[]>();
	ArrayList<Object[]> header = new ArrayList<Object[]>();
	private int genemutcount;

	VarNode hoverVar;
	private VarNode selectedVar;
	private int listAdd;
	private String[] selectedString;
	private int pointer;
	int headerHover;
	boolean dot = false;
	
	Polygon sortTriangle = new Polygon();
	private int hoverSample = -1;
	private int mutcount;
	ArrayList<SampleNode> vararray = new ArrayList<SampleNode>();
	SampleNode[] controlarray;
	StringBuffer[] bedarray;
	private Color textcolor;
	private double casefreq;
	VarNode varAdd;
	Map.Entry<String, ArrayList<SampleNode>>  entry;
	
	MethodLibrary.controlsorter ctrlsort = new MethodLibrary.controlsorter();

	private Color linecolor;
	private ClusterNode cluster;
	
	private int geneHeaderHover;
	private boolean mouseDrag;
	private int resizeColumn;
	private int dragX;
	private int cases;
	private int firstrow;
	
	ClusterTable(int width, int height, JScrollPane tablescroll) {			
		this.width = width;
		this.height = height;
		this.tablescroll = tablescroll;
		Object[] obj = new Object[3]; obj[0] = "Sample"; 	obj[1] = 10; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
		 obj = new Object[3]; obj[0] = "Mut. count";obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
		 obj = new Object[3]; obj[0] = "Position"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
		 obj = new Object[3]; obj[0] = "Base change"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
		 obj = new Object[3]; obj[0] = "Genotype"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
		 obj = new Object[3]; obj[0] = "Quality"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
		 obj = new Object[3]; obj[0] = "rs-code"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
		 geneheaderlength = geneheader.size();
		 obj = new Object[3]; obj[0] = "ClusterId"; 	obj[1] = 10; obj[2] = (int)((width-10)/7.0); header.add(obj);
		 obj = new Object[3]; obj[0] = "Mut. count";obj[1] = (int)header.get(header.size()-1)[1] + (int)header.get(header.size()-1)[2]; obj[2] = (int)((width-10)/7.0); header.add(obj);
		 obj = new Object[3]; obj[0] = "Width"; 	obj[1] = (int)header.get(header.size()-1)[1] + (int)header.get(header.size()-1)[2]; obj[2] = (int)((width-10)/7.0); header.add(obj);
		 obj = new Object[3]; obj[0] = "Start position"; 	obj[1] = (int)header.get(header.size()-1)[1] + (int)header.get(header.size()-1)[2]; obj[2] = (int)((width-10)/7.0); header.add(obj);
		 obj = new Object[3]; obj[0] = "Variant freq."; 	obj[1] = (int)header.get(header.size()-1)[1] + (int)header.get(header.size()-1)[2]; obj[2] = (int)((width-10)/7.0); header.add(obj);
		 obj = new Object[3]; obj[0] = "Flanking genes"; 	obj[1] = (int)header.get(header.size()-1)[1] + (int)header.get(header.size()-1)[2]; obj[2] = (int)((width-10)/7.0); header.add(obj);
		
		 
		bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width, (int)height, BufferedImage.TYPE_INT_ARGB));	
		buf = (Graphics2D)bufImage.getGraphics();
		this.addMouseListener(this);
		this.addMouseMotionListener(this);
			
		
	}	
	
void resizeTable() {

	if(bufImage.getWidth() < (int)header.get(header.size()-1)[1]+ (int)header.get(header.size()-1)[2]) {
		bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
		buf = (Graphics2D)bufImage.getGraphics();
	}
	

}

void resizeTable(int width) {
	header.get(0)[1] = 0;
	header.get(0)[2] = (int)(width/(double)header.size());
	for(int i = 1 ; i<header.size(); i++) {
		header.get(i)[1] = (int)header.get(i-1)[1]+(int)header.get(i-1)[2];	
		header.get(i)[2] = (int)(width/(double)header.size());		
	}
	
	geneheader.get(0)[2] = (int)((width-10)/geneheader.size());
	
	for(int i = 1; i<geneheader.size(); i++) {
		geneheader.get(i)[1] = (int)geneheader.get(i-1)[1] + (int)geneheader.get(i-1)[2];
		geneheader.get(i)[2] = (int)((width-10)/geneheader.size());
	}
	this.setPreferredSize(new Dimension(width, this.getHeight()));
	this.revalidate();
}
void resizeTable(int column, int amount) {
	
	if(headerHover != -1) {
		if((int)header.get(column-1)[2]+amount > 20) {
			header.get(column-1)[2]=(int)header.get(column-1)[2]+amount;
			for(int i = column;i<header.size(); i++) {
				header.get(i)[1]=(int)header.get(i)[1]+amount;
			}
		}
		if((int)header.get(header.size()-1)[1]+(int)header.get(header.size()-1)[2] > this.getWidth()) {
			if(bufImage.getWidth() < (int)header.get(header.size()-1)[1]+(int)header.get(header.size()-1)[2]) {
				bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
				buf = (Graphics2D)bufImage.getGraphics();
			}
			this.setPreferredSize(new Dimension((int)header.get(header.size()-1)[1]+(int)header.get(header.size()-1)[2], this.getHeight()));
			this.revalidate();
		}
	}
	else {
		if((int)geneheader.get(column-1)[2]+amount > 20) {
			geneheader.get(column-1)[2]=(int)geneheader.get(column-1)[2]+amount;
			for(int i = column;i<geneheader.size(); i++) {
				geneheader.get(i)[1]= (int)geneheader.get(i)[1]+amount;
			}
		}
		if((int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2] > this.getWidth()) {
			if(bufImage.getWidth() < (int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2]) {
				bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
				buf = (Graphics2D)bufImage.getGraphics();
			}
			this.setPreferredSize(new Dimension((int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2], this.getHeight()));
			this.revalidate();
		}
	}
	createPolygon();	
	repaint();
}
void addColumnGeneheader(Object column) {
	
	Object[] obj = new Object[3]; 
	obj[0] = column; 	
	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2];
	obj[2] = 100;
	geneheader.add(obj);
	
	if((int)obj[1]+(int)obj[2] > this.getWidth()) {
		if(bufImage.getWidth() < (int)obj[1]+(int)obj[2]) {
			bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
			buf = (Graphics2D)bufImage.getGraphics();
		}
		this.setPreferredSize(new Dimension((int)obj[1]+(int)obj[2], this.getHeight()));
		this.revalidate();
	}
}

void addHeaderColumn(Object column) {
	Object[] obj = new Object[3]; 
	obj[0] = column; 	
	obj[1] = (int)header.get(header.size()-1)[1] + (int)header.get(header.size()-1)[2];
	obj[2] = 100;
	header.add(obj);
	
	if((int)obj[1]+(int)obj[2] > this.getWidth()) {
		if(bufImage.getWidth() < (int)obj[1]+(int)obj[2]) {
			bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
			buf = (Graphics2D)bufImage.getGraphics();
		}
		
	}
	this.setPreferredSize(new Dimension((int)obj[1]+(int)obj[2], this.getHeight()));
	this.revalidate();
}

void drawScreen(Graphics g) {	
	try {
		
	buf.setColor(Color.black);	
	buf.fillRect(0, 0, this.getWidth(), tablescroll.getViewport().getHeight());		
	genemutcount = 0;	
	hoverVar = null;
	hoverSample = -1;
	headerHover = -1;
	geneHeaderHover = -1;
	
	if(!mouseDrag) {
		resizeColumn = -1;
	}
	if(Main.drawCanvas.clusterNodes != null) {		
	firstrow = tablescroll.getVerticalScrollBar().getValue()/rowHeight -1-Main.drawCanvas.clusterNodes.size();
	
	if(firstrow < 0) {
		firstrow = 0;
	}
	for(int i = 0; i<Main.drawCanvas.clusterNodes.size(); i++) {
		dot = false;
					
		if((i+1+samplecount+Main.drawCanvas.clusterNodes.size())*rowHeight < tablescroll.getVerticalScrollBar().getValue()) {			
			continue;
		}
		
		if(i*rowHeight > tablescroll.getVerticalScrollBar().getValue() + tablescroll.getViewport().getHeight()) {
			break;
		}
		
		if(mouseY >= (rowHeight*(i+genemutcount+1)) && mouseY < (rowHeight*(i+genemutcount+2))) {
			hoverNode = Main.drawCanvas.clusterNodes.get(i);	
		}	
		try {
			buf.setColor(Color.darkGray);
			buf.drawLine(4, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3);	
			
		
			buf.setColor(linecolor);
			cluster = Main.drawCanvas.clusterNodes.get(i);
			if(cluster.varnodes.get(0).getTranscripts() != null) {
				if(!chrom.equals(cluster.varnodes.get(0).getTranscripts().get(0).getChrom())) {
					chrom = cluster.varnodes.get(0).getTranscripts().get(0).getChrom();
				}
			}
			else {
				if(!chrom.equals(cluster.varnodes.get(0).getExons().get(0).transcript.getChrom())) {
					chrom = cluster.varnodes.get(0).getExons().get(0).transcript.getChrom();
				}
			}
			
			for(int c=0; c<header.size(); c++) {
				
				if(Main.drawCanvas.clusterNodes.get(i).equals(hoverNode) || Main.drawCanvas.clusterNodes.get(i).equals(selectedNode)) {
					
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				
				if(c == 0) {
					
					buf.drawString(""+cluster.ID, (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}
				else if(c == 1) { 
					buf.drawString("" +cluster.nodecount, (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}
				else if(c == 2) { 
					buf.drawString("" +cluster.width, (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}
				else if(c == 3) { 
					buf.drawString(chrom +":" +MethodLibrary.formatNumber(cluster.varnodes.get(0).getPosition()), (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}
				else if(c == 4) { 
					buf.drawString("" +MethodLibrary.round((cluster.nodecount/(double)cluster.width),4), (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}	
				else if(c == 5) {
					if(cluster.varnodes.get(0).getExons() != null) {
						
						if(cluster.varnodes.get(0).coding) {
							buf.setColor(Color.red);
							buf.drawString(cluster.varnodes.get(0).getExons().get(0).getTranscript().getGenename() +" (Coding)", (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
							
						}
						else {
							buf.setColor(Color.lightGray);
							buf.drawString(cluster.varnodes.get(0).getExons().get(0).getTranscript().getGenename() +" (UTR)", (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
							
						}
							}
					else if(cluster.varnodes.get(0).isInGene()) {
						buf.setColor(Color.lightGray);
						buf.drawString(cluster.varnodes.get(0).getTranscripts().get(0).getGenename() +" (Intronic)", (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
					else {
						buf.setColor(Color.gray);
						if(!cluster.varnodes.get(0).getTranscripts().get(0).equals(cluster.varnodes.get(0).getTranscripts().get(1))) {
							buf.drawString(cluster.varnodes.get(0).getTranscripts().get(0).getGenename() +" ... " +cluster.varnodes.get(0).getTranscripts().get(1).getGenename(), (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
							
						}
						else {
							if(cluster.varnodes.get(0).getTranscripts().get(0).getEnd() > cluster.varnodes.get(0).getPosition()) {
								buf.drawString(" ... " +cluster.varnodes.get(0).getTranscripts().get(0).getGenename(), (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
							}
							else {
								buf.drawString(cluster.varnodes.get(0).getTranscripts().get(0).getGenename() +" ... ", (int)header.get(c)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
								
							}
						}
					}					
				}
				else if (c == 6) {
					
					
					
					if(cluster.varnodes.get(0).getBedHits() != null) {
						bedarray = MethodLibrary.makeTrackArray(cluster.varnodes.get(0),null);
						for(int b = 0 ; b<bedarray.length; b++) {
							buf.setColor(Color.black);
							if(b == bedarray.length-1) {
								buf.fillRect((int)header.get(c+b)[1]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth()-(int)header.get(c+b)[1], rowHeight-1);	
							}
							else {
								buf.fillRect((int)header.get(c+b)[1]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, (int)header.get(c+b)[2], rowHeight-1);	
								
							}
							buf.setColor(Color.white);
							if(bedarray[b] != null) {
								buf.drawString(bedarray[b].toString(), (int)header.get(c+b)[1]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
							
							}							
						}
					}					
				}
				if(c < header.size()-1-Main.bedCanvas.bedTrack.size()) {
					buf.setColor(Color.black);
					buf.fillRect((int)header.get(c+1)[1]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, (int)header.get(c+1)[2], rowHeight-1);	
								
				}
			}		
			buf.setColor(Color.darkGray);
			buf.drawLine(3, rowHeight+3, 3, (rowHeight*(i+genemutcount+2))-tablescroll.getVerticalScrollBar().getValue()+3);	
			
			for(int r = 0; r <header.size();r++) {					
				buf.drawLine((int)header.get(r)[1], (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, (int)header.get(r)[1], (rowHeight*(i+genemutcount+2))-tablescroll.getVerticalScrollBar().getValue()+3);	
			}
			if(selectedNode != null && selectedNode.equals(cluster)) {
				
				hoverSample = -1;
				genemutcount = aminoarray.size()+1;						
				listAdd = 1;
				buf.drawLine(20, (rowHeight*(i+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(i+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+3);	
				drawGeneheader((rowHeight*(i+listAdd+1))-tablescroll.getVerticalScrollBar().getValue()+3);
				
				for(int s = 0; s<aminoarray.size(); s++) {
					
					buf.setColor(Color.darkGray);
					buf.drawLine(21, (rowHeight*(i+s+listAdd+3))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(i+s+listAdd+3))-tablescroll.getVerticalScrollBar().getValue()+3);	
					
					if(mouseY >= (rowHeight*(i+s+listAdd+2)) && mouseY < (rowHeight*(i+s+listAdd+3))) {
						hoverNode = null;
						hoverVar = aminoarray.get(s).getNode();		
						hoverString = aminoarray.get(s).getRow();
						buf.setColor(Color.white);
			
					for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
							if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(aminoarray.get(s).getRow()[5])) {
								hoverSample = aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getSample().getIndex();				
								break;
							}
						}									
					}	
					else {			
						if(MethodLibrary.aminoEffect(aminoarray.get(s).getRow()[3]).equals("nonsense")) {
							buf.setColor(Color.red);
						}
						else if (MethodLibrary.aminoEffect(aminoarray.get(s).getRow()[3]).equals("missense")) {
							buf.setColor(Color.yellow);
						}
						else if (MethodLibrary.aminoEffect(aminoarray.get(s).getRow()[3]).equals("synonymous")) {
							buf.setColor(Color.green);
						}
						else if(MethodLibrary.aminoEffect(aminoarray.get(s).getRow()[3]).equals("UTR")) {
							buf.setColor(Color.lightGray);
						}
						else {
							buf.setColor(Color.gray);
						}
					}
					if(!aminoarray.get(s).getRow()[1].equals("1")) {
						
						buf.drawString("Multiple", 24, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
					}
					else {			
						for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
							
							if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(aminoarray.get(s).getRow()[5])) {
								
								buf.drawString(aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getSample().getName(), 24, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
								break;
							}
						}
							
					}						
				
					
						if(hoverVar != null && hoverString.equals(aminoarray.get(s).getRow()) ) {
							//TODO
							textcolor = Color.white;
							
						}
						else {
							if(MethodLibrary.aminoEffect(aminoarray.get(s).getRow()[3]).equals("nonsense")) {
								textcolor = Color.red;
							}
							else if (MethodLibrary.aminoEffect(aminoarray.get(s).getRow()[3]).equals("missense")) {
								textcolor = Color.yellow;
							}
							else if (MethodLibrary.aminoEffect(aminoarray.get(s).getRow()[3]).equals("synonymous")) {
								textcolor = Color.green;
							}
							else if(aminoarray.get(s).getRow()[3].contains("UTR")) {
								textcolor = Color.lightGray;
							}
							else {
								textcolor = Color.gray;
							}
						}		
						for(int h = 1; h< 4; h++) {
							buf.setColor(Color.black);
							buf.fillRect((int)geneheader.get(h)[1]+10, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(h)[2], rowHeight-1);	
							buf.setColor(textcolor);
							buf.drawString(aminoarray.get(s).getRow()[h],  (int)geneheader.get(h)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
						
						}							
						
						
						if(aminoarray.get(s).getRow()[1].equals("1")) {
							buf.setColor(Color.black);
							buf.fillRect((int)geneheader.get(4)[1]+10, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
							buf.setColor(textcolor);
							
							for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
								if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(aminoarray.get(s).getRow()[5])) {
									if(aminoarray.get(s).getNode().vars.get(v).getValue().get(0).isHomozygous()) {
									buf.drawString("Hom (" +aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getCalls() +"/" +aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getCoverage() +")", (int)geneheader.get(4)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
									if(Control.controlData.controlsOn) {
										cases = 2;
										casefreq = 2/(double)(Main.varsamples*2);
									}
								}
								else {
									buf.drawString("Het ("+aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getCalls() +"/" +aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getCoverage() +")", (int)geneheader.get(4)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
									if(Control.controlData.controlsOn) {
										cases = 1;
										casefreq = 1/(double)(Main.varsamples*2);
									}
								
								}
									buf.setColor(Color.black);
									buf.fillRect((int)geneheader.get(5)[1]+1, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
									buf.setColor(textcolor);
									buf.drawString(""+aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getQuality(), (int)geneheader.get(5)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
								
								}
							}
						}
						
						if(Control.controlData.controlsOn) {
							cases = 0;
							for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
								if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(aminoarray.get(s).getRow()[5])) {
									if(aminoarray.get(s).getNode().vars.get(v).getValue().get(0).isHomozygous()) {											
										cases += Integer.parseInt(aminoarray.get(s).getRow()[1])*2;
									}
									else {												
										cases += Integer.parseInt(aminoarray.get(s).getRow()[1]);																						
									}	
								}																				
							}
							casefreq = cases/(double)(Main.varsamples*2);									
						}								
						
						buf.setColor(textcolor);
						buf.drawString(aminoarray.get(s).getRow()[4],  (int)geneheader.get(6)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
													
						if(Control.controlData.controlsOn) {
							buf.setColor(textcolor);
							
							for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
								if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(aminoarray.get(s).getRow()[5])) {
									vararray = aminoarray.get(s).getNode().vars.get(v).getValue();
									controlarray = new SampleNode[Control.controlData.fileArray.size()];
									if(vararray.get(vararray.size()-1).alleles != null) {											
										
										for(int e = vararray.size()-1; e> 0;e-- ) {
											
											if(vararray.get(e).alleles == null) {													
												break;
											}
											controlarray[vararray.get(e).getControlSample().getIndex()] = vararray.get(e);
									
										}										
									}
									for(int e = 0; e<controlarray.length;e++ ) {
										if(Control.controlData.fileArray.get(e).controlOn) {
											if(controlarray[e] == null) {
												buf.setColor(Color.black);
												buf.fillRect((int)geneheader.get(this.geneheaderlength+e*2)[1]+11, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
												buf.setColor(textcolor);											
												buf.drawString("0", (int)geneheader.get(this.geneheaderlength+e*2)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
												buf.setColor(Color.black);
												buf.fillRect((int)geneheader.get(this.geneheaderlength+e*2+1)[1]+11, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
												buf.setColor(textcolor);
												buf.drawString("-", (int)geneheader.get(this.geneheaderlength+e*2+1)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
										
											}
											else {
												buf.setColor(Color.black);
												buf.fillRect((int)geneheader.get(this.geneheaderlength+e*2)[1]+11, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
												buf.setColor(textcolor);											
												buf.drawString(""+MethodLibrary.round(controlarray[e].alleles/(double)controlarray[e].allelenumber,2), (int)geneheader.get(this.geneheaderlength+e*2)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
												buf.setColor(Color.black);
												buf.fillRect((int)geneheader.get(this.geneheaderlength+e*2+1)[1]+11, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
												buf.setColor(textcolor);
												buf.drawString(""+MethodLibrary.round(casefreq/(controlarray[e].alleles/(double)controlarray[e].allelenumber),2) +" (p=" +MethodLibrary.round(VariantHandler.table.fe.getRightTailedP(cases, Main.varsamples*2-cases, controlarray[e].alleles, controlarray[e].allelenumber-controlarray[e].alleles) ,2) +")", (int)geneheader.get(this.geneheaderlength+e*2+1)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
											}													
										}
										else {
											buf.setColor(Color.black);
											buf.fillRect((int)geneheader.get(this.geneheaderlength+e*2)[1]+11, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
											buf.setColor(Color.darkGray);											
											buf.drawString("Apply controls", (int)geneheader.get(this.geneheaderlength+e*2)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
											buf.setColor(Color.black);
											buf.fillRect((int)geneheader.get(this.geneheaderlength+e*2+1)[1]+11, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
											buf.setColor(Color.darkGray);			
											buf.drawString("-", (int)geneheader.get(this.geneheaderlength+e*2+1)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
									
										}									
									}									
								}
							}							
						}
						else {
							buf.setColor(Color.darkGray);
							
							for(int e = geneheaderlength; e<geneheader.size();e++ ) {
								if(geneheader.get(e)[0] instanceof ControlFile) {
								buf.drawString("Apply controls",  (int)geneheader.get(e)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
								}
							}
							buf.setColor(Color.lightGray);
						}
						vararray = null;
						if(Main.bedCanvas.bedOn) {								
							for(int a =0;a<aminoarray.size(); a++) {									
								bedarray = MethodLibrary.makeTrackArray(aminoarray.get(a).getNode(),aminoarray.get(a).getRow()[5]);
								if(bedarray != null) {
								for(int b = 0 ; b<bedarray.length; b++) {
										buf.setColor(Color.black);
										if(b == bedarray.length-1) {
											buf.fillRect((int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+b)[1]+12, (rowHeight*(i+a+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth()-(int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+b)[1], rowHeight-1);	
										}
										else {
											buf.fillRect((int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+b)[1]+12, (rowHeight*(i+a+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+b)[2], rowHeight-1);												
										}
										buf.setColor(Color.white);
										if(bedarray[b] != null) {
											buf.drawString(bedarray[b].toString(), (int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+b)[1]+14, (rowHeight*(i+a+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
										
										}											
									}
								}									
							}
						}
						buf.setColor(Color.darkGray);
						for(int j = 0; j<geneheader.size(); j++) {								
							
							buf.drawLine((int)geneheader.get(j)[1]+11, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(j)[1]+11, (rowHeight*(i+s+listAdd+3))-tablescroll.getVerticalScrollBar().getValue()+3);	
						}
					if(selectedVar!=null && selectedString.equals(aminoarray.get(s).getRow()) && Integer.parseInt(selectedString[1]) > 1) {
						//hoverSample = -1;
						pointer = 0;
						//TODO
						
						for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
								if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(selectedString[5])) {
									
									for(int l = 0; l<aminoarray.get(s).getNode().vars.get(v).getValue().size();l++) {
										if(aminoarray.get(s).getNode().vars.get(v).getValue().get(l).alleles != null) {
											break;
										}
									if(mouseY > (rowHeight*(i+s+pointer+4)) && mouseY < (rowHeight*(i+s+pointer+5))) {
										textcolor = Color.white;
										
										hoverVar = aminoarray.get(s).getNode();		
										hoverString = aminoarray.get(s).getRow();
										hoverSample = aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getSample().getIndex();
									}
									else {
										textcolor = Color.lightGray;
									}
									
								//	if(aminoarray.get(s).getNode().getSamples().get(l).getVariation().equals(selectedString[5])) {									
										buf.setColor(textcolor);
										buf.drawString(aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getSample().getName(), 30, (rowHeight*(i+s+pointer+4))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
										pointer++;
								//	}	
									
										buf.setColor(Color.black);
										buf.fillRect((int)geneheader.get(4)[1]+10, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
										buf.setColor(textcolor);
									if(aminoarray.get(s).getNode().vars.get(v).getValue().get(l).isHomozygous()) {
										buf.drawString("Hom (" +aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getCalls() +"/" +aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getCoverage() +")",  (int)geneheader.get(4)[1]+14, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
										
									}
									else {
										buf.drawString("Het ("+aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getCalls() +"/" +aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getCoverage() +")",  (int)geneheader.get(4)[1]+14, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
									}
									
										buf.setColor(Color.black);
										buf.fillRect((int)geneheader.get(5)[1]+10, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
										buf.setColor(textcolor);
										buf.drawString(""+aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getQuality(),  (int)geneheader.get(5)[1]+14, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
										buf.setColor(Color.darkGray);
										for(int j = 4; j<7; j++) {								
											
											buf.drawLine((int)geneheader.get(j)[1]+11, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue(), (int)geneheader.get(j)[1]+11, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+rowHeight+2);	
										}
									
									}
							}
						}
						listAdd = Integer.parseInt(selectedString[1])+1;
						genemutcount = aminoarray.size() +listAdd;
						buf.setColor(Color.darkGray);
						buf.drawLine(21, (rowHeight*(i+s+listAdd+3))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(i+s+listAdd+3))-tablescroll.getVerticalScrollBar().getValue()+3);	
					
						
					}
					
				
				}
			}
			
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
		}	
		
	}
	buf.setColor(Color.darkGray);
	buf.drawLine(4, (rowHeight*(Main.drawCanvas.clusterNodes.size()+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(Main.drawCanvas.clusterNodes.size()+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3);	

}
drawHeader();

if(headerHover == -1 && geneHeaderHover == -1) {		
	
	setCursor( Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));			
}
else {
	if(resizeColumn == -1) {
		setCursor( Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));		
	}
	else {
		setCursor( Cursor.getPredefinedCursor(Cursor.E_RESIZE_CURSOR));		
	}
}


g.drawImage(bufImage, 0, tablescroll.getVerticalScrollBar().getValue(), null);

}
catch(Exception e) {
	e.printStackTrace();
}
}	
void drawHeader() { 
	for(int i = 0; i<header.size(); i++) {
		buf.setColor(Color.darkGray);
		buf.fillRect((int)header.get(i)[1], 0, (int)(header.get(i)[2])+1, rowHeight);
		
		if(mouseY-tablescroll.getVerticalScrollBar().getValue() <= rowHeight) {
			if(mouseX >= (int)((int)header.get(i)[1]) && mouseX <= (int)((int)header.get(i)[1])+(int)(header.get(i)[2])) {
				headerHover = i;
				buf.setColor(Color.yellow);			
			}
			else if( mouseX >= (int)(header.get(header.size()-1)[1])+(int)(header.get(header.size()-1)[2])) {
				headerHover = header.size();
				buf.setColor(Color.white);			
			}
			else {
				buf.setColor(Color.white);		
			}		
		}
		else {
			
			buf.setColor(Color.white);
		}
		if(!mouseDrag && headerHover > -1 && i > 0 && i < header.size()) {
			if(mouseX > (int)((int)header.get(i)[1])-5 && mouseX < (int)((int)header.get(i)[1])+5) {
				resizeColumn = i;
				
			}
		}
		if(header.get(i)[0] instanceof String) {
			buf.drawString((String)header.get(i)[0], (int)((int)header.get(i)[1])+4, rowHeight-2);
		}
		else if (header.get(i)[0] instanceof BedTrack) {
			BedTrack track =  (BedTrack)header.get(i)[0];
			
			if(track.file != null) {
				buf.drawString(track.file.getName(), (int)((int)header.get(i)[1])+4, rowHeight-2);
			}
			else {
				buf.drawString(FilenameUtils.getName(track.url.getFile()),  (int)((int)header.get(i)[1])+4, rowHeight-2);
			}
			track = null;
		}
		
		buf.setColor(Color.black);
		buf.drawLine((int)((int)header.get(i)[1]), 0, (int)((int)header.get(i)[1]), rowHeight);
		if(sorter.index > -1) {
			if(sorter.ascending) {
				buf.setColor(Color.white);
				buf.fillPolygon(this.sortTriangle);
			}
			else {
				buf.setColor(Color.white);
				buf.fillPolygon(this.sortTriangle);
			}
		}
	}
	buf.setColor(Color.darkGray);
	buf.fillRect((int)header.get(header.size()-1)[1]+(int)header.get(header.size()-1)[2]+2, 0, this.width, rowHeight);
	buf.setColor(Color.black);
	buf.drawLine((int)header.get(header.size()-1)[1]+(int)header.get(header.size()-1)[2]+1,0, (int)header.get(header.size()-1)[1]+(int)header.get(header.size()-1)[2]+1, rowHeight);
	if(!mouseDrag && headerHover > -1 && resizeColumn == -1) {
		
		if(mouseX > ((int)header.get(header.size()-1)[1]+(int)header.get(header.size()-1)[2])-5 && mouseX < ((int)header.get(header.size()-1)[1]+(int)header.get(header.size()-1)[2])+5) {
			resizeColumn = header.size();
			
		}
	}
	
}
void drawGeneheader(int y) {
	for(int i = 0; i<geneheader.size(); i++) {
		buf.setColor(Color.darkGray);
		buf.fillRect((int)((int)geneheader.get(i)[1])+10, y, (int)geneheader.get(i)[2], rowHeight);
		
		if(mouseY-tablescroll.getVerticalScrollBar().getValue() <= y+rowHeight && mouseY-tablescroll.getVerticalScrollBar().getValue() >= y) {
			
			if(mouseX >= (int)geneheader.get(i)[1]+10 && mouseX <= (int)geneheader.get(i)[1]+(int)geneheader.get(i)[2]+10) {				
				
				geneHeaderHover = i;				
				buf.setColor(Color.yellow);
			}
			else if(mouseX >= (int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2]+10) {
				geneHeaderHover = geneheader.size();
				buf.setColor(Color.white);
			}
			else {
				buf.setColor(Color.white);
			}
			
		}
		else {
			
			buf.setColor(Color.white);
		}
		if(!mouseDrag && geneHeaderHover > -1 && i > 0 && i < geneheader.size()) {
			if(mouseX > (int)geneheader.get(i)[1]+5 && mouseX < (int)geneheader.get(i)[1]+15) {
				if(resizeColumn != i) {
					resizeColumn = i;						
				}
			}
		}
		if(geneheader.get(i)[0] instanceof String) {
			buf.drawString((String)geneheader.get(i)[0], (int)geneheader.get(i)[1]+14, y+rowHeight-3);
		}
		else if (geneheader.get(i)[0] instanceof ControlFile) {
			ControlFile ctrlfile = (ControlFile)geneheader.get(i)[0];
			buf.drawString("AF: " +ctrlfile.getName(), (int)geneheader.get(i)[1]+14, y+rowHeight-3);
			ctrlfile = null;
		}
		else  {
			BedTrack track = (BedTrack)geneheader.get(i)[0];
			buf.drawString(track.file.getName(), (int)geneheader.get(i)[1]+14, y+rowHeight-3);
			track = null;
		}
		buf.setColor(Color.black);
		buf.drawLine((int)geneheader.get(i)[1]+10, y, (int)((int)geneheader.get(i)[1])+10, y+rowHeight);
		
	}

	if(!mouseDrag && geneHeaderHover > -1 && resizeColumn == -1) {		
		if(mouseX > (int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2]-5 && mouseX < (int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2]+15) {
			resizeColumn = geneheader.size();			
		}
	}
	buf.setColor(Color.darkGray);
	buf.fillRect((int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2]+11, y, this.width, rowHeight);
	buf.setColor(Color.black);
	buf.drawLine((int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2]+10, y, (int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2]+10, y+rowHeight);
}



void clear() {
	samples.clear();
	Main.drawCanvas.clusterNodes.clear();
	
	varAdd = null;
	hoverNode = null;
	selectedNode = null;
	entry = null;
	vararray.clear();
	hoverVar = null;
	selectedVar = null;

}
int getTableSize() {
	if(Main.drawCanvas.clusterNodes == null) {
		return 0;
	}
	return Main.drawCanvas.clusterNodes.size();
}

void createPolygon() {
	if(sorter.index > header.size()-1) {
		return;
	}
	if(sorter.ascending) {
		int[] x = { (int)header.get(sorter.index)[1]+(int)(header.get(sorter.index)[2])-15, (int)header.get(sorter.index)[1]+(int)(header.get(sorter.index)[2])-10, (int)header.get(sorter.index)[1]+(int)(header.get(sorter.index)[2])-5};
		int[] y = { 12, 4, 12 };
		int n = 3;
		this.sortTriangle = new Polygon(x,y,n);
	}
	else {
		int[] x = { (int)header.get(sorter.index)[1]+(int)(header.get(sorter.index)[2])-15, (int)header.get(sorter.index)[1]+(int)(header.get(sorter.index)[2])-10, (int)header.get(sorter.index)[1]+(int)(header.get(sorter.index)[2])-5};
		int[] y = { 4, 12, 4 };
		int n = 3;
		this.sortTriangle = new Polygon(x,y,n);
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
public static class ListSorter implements Comparator<ClusterNode> {
	public int index;
	boolean ascending = true;
	
	
	public int compare(ClusterNode o1, ClusterNode o2) {  
		
		
		if(index < 3) {
			int f1, f2;
			if(index == 0) {
				f1 = o1.ID;
				f2 = o2.ID;
			}
			else if(index == 1) {
				f1 = o1.nodecount;
				f2 = o2.nodecount;
			}
			else {
				f1 = o1.width;
				f2 = o2.width;
			}
			 if ( f1 < f2) {  
		        	if(ascending) {
		        		return -1;  
		        	}
		        	else {
		        		return 1;
		        	}
		                
		        } else if ( f1 > f2 ) {  
		        	if(ascending) {
		                return 1;  
		        	}
		        	else {
		        		return -1;
		        	}
		        }
		        else {
		        	return 0;
		        }      
		}
		else if(index == 3) {
			String chrom1, chrom2;
			if(o1.varnodes.get(0).getExons() != null) {
				chrom1 = o1.varnodes.get(0).getExons().get(0).getTranscript().getChrom();
			}
			else {
				chrom1 = o1.varnodes.get(0).getTranscripts().get(0).getChrom();
			}
			if(o2.varnodes.get(0).getExons() != null) {
				chrom2 = o2.varnodes.get(0).getExons().get(0).getTranscript().getChrom();
			}
			else {
				chrom2 = o2.varnodes.get(0).getTranscripts().get(0).getChrom();
			}
			
			 if ( Main.chromIndex.get(chrom1)[0].longValue()+o1.varnodes.get(0).getPosition() < Main.chromIndex.get(chrom2)[0].longValue()+o2.varnodes.get(0).getPosition()) {  
		        	if(ascending) {
		        		return -1;  
		        	}
		        	else {
		        		return 1;
		        	}
		                
		        } else if ( Main.chromIndex.get(chrom1)[0].longValue()+o1.varnodes.get(0).getPosition() > Main.chromIndex.get(chrom2)[0].longValue()+o2.varnodes.get(0).getPosition()) {  
		        	if(ascending) {
		                return 1;  
		        	}
		        	else {
		        		return -1;
		        	}
		        }
		        else {
		        	return 0;
		        }  
			
			
		}
		else if(index == 4){
			double f1 = o1.nodecount/(double)o1.width, f2 = o2.nodecount/(double)o2.width;
			if ( f1 < f2) {  
	        	if(ascending) {
	        		return -1;  
	        	}
	        	else {
	        		return 1;
	        	}
	                
	        } else if ( f1 > f2 ) {  
	        	if(ascending) {
	                return 1;  
	        	}
	        	else {
	        		return -1;
	        	}
	        }
	        else {
	        	return 0;
	        }      
		}
		else {
			return 0;
		}
			
			
		
			
}  
}

public void mouseClicked(MouseEvent event) {
	switch(event.getModifiers()) {
		case InputEvent.BUTTON1_MASK: {	
			if(!this.isEnabled()) {
				break;
			}
			if(this.headerHover == -1) {
				if(event.getClickCount() == 2 ) {
					FileRead.novars = true;
							
					if(hoverSample > -1) {
						Main.drawCanvas.drawVariables.visiblestart = (short)(hoverSample);
						Main.drawCanvas.drawVariables.visiblesamples = (short)(1);						
						Main.drawCanvas.resizeCanvas(Main.drawScroll.getViewport().getWidth(), (int)(Main.samples*Main.drawCanvas.drawVariables.sampleHeight));					
						Draw.setScrollbar((int)(hoverSample*Main.drawCanvas.drawVariables.sampleHeight));	
					}
					FileRead.search = true;	
					if(hoverVar != null) {						
						
						VarNode searchHead = hoverVar;
						
						while(searchHead.getPrev() != null) {
							if(searchHead.getPrev().getPosition() == 0) {
								searchHead.getPrev().putNext(searchHead);
							}
							searchHead = searchHead.getPrev();		
						
						}
						
						FileRead.head = searchHead;						
						searchHead = null;
						Main.drawCanvas.current = hoverVar;	
									
						if(hoverVar.getExons() != null) {
							Main.drawCanvas.gotoPos(hoverVar.getExons().get(0).getTranscript().getChrom(), hoverVar.getPosition()+1-50, hoverVar.getPosition()+1+50);
						}
						else {
							if(hoverVar.getTranscripts() != null) {
								Main.drawCanvas.gotoPos(hoverVar.getTranscripts().get(0).getChrom(), hoverVar.getPosition()+1-50, hoverVar.getPosition()+1+50);
							}
						}						
					}
					else if(hoverNode != null) {						
						VarNode searchHead = hoverNode.varnodes.get(0);					
						if(hoverNode.varnodes.get(0).getTranscripts() != null) {
							if(!chrom.equals(hoverNode.varnodes.get(0).getTranscripts().get(0).getChrom())) {
								chrom = hoverNode.varnodes.get(0).getTranscripts().get(0).getChrom();
							}
						}
						else {
							if(!chrom.equals(hoverNode.varnodes.get(0).getExons().get(0).transcript.getChrom())) {
								chrom = hoverNode.varnodes.get(0).getExons().get(0).transcript.getChrom();
							}
						}
						while(searchHead.getPrev() != null) {
							if(searchHead.getPrev().getPosition() == 0) {
								searchHead.getPrev().putNext(searchHead);
							}
							searchHead = searchHead.getPrev();							
						}
						
						FileRead.head = searchHead;					
						searchHead = null;
						Main.drawCanvas.current = hoverNode.varnodes.get(0);
						Main.drawCanvas.gotoPos(chrom, hoverNode.varnodes.get(0).getPosition(),  hoverNode.varnodes.get(0).getPosition()+hoverNode.width);						
					}
				
					break;
				}
				else if(event.getClickCount() == 1 ) {
					
					if(hoverVar != null && (selectedVar == null || !selectedVar.equals(hoverVar) )){
						
						selectedVar = hoverVar;
						selectedString = hoverString;
					
						if(selectedVar.isRscode() !=null) {
							hoverString[4] = selectedVar.rscode;
						}
						repaint();
					
						this.setPreferredSize(new Dimension(this.getWidth(), (this.getTableSize()+2+samplecount+Integer.parseInt(selectedString[1]))*rowHeight));
						
						this.revalidate();
						
					}
					else if(hoverVar != null && selectedVar.equals(hoverVar)){
					
						if(hoverSample == -1) {
							
						
							selectedVar = null;
							
							repaint();
							this.setPreferredSize(new Dimension(this.getWidth(), (this.getTableSize()+2+samplecount)*rowHeight));
							
							this.revalidate();
						}
					}
					else {					
					
						if(hoverVar != null && selectedNode != null && hoverNode.equals(selectedNode)) {
							
							selectedString = hoverString;
							
							samplecount = 0;
							repaint();
						
						}
						else {
							
							if(hoverSample == -1) {
								
								if(hoverNode != null) {
									
									if(selectedNode != null && selectedNode.equals(hoverNode)) {
										
										selectedString = null;
										selectedNode = null;
										hoverVar = null;
										selectedVar = null;
										aminoarray.clear();
										this.setPreferredSize(new Dimension(this.getWidth(), (this.getTableSize()+2)*rowHeight+10));
										
										this.revalidate();
										
										repaint();
										break;
									}
									
									selectedNode = hoverNode;
									selectedString = hoverString;
								//	samplecount = selectedNode.mutations;
									samplecount = selectedNode.varnodes.size();
														
									this.setPreferredSize(new Dimension(this.getWidth(), (this.getTableSize()+2+samplecount)*rowHeight+10));
									this.revalidate();										
									getAminos(selectedNode);						
									repaint();
																		
								}	
							}
							
						}						
					}
				}
			}
			break;
			
		}
		case InputEvent.BUTTON3_MASK: {		
			
			selectedNode = null;
			selectedVar = null;
			repaint();
		}
	}
}
void getAminos(ClusterNode transcript) {
	
	try {
	aminoarray.clear();
	VarNode varnode = null;
	
	for(int t = 0; t<transcript.varnodes.size(); t++) {
	
			if(Main.drawCanvas.hideNode(transcript.varnodes.get(t))) {
				
				continue;
			}
		
			varnode = transcript.varnodes.get(t);
			
			
			for(int v = 0; v<varnode.vars.size(); v++) { // Map.Entry<String, ArrayList<SampleNode>> entry : varnode.vars) {
					entry = varnode.vars.get(v);
					mutcount = 0;
					if(Main.drawCanvas.hideNodeVar(varnode, entry)) {
						continue;
					}
					for(int m = 0; m<entry.getValue().size(); m++) {
						if(entry.getValue().get(m).alleles != null) {
							break;
						}
						if(!Main.drawCanvas.hideVar(entry.getValue().get(m))) {
							mutcount++;
						}
						else {
							entry.getValue().remove(m);
							m--;
						}
					}
					
					if(mutcount == 0) {
						continue;
					}
					if(varnode.getTranscripts() != null) {
						if(!chrom.equals(varnode.getTranscripts().get(0).getChrom())) {
							chrom = varnode.getTranscripts().get(0).getChrom();
						}
					}
					else {
						if(!chrom.equals(varnode.getExons().get(0).transcript.getChrom())) {
							chrom = varnode.getExons().get(0).transcript.getChrom();
						}
					}
					base = entry.getKey();
					String[] addrow = new String[6];						
					addrow[0] = ""+transcript.ID;						
					addrow[1] = ""+mutcount;
					addrow[2] = chrom +":"+MethodLibrary.formatNumber((varnode.getPosition()+1));
					if(base.length() == 1) {
						addrow[3] = "" +Main.getBase.get(varnode.getRefBase()) +">" +base;							
					}
					else {
						addrow[3] = base;	
					}
					if(varnode.isRscode() != null) {
						addrow[4] = varnode.rscode;
					}
					else {
						addrow[4] = "N/A";
					}
					addrow[5] = base;				
				
					AminoEntry aminoentry = new AminoEntry(addrow, varnode);
				
					aminoarray.add(aminoentry);					
				
				//	aminoarray.add(new Entry(vardraw.getExons().get(0).getTranscript().getGenename() +" " +vardraw.getPosition() +" " +baseTemp.get(base) +" " +Main.chromDraw.getAminoChange(vardraw,base,vardraw.getExons().get(exon), true), vardraw));
											
				
			}
			
		

	
	}
	varnode = null;
	}
	catch(Exception e) {
		ErrorLog.addError(e.getStackTrace());
		e.printStackTrace();
	}
	
}
@Override
public void mouseEntered(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}


@Override
public void mouseExited(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}




public void mousePressed(MouseEvent event) {
	
	switch(event.getModifiers()) {
	
		case InputEvent.BUTTON1_MASK: {		
			if(!this.isEnabled()) {
				break;
			}
			this.dragX = event.getX();
			if(headerHover > -1) {
				if(resizeColumn == -1) {
					if(sorter.ascending) {
						sorter.ascending = false;
					}
					else {
						sorter.ascending = true;
					}
					
					sorter.index = headerHover;
					Collections.sort(Main.drawCanvas.clusterNodes, sorter);
					createPolygon();
					
					repaint();
				}
			}
			
		}

		if(hoverNode != null || hoverVar != null) {
			Main.chromDraw.repaint();
		}
	}
	
}


@Override
public void mouseReleased(MouseEvent arg0) {
	mouseDrag = false;
}


@Override
public void mouseDragged(MouseEvent event) {
	mouseDrag = true;
	if(resizeColumn > 0) {
		
		resizeTable(resizeColumn, event.getX()-this.dragX);
		this.dragX = event.getX();	
	}
}

@Override
public void mouseMoved(MouseEvent event) {
	if(!this.isEnabled()) {
		return;
	}
	mouseY = event.getY();
	mouseX = event.getX();
	
	repaint();
}



@Override
public void mouseWheelMoved(MouseWheelEvent e) {
	if ( e.getWheelRotation() < 0 ) {	
		tablescroll.getVerticalScrollBar().setValue(tablescroll.getVerticalScrollBar().getValue()-16);
		
	}
	else {			
		tablescroll.getVerticalScrollBar().setValue(tablescroll.getVerticalScrollBar().getValue()+16);
		
	}
	repaint();
	
}

}
