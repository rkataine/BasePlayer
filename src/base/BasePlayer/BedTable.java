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


public class BedTable extends JPanel implements MouseMotionListener, MouseListener, MouseWheelListener {

	
		private static final long serialVersionUID = 1L;	
		BufferedImage bufImage;
		ListSorter sorter = new ListSorter();
		int rowHeight = 15, geneheaderlength = 0;
		Graphics2D buf;
		StringBuffer mutcountbuffer = new StringBuffer("");
		int width, height;
		ArrayList<AminoEntry> aminoarray = new ArrayList<AminoEntry>();
		ArrayList<BedNode> bedarray = new ArrayList<BedNode>();
		int mouseY=0, mouseX=0, pressX=0, pressY=0;
		final JScrollPane tablescroll;
		public BedNode hoverNode, selectedNode;
		String[] posSplit, hoverString;
		int samplecount = 0, variants = 0;
		Enumeration<String> e;
		String base;
	//	private FontMetrics fm;
		final BedTrack bedtrack;
	//	ArrayList<String> geneheader =  new ArrayList<String>();
		ArrayList<Object[]> geneheader =  new ArrayList<Object[]>();
		String[] header = {"Event", "Mut. count", "Event position", "Event length (bp)", "Mutation freq. (%)", "Flanking genes" };
		int[][] headerlengths = new int[header.length][2];
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
		private Color textcolor;
		private double casefreq;		
		Map.Entry<String, ArrayList<SampleNode>>  entry;
		//private int textWidth = 0;
		MethodLibrary.controlsorter ctrlsort = new MethodLibrary.controlsorter();
		private int geneHeaderHover;
		private boolean mouseDrag;
		private int resizeColumn;
		private int cases;
		private int dragX;
		private int firstrow = 0;
		private int firstvisible;
		
		BedTable(int width, int height, JScrollPane tablescroll, BedTrack bedtrack) {			
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
			/* for(int i = 0 ; i<Control.controlData.fileArray.size(); i++) {							
				  addRowGeneheader("AF: " +Control.controlData.fileArray.get(i).getName());
				  addRowGeneheader("OR");						 
			  }*/
			bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width, (int)height, BufferedImage.TYPE_INT_ARGB));	
			buf = (Graphics2D)bufImage.getGraphics();
			this.addMouseListener(this);
			this.addMouseMotionListener(this);
			this.controlarray = VariantHandler.table.controlarray;
	//		fm = buf.getFontMetrics();
			this.bedtrack = bedtrack;
			headerlengths[0][0] = 0;
			headerlengths[0][1] = (int)(0.15*width);
			headerlengths[1][0] = headerlengths[0][1]+headerlengths[0][0];	
			headerlengths[1][1] = (int)(0.15*width);
			headerlengths[2][0] = headerlengths[1][1]+headerlengths[1][0];
			headerlengths[2][1] = (int)(0.2*width);
			headerlengths[3][0] = headerlengths[2][1]+headerlengths[2][0];
			headerlengths[3][1] = (int)(0.5*width);
			
			this.addMouseWheelListener(this); /*
					new MouseWheelListener() {
						
						public void mouseWheelMoved(MouseWheelEvent mwe) {
							
							if ( mwe.getWheelRotation() < 0 ) {	
								tablescroll.getVerticalScrollBar().setValue(tablescroll.getVerticalScrollBar().getValue()-16);
								
							}
							else {			
								tablescroll.getVerticalScrollBar().setValue(tablescroll.getVerticalScrollBar().getValue()+16);
								
							}
						}
					}
					
				);	*/
		}	
		
void resizeTable() {
			
	if(bufImage.getWidth() < headerlengths[headerlengths.length-1][0]+headerlengths[headerlengths.length-1][1]) {
		bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
		buf = (Graphics2D)bufImage.getGraphics();
	}
	
}
void resizeTable(int width) {
	headerlengths[0][0] = 0;
	headerlengths[0][1] = (int)(width/(double)headerlengths.length);
	for(int i = 1 ; i<headerlengths.length; i++) {
		headerlengths[i][0] = headerlengths[i-1][1]+headerlengths[i-1][0];	
		headerlengths[i][1] = (int)(width/(double)headerlengths.length);		
	}
	
	geneheader.get(0)[2] = (int)((width-10)/geneheader.size());
	
	for(int i = 1; i<geneheader.size(); i++) {
		geneheader.get(i)[1] = (int)geneheader.get(i-1)[1] + (int)geneheader.get(i-1)[2];
		geneheader.get(i)[2] = (int)((width-10)/geneheader.size());
	}
	this.setPreferredSize(new Dimension(width, this.getHeight()));
	this.revalidate();
}
void addRowGeneheader(Object column) {
	
	Object[] obj = new Object[3]; 
	obj[0] = column; 	
	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2];
	obj[2] = 100;
	geneheader.add(obj);
	
	if((int)obj[1]+(int)obj[2] > this.getWidth()) {
		if(bufImage == null) {
			bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width, (int)height, BufferedImage.TYPE_INT_ARGB));	
			buf = (Graphics2D)bufImage.getGraphics();
		}
		if(bufImage.getWidth() < (int)obj[1]+(int)obj[2]) {
			bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
			buf = (Graphics2D)bufImage.getGraphics();
		}
		this.setPreferredSize(new Dimension((int)obj[1]+(int)obj[2], this.getHeight()));
		this.revalidate();
	}
}



void drawScreen(Graphics g) {	
	try {
	buf.setColor(Color.black);
	
	buf.fillRect(0, 0, this.getWidth(), this.getHeight());	
	if(width != this.getWidth()) {		
		width = this.getWidth();		
		createPolygon();
		resizeTable();
	}


	genemutcount = 0;
	if(!bedtrack.intersect) {
		buf.setColor(Color.white);
		buf.drawString("Press play on bed track to annotate variants", 5, 40);
	}
	else if(getTableSize() > 0) {		
		
		hoverVar = null;
		hoverSample = -1;
		headerHover = -1;
		geneHeaderHover = -1;
		if(!mouseDrag) {
			resizeColumn = -1;
		}
		
		if(aminoarray == null) {
			aminoarray = new ArrayList<AminoEntry>();
		}
		
		firstrow = tablescroll.getVerticalScrollBar().getValue()/rowHeight -samplecount-listAdd-aminoarray.size();
	
		if(firstrow < 0) {
			firstrow = 0;
		}
		for(int i = firstrow; i<bedarray.size(); i++) {
			dot = false;			
			
			if((i+1+samplecount+listAdd+aminoarray.size())*rowHeight < tablescroll.getVerticalScrollBar().getValue()) {
				
				continue;
			}
			
			if(i*rowHeight > tablescroll.getVerticalScrollBar().getValue() + tablescroll.getViewport().getHeight()) {
				
				break;
			}
			
			if(mouseY >= (rowHeight*(i+genemutcount+1)) && mouseY < (rowHeight*(i+genemutcount+2))) {
				hoverNode = bedarray.get(i);	
			}			
			
			try {
				buf.setColor(Color.darkGray);
				buf.drawLine(4, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3);	
				
				
				if(bedarray.get(i).equals(hoverNode) || bedarray.get(i).equals(selectedNode)) {
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				
				if(bedarray.get(i).getTrack().hasvalues) {
					buf.drawString((i+1) +".  " +MethodLibrary.shortString(bedarray.get(i).name,10) +"=" +MethodLibrary.round(bedarray.get(i).value,3), 5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}
				else {
					buf.drawString((i+1) +".  " +MethodLibrary.shortString(bedarray.get(i).name,10), 5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					
				}
				buf.setColor(Color.black);
				buf.fillRect(headerlengths[1][0]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, headerlengths[1][1], rowHeight-1);	
				
				if(bedarray.get(i).equals(hoverNode) || bedarray.get(i).equals(selectedNode)) {
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				mutcountbuffer = new StringBuffer(""+bedarray.get(i).mutations +" ");
				buf.drawString(mutcountbuffer.toString(),  headerlengths[1][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
			//		lastpos = Integer.toString(bedarray.get(i).mutations).length() +2;
		//TODO		textWidth = (int)fm.getStringBounds("", buf).getWidth();
				
				
			//	textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
			//	buf.drawString("  ", headerlengths[1][0]+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				
			//	buf.setColor(Color.gray);
			//	textWidth = (int)fm.getStringBounds(mutcountbuffer.toString() , buf).getWidth();
			//	buf.drawString(" " +bedarray.get(i).varnodes.size() +" samples",  headerlengths[1][0]+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				
				buf.setColor(Color.black);
				buf.fillRect(headerlengths[2][0]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
		
				if(bedarray.get(i).equals(hoverNode) || bedarray.get(i).equals(selectedNode)) {
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				
				buf.drawString(bedarray.get(i).getChrom() +":"+MethodLibrary.formatNumber(bedarray.get(i).getPosition()+1) +"-" +MethodLibrary.formatNumber(bedarray.get(i).getPosition()+1+bedarray.get(i).getLength()) , headerlengths[2][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				buf.setColor(Color.black);
				buf.fillRect(headerlengths[3][0]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
				
				if(bedarray.get(i).equals(hoverNode) || bedarray.get(i).equals(selectedNode)) {
					
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				
				buf.drawString(MethodLibrary.formatNumber(bedarray.get(i).getLength()) , headerlengths[3][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				buf.setColor(Color.black);
				buf.fillRect(headerlengths[4][0]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
				
				if(bedarray.get(i).equals(hoverNode) || bedarray.get(i).equals(selectedNode)) {
					
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				
				buf.drawString(""+MethodLibrary.round((bedarray.get(i).mutations/(double)bedarray.get(i).getLength())*100,4), headerlengths[4][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				buf.setColor(Color.black);
				buf.fillRect(headerlengths[5][0]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
				
				if(bedarray.get(i).equals(hoverNode) || bedarray.get(i).equals(selectedNode)) {
					
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				firstvisible = 0;
				if(bedarray.get(i).varnodes != null) {
					
				
					for(int f = 0 ;f<bedarray.get(i).varnodes.size(); f++) {
						if(!Main.drawCanvas.hideNode(bedarray.get(i).varnodes.get(f))) {
							firstvisible = f;
							break;
						}
					}
					if(bedarray.get(i).varnodes.get(firstvisible).getExons() != null) {
						
						if(bedarray.get(i).varnodes.get(firstvisible).coding) {
							buf.setColor(Color.red);
							buf.drawString(bedarray.get(i).varnodes.get(firstvisible).getExons().get(0).getTranscript().getGenename() +" (Coding)", headerlengths[5][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
							
						}
						else {
							buf.setColor(Color.lightGray);
							buf.drawString(bedarray.get(i).varnodes.get(firstvisible).getExons().get(0).getTranscript().getGenename() +" (UTR)", headerlengths[5][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
							
						}
					}
					else if(bedarray.get(i).varnodes.get(firstvisible).isInGene()) {
						
						buf.setColor(Color.lightGray);
						buf.drawString(bedarray.get(i).varnodes.get(firstvisible).getTranscripts().get(0).getGenename() +" (Intronic)", headerlengths[5][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
					else {
						buf.setColor(Color.gray);
						
						if(!bedarray.get(i).varnodes.get(firstvisible).getTranscripts().get(0).equals(bedarray.get(i).varnodes.get(firstvisible).getTranscripts().get(1))) {
							
							buf.drawString(bedarray.get(i).varnodes.get(firstvisible).getTranscripts().get(0).getGenename() +" ... " +bedarray.get(i).varnodes.get(firstvisible).getTranscripts().get(1).getGenename(), headerlengths[5][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
							
						}
						else {
							if(bedarray.get(i).varnodes.get(firstvisible).getTranscripts().get(0).getEnd() > bedarray.get(i).varnodes.get(firstvisible).getPosition()) {
								
								buf.drawString(" ... " +bedarray.get(i).varnodes.get(firstvisible).getTranscripts().get(0).getGenename(), headerlengths[5][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
							}
							else {
								buf.drawString(bedarray.get(i).varnodes.get(firstvisible).getTranscripts().get(0).getGenename() +" ... ", headerlengths[5][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
								
							}
						}
					}	
				}
				
				buf.setColor(Color.darkGray);
				buf.drawLine(3, rowHeight+3, 3, (rowHeight*(i+genemutcount+2))-tablescroll.getVerticalScrollBar().getValue()+3);	
				
				for(int r = 0; r <headerlengths.length;r++) {					
					buf.drawLine(headerlengths[r][0], (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, headerlengths[r][0], (rowHeight*(i+genemutcount+2))-tablescroll.getVerticalScrollBar().getValue()+3);	
				}
				
				if(selectedNode != null && selectedNode.equals(bedarray.get(i))) {
					
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
									
									StringBuffer[] bedarraytemp = MethodLibrary.makeTrackArray(aminoarray.get(a).getNode(),aminoarray.get(a).getRow()[5], true);
									if(bedarraytemp != null) {
										int h = 0;
										for(int b = 0 ; b<bedarraytemp.length; b++) {
											if(b == bedtrack.trackIndex) {
												continue;
											}
											buf.setColor(Color.black);
											if(b == bedarraytemp.length-1) {
												buf.fillRect((int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+h)[1]+12, (rowHeight*(i+a+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth()-(int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+h)[1], rowHeight-1);	
											}
											else {
												buf.fillRect((int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+h)[1]+12, (rowHeight*(i+a+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+h)[2], rowHeight-1);												
											}
											buf.setColor(Color.white);
											if(bedarraytemp[b] != null) {
												buf.drawString(bedarraytemp[b].toString(), (int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+h)[1]+14, (rowHeight*(i+a+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
											
											}	
											h++;																	
												
										//	buf.drawLine((int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+h)[1], (rowHeight*(i+a+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(geneheaderlength+Control.controlData.fileArray.size()*2+h)[1], (rowHeight*(i+a+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+10);	
										
										}
									}									
								}
							}
							
							buf.setColor(Color.darkGray);
							for(int j = 0; j<geneheader.size(); j++) {								
								buf.drawLine((int)geneheader.get(j)[1]+11, (rowHeight*(i+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(j)[1]+11, (rowHeight*(i+s+listAdd+3))-tablescroll.getVerticalScrollBar().getValue()+3);	
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
		buf.drawLine(4, (rowHeight*(bedarray.size()+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(bedarray.size()+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3);	
	
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
	for(int i = 0; i<header.length; i++) {
		buf.setColor(Color.darkGray);
		buf.fillRect((int)(headerlengths[i][0]), 0, (int)(headerlengths[i][1])+1, rowHeight);
		
		if(mouseY-tablescroll.getVerticalScrollBar().getValue() <= rowHeight) {
			if(mouseX >= (int)(headerlengths[i][0]) && mouseX <= (int)(headerlengths[i][0])+(int)(headerlengths[i][1])) {
				headerHover = i;
				buf.setColor(Color.yellow);			
			}
			else if( mouseX >= (int)(headerlengths[headerlengths.length-1][0])+(int)(headerlengths[headerlengths.length-1][1])) {
				headerHover = headerlengths.length;
				buf.setColor(Color.white);			
			}
			else {
				buf.setColor(Color.white);		
			}			
			
			
		}
		else {
			
			buf.setColor(Color.white);
		}
		if(!mouseDrag && headerHover > -1 && i > 0 && i < header.length) {
			if(mouseX > (int)(headerlengths[i][0])-5 && mouseX < (int)(headerlengths[i][0])+5) {
				resizeColumn = i;
				
			}
		}
		buf.drawString(header[i], (int)(headerlengths[i][0])+4, rowHeight-2);
		buf.setColor(Color.black);
		buf.drawLine((int)(headerlengths[i][0]), 0, (int)(headerlengths[i][0]), rowHeight);
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
	buf.fillRect((int)headerlengths[header.length-1][0]+(int)headerlengths[header.length-1][1]+2, 0, this.width, rowHeight);
	buf.setColor(Color.black);
	buf.drawLine((int)headerlengths[header.length-1][0]+(int)headerlengths[header.length-1][1]+1,0, (int)headerlengths[header.length-1][0]+(int)headerlengths[header.length-1][1]+1, rowHeight);
	if(!mouseDrag && headerHover > -1 && resizeColumn == -1) {
		
		if(mouseX > (int)((int)headerlengths[header.length-1][0]+(int)headerlengths[header.length-1][1])-5 && mouseX < ((int)headerlengths[header.length-1][0]+(int)headerlengths[header.length-1][1])+5) {
			resizeColumn = header.length;
			
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
		else {
			BedTrack track =  (BedTrack)geneheader.get(i)[0];
			
			
			if(track.file != null) {
				buf.drawString(track.file.getName(), (int)geneheader.get(i)[1]+14, y+rowHeight-3);
			}
			else {
				buf.drawString(FilenameUtils.getName(track.url.getFile()), (int)geneheader.get(i)[1]+14, y+rowHeight-3);
			}
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



void drawScreen2(Graphics g) {	
	
	buf.setColor(Color.black);
	buf.fillRect(0, 0, this.getWidth(), this.getHeight());	
	
	//Header Draw	
	if(aminoarray.size() > 0) {
		samplecount=0;
		genemutcount = 0;
		for(int i = 0; i<aminoarray.size(); i++) {
			
			if((i+1)*rowHeight < tablescroll.getVerticalScrollBar().getValue()) {				
				continue;
			}
			if(i*rowHeight > tablescroll.getVerticalScrollBar().getValue() + this.getHeight()) {
				break;
			}
			
			/*if(mouseY >= (rowHeight*(i+samplecount+1)) && mouseY < (rowHeight*(i+samplecount+2))) {
				hoverNode = aminoarray.get(i).getNode();
				hoverString = aminoarray.get(i).getRow();
			}*/
			
			for(int j = 0; j<aminoarray.get(i).getRow().length; j++) {	
				try {
				
								
			//	buf.setColor(Color.gray);
			//	buf.drawLine(this.getWidth()/header.length*(j+1), rowHeight, this.getWidth()/header.length*(j+1), this.getHeight());
				
				
					
					if(selectedNode != null && selectedNode.equals(aminoarray.get(i).getNode())) {
					//	samplecount = selectedNode.getSamples().size();
			/*			buf.setColor(Color.black);
						
						buf.fillRect(this.getWidth()/header.length*j, ((rowHeight)*((i+1)+samplecount))+2-tablescroll.getVerticalScrollBar().getValue(), this.getWidth()/header.length, rowHeight);
						buf.setColor(Color.darkGray);
						buf.drawRect(this.getWidth()/header.length*j, ((rowHeight)*((i+1)+samplecount))+2-tablescroll.getVerticalScrollBar().getValue(), this.getWidth()/header.length, rowHeight);
						*/
						buf.setColor(Color.yellow);
						if(j < 2) {
							buf.drawString(aminoarray.get(i).getRow()[j], 10+this.getWidth()/header.length*j, (rowHeight*(i+1))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
						}
						else {
							buf.setColor(Color.yellow);							
							
							for(int s = 0; s<samplecount; s++) {								
								if(j == 2) {			
									buf.setColor(Color.yellow);
						//			buf.drawString(selectedNode.getSamples().get(s).getSample().getName(), 10, (rowHeight*(i+s+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
								}								
								buf.setColor(Color.black);
								buf.fillRect(this.getWidth()/header.length*j, ((rowHeight)*((i+s+2)))+2-tablescroll.getVerticalScrollBar().getValue(), this.getWidth()/header.length+1, rowHeight);
								buf.setColor(Color.yellow);
								buf.drawString(aminoarray.get(i).getRow()[j], 10+this.getWidth()/header.length*j, (rowHeight*((i+s+2)))+2-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
																
							}
						}
					/*	
						for(int s = 0; s<samplecount; s++) {
							
							buf.setColor(Color.black);
							buf.fillRect(this.getWidth()/header.length*j, ((rowHeight)*((i+1)+samplecount))+2-tablescroll.getVerticalScrollBar().getValue(), this.getWidth()/header.length+1, rowHeight);
							
							buf.setColor(Color.yellow);
							buf.drawString(aminoarray.get(i).getRow()[j], 10+this.getWidth()/header.length*j, (rowHeight*((i+1)+samplecount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
							
							buf.drawString(selectedNode.getSamples().get(s).getSample().getName(), 10, (rowHeight*(i+s+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
							
						}
				*/		
						
					}
					else {
						buf.setColor(Color.black);
						
						buf.fillRect(this.getWidth()/header.length*j, ((rowHeight)*((i+1)+samplecount))+2-tablescroll.getVerticalScrollBar().getValue(), this.getWidth()/header.length, rowHeight);
						buf.setColor(Color.darkGray);
						buf.drawRect(this.getWidth()/header.length*j, ((rowHeight)*((i+1)+samplecount))+2-tablescroll.getVerticalScrollBar().getValue(), this.getWidth()/header.length, rowHeight);
						if(hoverNode != null && hoverNode.equals(aminoarray.get(i).getNode())) {						
							buf.setColor(Color.yellow);
						}
						else {
							buf.setColor(Color.white);	
						}
						buf.drawString(aminoarray.get(i).getRow()[j], 10+this.getWidth()/header.length*j, (rowHeight*((i+1)+samplecount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
						
					}
					
				}
				catch(Exception e) {
					e.printStackTrace();
					ErrorLog.addError(e.getStackTrace());
				}
			}
				
		}
	}
	for(int i = 0; i<header.length; i++) {
		buf.setColor(Color.darkGray);
		buf.fillRect(i*this.getWidth()/header.length, 0, this.getWidth()/header.length+1, rowHeight);
		buf.setColor(Color.white);
		buf.drawString(header[i], i*this.getWidth()/header.length+4, rowHeight-2);
		buf.setColor(Color.black);
		buf.drawLine(i*this.getWidth()/header.length, 0, i*this.getWidth()/header.length, rowHeight);
	}
	g.drawImage(bufImage, 0, tablescroll.getVerticalScrollBar().getValue(), null);
	
}	
void addEntry(BedNode entry) {
	bedarray.add(entry);	
}
/*void addEntry(AminoEntry entry) {
	aminoarray.add(entry);	
}*/
void clear() {
		
		bedarray.clear();
		aminoarray.clear();		
	
		
		hoverNode = null;
		selectedNode = null;
		entry = null;
		vararray.clear();
	//	vardraw = null;
		hoverVar = null;
		selectedVar = null;
	
}
int getTableSize() {
	if(bedarray == null) {
		return 0;
	}
	return bedarray.size();
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


@Override
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
						
						while(searchHead.getPrev() != null) {
							if(searchHead.getPrev().getPosition() == 0) {
								searchHead.getPrev().putNext(searchHead);
							}
							searchHead = searchHead.getPrev();							
						}
						
						FileRead.head = searchHead;					
						searchHead = null;
						Main.drawCanvas.current = hoverNode.varnodes.get(0);
						
						Main.drawCanvas.gotoPos(hoverNode.getChrom().replace("chr",""), hoverNode.getPosition(),  hoverNode.getPosition()+hoverNode.getLength());						
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
										aminoarray = null;
										this.setPreferredSize(new Dimension(this.getWidth(), (this.getTableSize()+2)*rowHeight));
										
										this.revalidate();
										
										repaint();
										break;
									}
									
									selectedNode = hoverNode;
									selectedString = hoverString;
								//	samplecount = selectedNode.mutations;
									samplecount = selectedNode.varnodes.size();
									if(VariantHandler.tabs.getSelectedIndex() == 0) {
										
										this.setPreferredSize(new Dimension(this.getWidth(), (this.getTableSize()+2+samplecount)*rowHeight));
										this.revalidate();										
										getAminos(selectedNode);						
										this.repaint();									
									}
									else {										
										this.setPreferredSize(new Dimension(this.getWidth(), (this.getTableSize()+2+samplecount)*rowHeight));
										this.revalidate();										
										getAminos(selectedNode);						
										repaint();
									}									
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
/*
public String getRscode(VarNode varnode) {
	
	if(varnode.rscode > -1) {
	
		try {
			tabixreader = new TabixReader(Main.drawCanvas.sampleList.get(varnode.rscode).getTabixFile());
		
				if(varnode.getExons().size() > 0) {
					iterator = tabixreader.query(varnode.getExons().get(0).getTranscript().getChrom() +":" +(varnode.getPosition()+1)+"-"+(varnode.getPosition()+1));
				}
				else {
					iterator = tabixreader.query(varnode.getTranscripts().get(0).getChrom() +":" +(varnode.getPosition()+1)+"-"+(varnode.getPosition()+1));
					
				}
		//	}
			return iterator.next().split("\t")[2];			
			
		}
		catch(Exception e) {
			System.out.println(varnode.getPosition() +" " +varnode.rscode);
			e.printStackTrace();
			ErrorLog.addError(e.getStackTrace());
			return "-";
		}
	}
	return "-";
}
*/
@Override
public void mouseEntered(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}


@Override
public void mouseExited(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}

void createPolygon() {
	
	if(sorter.ascending) {
		int[] x = { (int)((headerlengths[sorter.index][0]+headerlengths[sorter.index][1]))-15, (int)((headerlengths[sorter.index][0]+headerlengths[sorter.index][1]))-10, (int)((headerlengths[sorter.index][0]+headerlengths[sorter.index][1]))-5};
		int[] y = { 12, 4, 12 };
		int n = 3;
		this.sortTriangle = new Polygon(x,y,n);
	}
	else {
		int[] x = { (int)((headerlengths[sorter.index][0]+headerlengths[sorter.index][1]))-15, (int)((headerlengths[sorter.index][0]+headerlengths[sorter.index][1]))-10, (int)((headerlengths[sorter.index][0]+headerlengths[sorter.index][1]))-5};
		int[] y = { 4, 12, 4 };
		int n = 3;
		this.sortTriangle = new Polygon(x,y,n);
	}
}
@Override


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
					Collections.sort(bedarray, sorter);
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

public static class ListSorter implements Comparator<BedNode> {
	public int index;
	boolean ascending = true;
	
	
	public int compare(BedNode o1, BedNode o2) {  
		
		
		if(index == 0) {
			String f1 = o1.name;
			String f2 = o2.name;
			 if ( f1.compareTo(f2) < 0) {  
		        	if(ascending) {
		        		return -1;  
		        	}
		        	else {
		        		return 1;
		        	}
		                
		        } else if ( f1.compareTo(f2) > 0 ) {  
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
		else if(index == 1) {
			int f1 = o1.mutations;
			int f2 = o2.mutations;
			
			 if ( f1 < f2) {  
		        	if(ascending) {
		        		return -1;  
		        	}
		        	else {
		        		return 1;
		        	}
		                
		        } else if ( f1 > f2) {  
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
		else if(index == 2) {
			 if ( Main.chromIndex.get(o1.getChrom())[0].longValue()+o1.getPosition() < Main.chromIndex.get(o2.getChrom())[0].longValue()+o2.getPosition()) {  
		        	if(ascending) {
		        		return -1;  
		        	}
		        	else {
		        		return 1;
		        	}
		                
		        } else if ( Main.chromIndex.get(o1.getChrom())[0].longValue()+o1.getPosition() > Main.chromIndex.get(o2.getChrom())[0].longValue()+o2.getPosition()) {  
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
		else if (index == 3) {
			int f1 = o1.getLength();
			int f2 = o2.getLength();
			
			 if ( f1 < f2) {  
		        	if(ascending) {
		        		return -1;  
		        	}
		        	else {
		        		return 1;
		        	}
		                
		        } else if ( f1 > f2) {  
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
		else if (index == 4) {
			double f1 = o1.mutations/(double)o1.getLength();
			double f2 = o2.mutations/(double)o2.getLength();
			
			 if ( f1 < f2) {  
		        	if(ascending) {
		        		return -1;  
		        	}
		        	else {
		        		return 1;
		        	}
		                
		        } else if ( f1 > f2) {  
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
void resizeTable(int column, int amount) {
	
	if(headerHover != -1) {
		if(headerlengths[column-1][1]+amount > 20) {
			headerlengths[column-1][1]+=amount;
			for(int i = column;i<headerlengths.length; i++) {
				headerlengths[i][0]+=amount;
			}
		}
		if(headerlengths[headerlengths.length-1][0]+headerlengths[headerlengths.length-1][1] > this.getWidth()) {
			if(bufImage.getWidth() < headerlengths[headerlengths.length-1][0]+headerlengths[headerlengths.length-1][1]) {
				bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
				buf = (Graphics2D)bufImage.getGraphics();
			}
			this.setPreferredSize(new Dimension(headerlengths[headerlengths.length-1][0]+headerlengths[headerlengths.length-1][1], this.getHeight()));
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
/*
public ArrayList<AminoEntry> getMutations(VarNode vardraw) {
	
	if(vardraw.getExons().size() > 0) {	
		
		for(int exon = 0; exon<vardraw.getExons().size(); exon++) {
			
			for(int i = 0; i<vardraw.getSamples().size(); i++) {
				if(Main.drawCanvas.hideVar(vardraw.getSamples().get(i))) {
					continue;
				}
				if(baseTemp.isEmpty() || !baseTemp.containsKey(vardraw.getSamples().get(i).getVariation())) {
					baseTemp.put(vardraw.getSamples().get(i).getVariation(), 1);
				}
				else {
					baseTemp.put(vardraw.getSamples().get(i).getVariation(), baseTemp.get(vardraw.getSamples().get(i).getVariation())+1);
				}
			}
			
			if(!baseTemp.isEmpty()) {
				
				e = baseTemp.keys();
			
				while(e.hasMoreElements()) {
					
					base = e.nextElement();
			//		String[] addrow = new String[5];						
			//		addrow[0] = vardraw.getExons().get(0).getTranscript().getGenename();						
			//		addrow[1] = ""+baseTemp.get(base);
			//		addrow[2] = Main.chromosomeDropdown.getSelectedItem() +":"+(vardraw.getPosition()+1);
			//		addrow[3] = Main.chromDraw.getAminoChange(vardraw,base,vardraw.getExons().get(exon));
			/*		if(vardraw.isRscode() > -1) {
						addrow[4] = "Yes";
					}
					else {
						addrow[4] = "No";
					}
				*/	
				//	AminoEntry entry = new AminoEntry(addrow, vardraw);
					
				//	VariantHandler.table.addEntry(entry);
				/*	if(vardraw.getExons().get(0).getTranscript().mutations == 0) {
						VariantHandler.table.addEntry(vardraw.getExons().get(0).getTranscript());
					}
					vardraw.getExons().get(0).getTranscript().mutations+=baseTemp.get(base);
					
					VariantHandler.table.repaint();
				//	aminoarray.add(new Entry(vardraw.getExons().get(0).getTranscript().getGenename() +" " +vardraw.getPosition() +" " +baseTemp.get(base) +" " +Main.chromDraw.getAminoChange(vardraw,base,vardraw.getExons().get(exon), true), vardraw));
											
				}
			}
			baseTemp.clear();
		}
	}	
	
	return null;
}
*/
@Override
public void mouseMoved(MouseEvent event) {
	if(!this.isEnabled()) {
		return;
	}
	mouseY = event.getY();
	mouseX = event.getX();
	
	repaint();
}

void getAminos(BedNode transcript) {
	
	try {
	aminoarray = null;
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
						if(!Main.drawCanvas.hideVar(entry.getValue().get(m),entry.getKey().length() > 1)) {
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
					base = entry.getKey();
					String[] addrow = new String[6];						
					addrow[0] = transcript.name;						
					addrow[1] = ""+mutcount;
					addrow[2] = transcript.getChrom() +":"+MethodLibrary.formatNumber((varnode.getPosition()+1));
					
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
					if(aminoarray == null) {
						aminoarray = new ArrayList<AminoEntry>();
					}
					aminoarray.add(aminoentry);							
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
