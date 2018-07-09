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
import java.awt.FontMetrics;
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
import java.util.Map.Entry;

import javax.swing.JPanel;
import javax.swing.JScrollPane;

import org.apache.commons.io.FilenameUtils;
public class AminoTable extends JPanel implements MouseMotionListener, MouseListener, MouseWheelListener {
	
private static final long serialVersionUID = 1L;
	
		BufferedImage bufImage;
		ListSorter sorter = new ListSorter();
		int rowHeight = 15, geneheaderlength = 0;
		Graphics2D buf;
		StringBuffer mutcountbuffer = new StringBuffer("");
		int width, height;
		ArrayList<AminoEntry> aminoarray = new ArrayList<AminoEntry>();
		//ArrayList<Transcript> transarray = new ArrayList<Transcript>();
		ArrayList<Gene> genearray = new ArrayList<Gene>();
		int mouseY=0, mouseX=0, pressX=0, pressY=0, dragX =0;
		FisherExact fe = new FisherExact(2000000);
	//	VarNode hoverNode;
		StringBuffer[] bedarray;
		final JScrollPane tablescroll;
		Gene hoverNode, selectedNode;
		String[] posSplit, hoverString;
		int samplecount = 0;
		Enumeration<String> e;
		String base;
		FontMetrics fm;
		ArrayList<Object[]> geneheader =  new ArrayList<Object[]>();
		String[] header = {"Gene", "Mut. count", "Gene position", "Description" };
		int[][] headerlengths = new int[header.length][2];		
		int variants = 0;
		private int genemutcount;	
		VarNode hoverVar;
		private VarNode selectedVar;
		private int listAdd;
		private String[] selectedString;
		private int pointer;
		int headerHover = -1;
		boolean dot = false;		
		Polygon sortTriangle = new Polygon();
		private int hoverSample = -1;
		private int mutcount;
		ArrayList<SampleNode> vararray = new ArrayList<SampleNode>();
		SampleNode[] controlarray;
		private Color textcolor;
		private double casefreq;
	
		private int textWidth = 0;
		MethodLibrary.controlsorter ctrlsort = new MethodLibrary.controlsorter();
		private int cases;
		private int resizeColumn;
		private boolean mouseDrag = false;
		private int geneHeaderHover = -1;
		private int firstrow;

		private SampleNode hoverSampleNode;

		private String hoverBase;

		AminoTable(int width, int height, JScrollPane tablescroll) {			
			this.width = width;
			this.height = height;
			this.tablescroll = tablescroll;
			Object[] obj = new Object[3]; obj[0] = "Sample"; 	obj[1] = 10; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
					 obj = new Object[3]; obj[0] = "Mut. count";obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
					 obj = new Object[3]; obj[0] = "Position"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
					 obj = new Object[3]; obj[0] = "Base change";obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);						
					 obj = new Object[3]; obj[0] = "Effect in isoforms"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
					 obj = new Object[3]; obj[0] = "Genotype"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
					 obj = new Object[3]; obj[0] = "Quality"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
					 obj = new Object[3]; obj[0] = "rs-code"; 	obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2]; obj[2] = (int)((width-10)/7.0); geneheader.add(obj);
		
			geneheaderlength = geneheader.size();
			
			bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width, (int)height, BufferedImage.TYPE_INT_ARGB));	
			buf = (Graphics2D)bufImage.getGraphics();
			this.addMouseListener(this);
			this.addMouseMotionListener(this);
			fm = buf.getFontMetrics();
	//		controlarray = new ArrayList<SampleNode>();
			
			headerlengths[0][0] = 0;
			headerlengths[0][1] = (int)(0.15*width);
			headerlengths[1][0] = headerlengths[0][1]+headerlengths[0][0];	
			headerlengths[1][1] = (int)(0.15*width);
			headerlengths[2][0] = headerlengths[1][1]+headerlengths[1][0];
			headerlengths[2][1] = (int)(0.2*width);
			headerlengths[3][0] = headerlengths[2][1]+headerlengths[2][0];
			headerlengths[3][1] = (int)(0.5*width);
			
			this.addMouseWheelListener(this);
				/*
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

void resizeTable() {	
	if(bufImage.getWidth() < headerlengths[headerlengths.length-1][0]+headerlengths[headerlengths.length-1][1]) {
		
		bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
		buf = (Graphics2D)bufImage.getGraphics();
		buf.setFont(Main.menuFont);
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
				buf.setFont(Main.menuFont);
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
				buf.setFont(Main.menuFont);
			}
			this.setPreferredSize(new Dimension((int)geneheader.get(geneheader.size()-1)[1]+(int)geneheader.get(geneheader.size()-1)[2], this.getHeight()));
			this.revalidate();
		}
	}
	createPolygon();
	repaint();
}

void drawScreen(Graphics g) {	
	if(!isEnabled()) {		
		return;
	}
	
	buf.setColor(Color.black);	
	buf.fillRect(0, 0, VariantHandler.tableScroll.getViewport().getWidth(), tablescroll.getViewport().getHeight());	
	
	
	if( VariantHandler.writetofile.isSelected()) {
		buf.setColor(Color.white);
		if(FileRead.output !=null && Main.drawCanvas.loading && Draw.variantcalculator) {
			buf.drawString("Writing results to " +FileRead.outputName, 10, 20);
		}
		else {
			buf.drawString("Press annotate to write results straight to file", 10, 20);
		}
		g.drawImage(bufImage, 0, tablescroll.getVerticalScrollBar().getValue(), null);
		return;
	}
	//Header Draw	
	
		genemutcount = 0;		
		hoverVar = null;
		hoverSample = -1;
		headerHover = -1;
		geneHeaderHover = -1;
		
		if(!mouseDrag) {
			resizeColumn = -1;
		}
		firstrow = tablescroll.getVerticalScrollBar().getValue()/rowHeight -samplecount-listAdd-aminoarray.size();
		
		if(firstrow < 0) {
			firstrow = 0;
		}
		for(int i = 0; i<genearray.size(); i++) {
			dot = false;			
			
			if((i+1+samplecount+aminoarray.size()+listAdd)*rowHeight < tablescroll.getVerticalScrollBar().getValue()) {
				
				continue;
			}
			
			if(i*rowHeight > tablescroll.getVerticalScrollBar().getValue() + tablescroll.getViewport().getHeight()) {
				break;
			}
			
			if(mouseY >= (rowHeight*(i+genemutcount+1)) && mouseY < (rowHeight*(i+genemutcount+2))) {
				hoverNode = genearray.get(i);					
			}			
			
			try {
				buf.setColor(Color.darkGray);
				buf.drawLine(4, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3);	
				
				
				if(genearray.get(i).equals(hoverNode) || genearray.get(i).equals(selectedNode)) {
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				textWidth = (int)fm.getStringBounds(""+(i+1) +".  " +genearray.get(i).getName(), buf).getWidth();
				if(genearray.get(i).intergenic) {
					if(genearray.get(i).varnodes.get(0).getTranscripts() == null) {
						buf.drawString((i+1) +".  " +genearray.get(i).getName(), 5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
					else if(genearray.get(i).varnodes.get(0).getTranscripts().size() == 2) {
						
						buf.drawString((i+1) +".  " +genearray.get(i).getName() +" ... " +genearray.get(i).varnodes.get(0).getTranscripts().get(1).getGenename(), 5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					}
					else if(genearray.get(i).varnodes.get(0).getPosition() < genearray.get(i).getStart()) {
						buf.drawString((i+1) +".  " +" ... " +genearray.get(i).getName(), 5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
					else {
						buf.drawString((i+1) +".  " +genearray.get(i).getName() +" ... " , 5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
				}
				else {
					buf.drawString((i+1) +".  " +genearray.get(i).getName(), 5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}
				buf.setColor(Color.black);
				buf.fillRect((int)(headerlengths[1][0]+1), (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, (int)(headerlengths[1][1]), rowHeight-1);	
				
				if(genearray.get(i).equals(hoverNode) || genearray.get(i).equals(selectedNode)) {
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				
				mutcountbuffer = new StringBuffer(""+genearray.get(i).mutations +" (");
				buf.drawString(mutcountbuffer.toString(), (int)(headerlengths[1][0]+5), (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
		
				if(genearray.get(i).nonsense > 0) {
					buf.setColor(Color.red);
					textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
					buf.drawString(""+genearray.get(i).nonsense, (int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					mutcountbuffer.append(genearray.get(i).nonsense);
					dot = true;
				}
				if(genearray.get(i).missense > 0) {
					
					if(dot) {
						buf.setColor(Color.white);
						textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
						buf.drawString(", ",(int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						mutcountbuffer.append(", ");
					}
					
					textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
					buf.setColor(Color.yellow);
					buf.drawString(""+genearray.get(i).missense, (int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					mutcountbuffer.append(genearray.get(i).missense);
					dot = true;
				}
				if(genearray.get(i).synonymous > 0) {
					
					if(dot) {
						buf.setColor(Color.white);
						textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
						buf.drawString(", ",(int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						mutcountbuffer.append(", ");
					}
					textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
					buf.setColor(Color.green);
					buf.drawString(""+genearray.get(i).synonymous, (int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					mutcountbuffer.append(genearray.get(i).synonymous);
					dot = true;
				}
				if(genearray.get(i).utr > 0) {
					
					if(dot) {
						buf.setColor(Color.white);
						textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
						buf.drawString(", ",(int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						mutcountbuffer.append(", ");
					}
					buf.setColor(Color.lightGray);
					textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
					buf.drawString(""+genearray.get(i).utr, (int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					mutcountbuffer.append(genearray.get(i).utr);
					dot = true;
				}
				if(genearray.get(i).intronic > 0) {
					
					if(dot) {
						buf.setColor(Color.white);
						textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
						buf.drawString(", ",(int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						mutcountbuffer.append(", ");
					}
					buf.setColor(Color.gray);
					textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
					buf.drawString(""+genearray.get(i).intronic, (int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					mutcountbuffer.append(genearray.get(i).intronic);
					dot = true;
				}
				if(genearray.get(i).intergenic) {
										
					buf.setColor(Color.gray);
					textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
					buf.drawString(""+genearray.get(i).mutations, (int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					mutcountbuffer.append(""+genearray.get(i).mutations);
				}
				buf.setColor(Color.white);
				textWidth = (int)fm.getStringBounds(mutcountbuffer.toString(), buf).getWidth();
				buf.drawString(") ", (int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				
				buf.setColor(Color.gray);
				textWidth = (int)fm.getStringBounds(mutcountbuffer.toString() +") ", buf).getWidth();
				if(genearray.get(i).samples.size() == 1) {
					buf.drawString(" 1 sample",  (int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					
				}
				else {
					buf.drawString(" " +genearray.get(i).samples.size() +" samples",  (int)(headerlengths[1][0])+5+textWidth, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}
				buf.setColor(Color.black);
				buf.fillRect((int)(headerlengths[2][0])+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
		
				if(genearray.get(i).equals(hoverNode) || genearray.get(i).equals(selectedNode)) {
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				if(genearray.get(i).intergenic) {
					if(genearray.get(i).varnodes.get(0).getTranscripts()== null) {
						buf.drawString(genearray.get(i).getChrom(), (int)(headerlengths[2][0])+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
					else if(genearray.get(i).varnodes.get(0).getTranscripts().size() == 2) {
						buf.drawString(genearray.get(i).getChrom() +":" +MethodLibrary.formatNumber(genearray.get(i).getEnd()) +"-" +MethodLibrary.formatNumber(genearray.get(i).varnodes.get(0).getTranscripts().get(1).getStart()), (int)(headerlengths[2][0])+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					}
					else if(genearray.get(i).varnodes.get(0).getPosition() < genearray.get(i).getStart()) {
						buf.drawString(genearray.get(i).getChrom() +":1-" +MethodLibrary.formatNumber(genearray.get(i).varnodes.get(0).getTranscripts().get(1).getStart()), (int)(headerlengths[2][0])+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
					else {
						buf.drawString(genearray.get(i).getChrom() +":" +MethodLibrary.formatNumber(genearray.get(i).getEnd()) +"-end", (int)(headerlengths[2][0])+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
				}
				else {
					buf.drawString(genearray.get(i).getChrom() +":"+MethodLibrary.formatNumber(genearray.get(i).getStart()) +"-" +MethodLibrary.formatNumber(genearray.get(i).getEnd()), (int)(headerlengths[2][0])+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}
				buf.setColor(Color.black);
				buf.fillRect((int)(headerlengths[3][0])+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
				
				if(genearray.get(i).equals(hoverNode) || genearray.get(i).equals(selectedNode)) {
					
					buf.setColor(Color.yellow);
				}
				else {
					buf.setColor(Color.white);
				}
				if(genearray.get(i).intergenic) {
					if (genearray.get(i).varnodes.get(0).getTranscripts() == null) {
						buf.drawString("-", (int)(headerlengths[3][0])+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
					else if(genearray.get(i).varnodes.get(0).getTranscripts().size() == 2) {
						buf.drawString(genearray.get(i).getDescription() +";" +genearray.get(i).varnodes.get(0).getTranscripts().get(1).getGene().getDescription(), (int)(headerlengths[3][0])+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
				}
				else {
					buf.drawString(genearray.get(i).getDescription(), (int)(headerlengths[3][0])+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
				}
				buf.setColor(Color.darkGray);
				buf.drawLine(3, rowHeight+3, 3, (rowHeight*(i+genemutcount+2))-tablescroll.getVerticalScrollBar().getValue()+3);	
				
				for(int r = 0; r <headerlengths.length;r++) {					
					buf.drawLine((int)(headerlengths[r][0]), (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, (int)(headerlengths[r][0]), (rowHeight*(i+genemutcount+2))-tablescroll.getVerticalScrollBar().getValue()+3);	
				}
				
				if(selectedNode != null && selectedNode.equals(genearray.get(i))) {
					
					hoverSample = -1;
					genemutcount = aminoarray.size()+1;						
					listAdd = 1;
			//		buf.drawLine(10, (rowHeight*(i+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(i+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+3);	
					drawGeneheader((rowHeight*(i+listAdd+1))-tablescroll.getVerticalScrollBar().getValue()+3);
					
					for(int s = 0; s<aminoarray.size(); s++) {
						
						buf.setColor(Color.darkGray);
						buf.drawLine(21, (rowHeight*(i+s+listAdd+3))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(i+s+listAdd+3))-tablescroll.getVerticalScrollBar().getValue()+3);	
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
						buf.setColor(textcolor);
						if(mouseY >= (rowHeight*(i+s+listAdd+2)) && mouseY < (rowHeight*(i+s+listAdd+3))) {
							hoverNode = null;
							hoverVar = aminoarray.get(s).getNode();		
							hoverString = aminoarray.get(s).getRow();
							buf.setColor(Color.white);						
							hoverSample = -1;
							
							if(aminoarray.get(s).getRow()[1].equals("1")) {
								for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
									if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(aminoarray.get(s).getRow()[5])) {
										
										hoverSample = aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getSample().getIndex();	
										hoverSampleNode = aminoarray.get(s).getNode().vars.get(v).getValue().get(0);
										hoverBase = aminoarray.get(s).getRow()[5];
										break;
									}
								}
							}
					//		hoverSample = -1;
										
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
								
							for(int h = 1; h< 4; h++) {
								buf.setColor(Color.black);
								buf.fillRect((int)geneheader.get(h)[1]+10, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(h)[2], rowHeight-1);	
								buf.setColor(textcolor);
								if(h == 3) {
									if(aminoarray.get(s).getRow()[5].length() == 1) {
										buf.drawString(Main.getBase.get(aminoarray.get(s).getNode().getRefBase()) +">" +aminoarray.get(s).getRow()[5],  (int)geneheader.get(h)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
									}
									else {
										buf.drawString(aminoarray.get(s).getRow()[5],  (int)geneheader.get(h)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
										
									}
									buf.setColor(Color.black);
									buf.fillRect((int)geneheader.get(4)[1]+10, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(4)[2], rowHeight-1);	
									buf.setColor(textcolor);
									buf.drawString(aminoarray.get(s).getRow()[h],  (int)geneheader.get(4)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
									
								}
								else {
									
									buf.drawString(aminoarray.get(s).getRow()[h],  (int)geneheader.get(h)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
									
								}
							}
							
							if(aminoarray.get(s).getRow()[1].equals("1")) {
								buf.setColor(Color.black);
								buf.fillRect((int)geneheader.get(5)[1]+10, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(5)[2], rowHeight-1);	
								buf.setColor(textcolor);
								
								for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
									if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(aminoarray.get(s).getRow()[5])) {
										if(aminoarray.get(s).getNode().vars.get(v).getValue().get(0).isHomozygous()) {
										buf.drawString("Hom (" +aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getCalls() +"/" +aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getCoverage() +")", (int)geneheader.get(5)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
										if(Control.controlData.controlsOn) {
											cases = 2;
											casefreq = 2/(double)(Main.varsamples*2-2);
										}
									}
									else {
										buf.drawString("Het ("+aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getCalls() +"/" +aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getCoverage() +")",  (int)geneheader.get(5)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
										if(Control.controlData.controlsOn) {
											cases = 1;
											casefreq = 1/(double)(Main.varsamples*2-1);
										}
									
									}
										buf.setColor(Color.black);
										buf.fillRect((int)geneheader.get(6)[1]+10, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
										buf.setColor(textcolor);
										buf.drawString(""+aminoarray.get(s).getNode().vars.get(v).getValue().get(0).getQuality(), (int)geneheader.get(6)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
									
									}
								}
							}
							else {
								//TODO piirra mustat boksit
								buf.setColor(Color.black);
								buf.fillRect((int)geneheader.get(5)[1]+10, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
								
								if(Control.controlData.controlsOn) {
									cases = 0;
									
									for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
										if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(aminoarray.get(s).getRow()[5])) {
											for(int j = 0; j<aminoarray.get(s).getNode().vars.get(v).getValue().size(); j++) {
												if(aminoarray.get(s).getNode().vars.get(v).getValue().get(j).alleles != null) {
													continue;
												}
												if(aminoarray.get(s).getNode().vars.get(v).getValue().get(j).isHomozygous()) {											
													cases += 2;
												}
												else {												
													cases += 1;																						
												}	
											}
										}																				
									}
									casefreq = cases/(double)(Main.varsamples*2-cases);			
									
								}								
							}
							buf.setColor(textcolor);
							buf.drawString(aminoarray.get(s).getRow()[4],  (int)geneheader.get(7)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
					//		buf.setColor(Color.black);
							
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
														buf.drawString(""+MethodLibrary.round(controlarray[e].alleles/(double)controlarray[e].allelenumber,5), (int)geneheader.get(this.geneheaderlength+e*2)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
														buf.setColor(Color.black);
														buf.fillRect((int)geneheader.get(this.geneheaderlength+e*2+1)[1]+11, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
														buf.setColor(textcolor);
														
														buf.drawString(""+MethodLibrary.round(casefreq/(controlarray[e].alleles/(double)(controlarray[e].allelenumber-controlarray[e].alleles)),2) +" (p=" +MethodLibrary.round(fe.getRightTailedP(cases, Main.varsamples*2-cases, controlarray[e].alleles, controlarray[e].allelenumber-controlarray[e].alleles) ,2) +")", (int)geneheader.get(this.geneheaderlength+e*2+1)[1]+14, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
												
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
							//if(Main.bedCanvas.bedOn) {					
								
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
								
								
							
							/*if(c < header.size()-1-Main.bedCanvas.bedTrack.size()) {
								buf.setColor(Color.black);
								buf.fillRect((int)header.get(c+1)[1]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, (int)header.get(c)[2], rowHeight-1);	
											
							}*/
						//	}
							
							buf.setColor(Color.darkGray);
							for(int j = 0; j<geneheader.size(); j++) {								
								
								buf.drawLine((int)geneheader.get(j)[1]+11, (rowHeight*(i+s+listAdd+2))-tablescroll.getVerticalScrollBar().getValue()+4, (int)geneheader.get(j)[1]+11, (rowHeight*(i+s+listAdd+3))-tablescroll.getVerticalScrollBar().getValue()+3);	
							}
						if(selectedVar!=null && selectedString.equals(aminoarray.get(s).getRow()) && Integer.parseInt(selectedString[1]) > 1) {
							pointer = 0;
							//TODO
							
							for(int v =0;v<aminoarray.get(s).getNode().vars.size(); v++) {
									if(aminoarray.get(s).getNode().vars.get(v).getKey().equals(selectedString[5])) {
										
										for(int l = 0; l<aminoarray.get(s).getNode().vars.get(v).getValue().size();l++) {
											if(aminoarray.get(s).getNode().vars.get(v).getValue().get(l).alleles != null) {
												break;
											}
											if(aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getSample().annotation) {
												continue;
											}
										if(mouseY > (rowHeight*(i+s+pointer+4)) && mouseY < (rowHeight*(i+s+pointer+5))) {
											textcolor = Color.white;
											
											hoverVar = aminoarray.get(s).getNode();		
											hoverString = aminoarray.get(s).getRow();
											hoverSample = aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getSample().getIndex();
											hoverSampleNode = aminoarray.get(s).getNode().vars.get(v).getValue().get(l);
											hoverBase = aminoarray.get(s).getRow()[5];
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
											buf.fillRect((int)geneheader.get(5)[1]+10, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
											buf.setColor(textcolor);
										if(aminoarray.get(s).getNode().vars.get(v).getValue().get(l).isHomozygous()) {
											buf.drawString("Hom (" +aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getCalls() +"/" +aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getCoverage() +")",  (int)geneheader.get(5)[1]+14, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
											
										}
										else {
											buf.drawString("Het ("+aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getCalls() +"/" +aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getCoverage() +")",  (int)geneheader.get(5)[1]+14, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
										}
											buf.setColor(Color.black);
											buf.fillRect((int)geneheader.get(6)[1]+10, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+4, this.getWidth(), rowHeight-1);	
											buf.setColor(textcolor);
											buf.drawString(""+aminoarray.get(s).getNode().vars.get(v).getValue().get(l).getQuality(),  (int)geneheader.get(6)[1]+14, (rowHeight*(i+s+pointer+3))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);				
											buf.setColor(Color.darkGray);
											for(int j = 5; j<7; j++) {								
												
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
		buf.drawLine(4, (rowHeight*(genearray.size()+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3, this.getWidth(), (rowHeight*(genearray.size()+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3);	
	
	
	
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

void addRowGeneheader(Object column) {
		Object[] obj = new Object[3]; 
		if(Main.bedCanvas.bedTrack.size() == 0 || column instanceof BedTrack) {			
			obj[0] = column; 	
			obj[1] = (int)geneheader.get(geneheader.size()-1)[1] + (int)geneheader.get(geneheader.size()-1)[2];
			obj[2] = 100;
			geneheader.add(obj);
		}
		else {
			
		}
		
		
		if((int)obj[1]+(int)obj[2] > this.getWidth()) {
			if(bufImage.getWidth() < (int)obj[1]+(int)obj[2]) {
				bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width*2, (int)height, BufferedImage.TYPE_INT_ARGB));	
				buf = (Graphics2D)bufImage.getGraphics();
			}
			this.setPreferredSize(new Dimension((int)obj[1]+(int)obj[2], this.getHeight()));
			this.revalidate();
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
			BedTrack track = (BedTrack)geneheader.get(i)[0];
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

void addEntry(Gene entry) {
	genearray.add(entry);	
}
/*void addEntry(AminoEntry entry) {
	aminoarray.add(entry);	
}*/
void clear() {
	
		genearray.clear();
		aminoarray.clear();	
		this.controlarray = null;
	
		hoverNode = null;
		selectedNode = null;
	//	entry = null;
	//	vararray.clear();
	//	vardraw = null;
		hoverVar = null;
		selectedVar = null;
}
int getTableSize() {
	if(genearray == null) {
		return 0;
	}
	return genearray.size();
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
					if(Draw.variantcalculator) {
						return;
					}
					FileRead.novars = true;
					
					if(hoverSample > -1) {
						Main.drawCanvas.drawVariables.visiblestart = (short)(hoverSample);
						Main.drawCanvas.drawVariables.visiblesamples = (short)(1);						
						Main.drawCanvas.resizeCanvas(Main.drawScroll.getViewport().getWidth(), (int)(Main.samples*Main.drawCanvas.drawVariables.sampleHeight));					
						Draw.setScrollbar((int)(hoverSample*Main.drawCanvas.drawVariables.sampleHeight));	
					}
					else {
						Main.drawCanvas.drawVariables.visiblestart = (short)(0);
						Main.drawCanvas.drawVariables.visiblesamples = (short)(Main.samples);						
						Main.drawCanvas.resizeCanvas(Main.drawScroll.getViewport().getWidth(), (int)(Main.samples*Main.drawCanvas.drawVariables.sampleHeight));					
						Draw.setScrollbar(0);	
					}
					
					if(hoverVar != null) {
						
						FileRead.search = true;
						VarNode searchHead = hoverVar;
					
						while(searchHead.getPrev() != null) {
							/*if(searchHead.getPrev().getPosition() == 0) {
								searchHead.getPrev().putNext(searchHead);
								
							}*/
							if(searchHead.getPrev().getChrom() == null || !searchHead.getPrev().getChrom().equals(hoverVar.getChrom())) {
								if(searchHead.getPrev().getNext() == null) {
									searchHead.getPrev().putNext(searchHead);
									
								}
								FileRead.head = searchHead.getPrev();
								
								break;
							}
							searchHead = searchHead.getPrev();							
						}					
						
						searchHead = null;
						
						Main.drawCanvas.current = hoverVar;	
						
						if(hoverVar.getExons() != null) {
							Main.drawCanvas.gotoPos(hoverVar.getChrom(), hoverVar.getPosition()+1-50, hoverVar.getPosition()+1+50);
						}
						else {
							if(hoverVar.getTranscripts().size() > 0) {
								Main.drawCanvas.gotoPos(hoverVar.getChrom(), hoverVar.getPosition()+1-50, hoverVar.getPosition()+1+50);
							}
						}
					}
					else if(hoverNode != null) {
						
						
						
						FileRead.search = true;					
						
						if(hoverNode.varnodes.get(0).getTranscripts() != null && hoverNode.varnodes.get(0).getTranscripts().size() == 2) {
							Main.drawCanvas.gotoPos(hoverNode.getChrom(), hoverNode.getEnd(), hoverNode.varnodes.get(0).getTranscripts().get(1).getGene().getStart());		
						}
						else {
							VarNode searchHead = hoverNode.varnodes.get(0);
							
							/*while(searchHead.getPrev() != null) {
								if(searchHead.getPrev().getPosition() == 0) {
									searchHead.getPrev().putNext(searchHead);
								}
								searchHead = searchHead.getPrev();							
							}*/
							
							while(searchHead.getPrev() != null) {
								/*if(searchHead.getPrev().getPosition() == 0) {
									searchHead.getPrev().putNext(searchHead);
									
								}*/
								if(searchHead.getPrev().getChrom() == null || !searchHead.getPrev().getChrom().equals(searchHead.getChrom())) {
									
									FileRead.head = searchHead.getPrev();
									break;
								}
								searchHead = searchHead.getPrev();							
							}					
							//FileRead.head = searchHead;					
							searchHead = null;
							Main.drawCanvas.current = hoverNode.varnodes.get(0);	
							Main.drawCanvas.gotoPos(hoverNode.getChrom(), hoverNode.getStart(), hoverNode.getEnd());		
						}
					}			
					Main.chromDraw.updateExons = true;
					Main.chromDraw.repaint();
					
					break;
				}
				else if(event.getClickCount() == 1 ) {
					
				/*	if(hoverVar != null && selectedNode.equals(selectedVar)) {
						selectedString = hoverString;
						
					}
					else*/
					
					//MethodLibrary.showVariantMenu(this, varOverLap, sampleOverLap, moveX+(int)selectedSplit.pixel, moveY);
					if(hoverVar != null && (selectedVar == null || !selectedVar.equals(hoverVar) )){
						
						selectedVar = hoverVar;
						selectedString = hoverString;
						
						if(selectedVar.isRscode() !=null) {
							hoverString[4] = selectedVar.rscode;
						}
					
						//MethodLibrary.showVariantMenu(this, hoverVar, null, mouseX+(int)Main.defaultFontSize*2, mouseY,hoverBase);
					/*	if(hoverBase != null) {
							MethodLibrary.showVariantMenu(this, hoverVar, null, 0, mouseY-Main.defaultFontSize*6,hoverBase);
						}
						*/
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
						else {
							if(hoverSampleNode != null) {
								MethodLibrary.showVariantMenu(this, hoverVar, hoverSampleNode, mouseX+(int)Main.defaultFontSize*2, mouseY,hoverBase);
							}
							
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
			selectedString = null;
			selectedNode = null;
			hoverVar = null;
			selectedVar = null;
			aminoarray.clear();
			this.setPreferredSize(new Dimension(this.getWidth(), (this.getTableSize()+2)*rowHeight));
			
			this.revalidate();
			
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
	if(sorter.index > headerlengths.length-1) {
		return;
	}
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
					Collections.sort(genearray, sorter);
					createPolygon();					
				}				
					repaint();
			}
			
		}

		if(hoverNode != null || hoverVar != null) {
			Main.chromDraw.repaint();
		}
	}
	
}

public static class ListSorter implements Comparator<Gene> {
	public int index;
	boolean ascending = true;
	
	
	public int compare(Gene o1, Gene o2) {  
		
		
		if(index == 0) {
			String f1 = o1.getName();
			String f2 = o2.getName();
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
		//	long f1,f2;

			
			//	f1 = Main.chromIndex.get(o1.getChrom())[0].longValue()+o1.getStart();
			//	f2 = Main.chromIndex.get(o2.getChrom())[0].longValue()+o2.getStart();
			 
			 if ( Main.chromIndex.get(o1.getChrom())[0].longValue()+o1.getStart() < Main.chromIndex.get(o2.getChrom())[0].longValue()+o2.getStart()) {  
		        	if(ascending) {
		        		return -1;  
		        	}
		        	else {
		        		return 1;
		        	}
		                
		        } else if ( Main.chromIndex.get(o1.getChrom())[0].longValue()+o1.getStart() > Main.chromIndex.get(o2.getChrom())[0].longValue()+o2.getStart()) {  
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
	// TODO Auto-generated method stub
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

void getAminos(Gene gene) {
	
	try {
	aminoarray.clear();
	
	VarNode varnode = null;
	Map.Entry<String, ArrayList<SampleNode>>  entry;
	
	for(int t = 0; t<gene.varnodes.size(); t++) {
		
		varnode = gene.varnodes.get(t);
		
		if(gene.intergenic) {
			
			for(int v = 0; v<varnode.vars.size(); v++) {
				entry = varnode.vars.get(v);
				mutcount = 0;
					for(int m = 0; m<entry.getValue().size(); m++) {
						if(entry.getValue().get(m).alleles != null) {							
							break;
						}
						if(entry.getValue().get(m).getSample().annotation) {
							entry.getValue().remove(m);
							m--;
							continue;
						}
						if(!Main.drawCanvas.hideVar(entry.getValue().get(m), entry.getKey().length() > 1) && !entry.getValue().get(m).getSample().annotation) {
							if(!VariantHandler.none.isSelected()) {
								if(!entry.getValue().get(m).inheritance) {
									if(!entry.getValue().get(m).getSample().equals(Main.drawCanvas.selectedSample)) {
										entry.getValue().remove(m);
										m--;
										continue;
									}
								}
								//mutcount++;
							}
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
					
				if(mutcount == 0) {
					continue;
				}
				
				base = entry.getKey();
				String[] addrow = new String[9];
				
				if(varnode.getTranscripts() == null) {
					addrow[0] = gene.getName();		
				}
				else if(varnode.getTranscripts().size() == 2) {
					addrow[0] = gene.getName() +" ... " +varnode.getTranscripts().get(1).getGenename();		
				}
				else if(varnode.getPosition() < gene.getStart()) {
					addrow[0] = "... " +varnode.getTranscripts().get(1).getGenename();	
				}
				else {
					addrow[0] = gene.getName() +" ...";	
				}
												
				addrow[1] = ""+mutcount;
				addrow[2] = gene.getChrom() +":"+MethodLibrary.formatNumber((varnode.getPosition()+1));
				addrow[3] = "Intergenic";
				if(varnode.isRscode() != null) {
					addrow[4] = varnode.rscode;
				}
				else {
					addrow[4] = "N/A";
				}
				addrow[5] = base;
				addrow[6] = "Intergenic";
				addrow[7] = "Intergenic";
				addrow[8] = "Intergenic";				
				AminoEntry aminoentry = new AminoEntry(addrow, varnode);				
				aminoarray.add(aminoentry);		
			}			
			continue;
		}
		if(varnode.getExons() != null) {
			
	//	for(int exon = 0; exon<varnode.getExons().size(); exon++) {
		/*	if(!VariantHandler.allIsoforms.isSelected() && varnode.getExons().get(exon).getTranscript().getGene().getCanonical() != null && !varnode.getExons().get(exon).getTranscript().isCanonical()) {
				continue;
			}*/
			/*if(!gene.equals(varnode.getExons().get(exon).getTranscript().getGene())) {
				continue;
			}
			*/
			if(!varnode.coding && !VariantHandler.utr.isSelected()) {	
			
				continue;
			}
			
			for(int v = 0; v<varnode.vars.size(); v++) {
				
					entry = varnode.vars.get(v);
					mutcount = 0;
					for(int m = 0; m<entry.getValue().size(); m++) {
						if(entry.getValue().get(m).alleles != null) {							
							continue;
						}
						if(entry.getValue().get(m).getSample().annotation) {
							
							
							continue;
						}
						if(!VariantHandler.none.isSelected()) {
							
							if(!entry.getValue().get(m).inheritance) {
								entry.getValue().remove(m);
								m--;
								continue;
							}
						}
						
						if(!Main.drawCanvas.hideVar(entry.getValue().get(m), entry.getKey().length() > 1) ) {
							
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
					
					if(mutcount == 0) {
						continue;
					}
					
					base = entry.getKey();
					String[] addrow = new String[9];						
					StringBuffer aminos = new StringBuffer(""), transcripts= new StringBuffer(""), exons =new StringBuffer(""), biotypes =new StringBuffer("");
					addrow[0] = gene.getName();						
					addrow[1] = ""+mutcount;
					addrow[2] = gene.getChrom() +":"+MethodLibrary.formatNumber((varnode.getPosition()+1));
					String aminochange;
					for(int exon = 0; exon<varnode.getExons().size(); exon++) {
						
						if(!varnode.getExons().get(exon).getTranscript().getGene().equals(gene)) {
							continue;
						}
						aminochange = Main.chromDraw.getAminoChange(varnode,base,varnode.getExons().get(exon));
						if(aminochange.contains("UTR") && !VariantHandler.utr.isSelected()) {				
							continue;
						}
						if(VariantHandler.nonsense.isSelected()) {
							if(!MethodLibrary.aminoEffect(aminochange).contains("nonsense")) {							
								continue;
							}					
						}
						else if(VariantHandler.synonymous.isSelected()) {
						
							if(MethodLibrary.aminoEffect(aminochange).contains("synonymous")) {							
								continue;
							}						
						}						
						
						if(aminos.length() == 0) {
							aminos.append(aminochange);
							transcripts.append(varnode.getExons().get(exon).getTranscript().getENST());
							biotypes.append(varnode.getExons().get(exon).getTranscript().getBiotype());
							exons.append(varnode.getExons().get(exon).getNro());
						}
						else {
							aminos.append(";" +aminochange);
							transcripts.append(";" +varnode.getExons().get(exon).getTranscript().getENST());
							biotypes.append(";" +varnode.getExons().get(exon).getTranscript().getBiotype());
							exons.append(";" +varnode.getExons().get(exon).getNro());
						}
					}
					addrow[3]=aminos.toString();
					
					
					if(varnode.isRscode() != null) {
						addrow[4] = varnode.rscode;
					}
					else {
						addrow[4] = "N/A";
					}
					addrow[5] = base;
					addrow[6] = transcripts.toString();
					addrow[7] = biotypes.toString();
					addrow[8] = exons.toString();
			//		varAdd.putNext(null);
			//		varAdd.putPrev(null);
					
					AminoEntry aminoentry = new AminoEntry(addrow, varnode);
					/*if(varnode.getPosition() == 127796784) {
						System.out.println("Jou");
					}*/
					aminoarray.add(aminoentry);					
			
			}
			
	//	}
	}
		if(VariantHandler.intronic.isSelected() && varnode.isInGene() &&  varnode.getTranscripts() != null && varnode.getExons() == null) {
			
								
				for(int v = 0; v<varnode.vars.size(); v++) {
						entry = varnode.vars.get(v);
						base = entry.getKey();
						mutcount = 0;
						for(int m = 0; m<entry.getValue().size(); m++) {
							if(entry.getValue().get(m).alleles != null) {
								
								break;
							}
							if(entry.getValue().get(m).getSample().annotation) {
								entry.getValue().remove(m);
								m--;
								continue;
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
						if(mutcount == 0) {
							continue;
						}
						StringBuffer transcripts= new StringBuffer(""), biotypes =new StringBuffer("");
						
						String[] addrow = new String[9];						
						addrow[0] = gene.getName();			
						addrow[1] = ""+mutcount;
						addrow[2] = gene.getChrom() +":"+MethodLibrary.formatNumber((varnode.getPosition()+1));					
						addrow[3] = Main.getBase.get(varnode.getRefBase()) +"->" +base +" (intronic)";							
						
						if(varnode.isRscode() != null) {
							addrow[4] = varnode.rscode;
						}
						else {
							addrow[4] = "N/A";
						}
						addrow[5] = base;	
						for(int trans = 0; trans<varnode.getTranscripts().size(); trans++) {
							if(!varnode.getTranscripts().get(trans).getGene().equals(gene)) {
								continue;
							}
								
							if(transcripts.length() == 0) {
								
								transcripts.append(varnode.getTranscripts().get(trans).getENST());
								biotypes.append(varnode.getTranscripts().get(trans).getBiotype());
								
							}
							else {
								
								transcripts.append(";" +varnode.getTranscripts().get(trans).getENST());
								biotypes.append(";" +varnode.getTranscripts().get(trans).getBiotype());
								
							}
					}	
					
					
					if(varnode.isRscode() != null) {
						addrow[4] = varnode.rscode;
					}
					else {
						addrow[4] = "N/A";
					}
					addrow[5] = base;
					addrow[6] = transcripts.toString();
					addrow[7] = biotypes.toString();
					addrow[8] = "Intronic";
						AminoEntry aminoentry = new AminoEntry(addrow, varnode);						
						aminoarray.add(aminoentry);						
												
					}
			}
		
		}
	
	
		varnode = null;
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
		}
	
}
/*boolean checkCompounds(Gene gene) {
	
	boolean mother = false;
	boolean father = false;
	
	Entry<String, ArrayList<SampleNode>> entry;
	for(int i = 0 ;i<gene.varnodes.size(); i++) {
		for(int v = 0; v<gene.varnodes.get(i).vars.size(); v++) {
			
			entry = gene.varnodes.get(i).vars.get(v);
			for(int m = 0; m<entry.getValue().size(); m++) {
				if(!entry.getValue().get(m).inheritance) {
					continue;
				}
				if(entry.getValue().size() != 2) continue;
				if(entry.getValue().get(m).getSample().children != null) {
			 		if(entry.getValue().get(m).getSample().female) {
			 			mother = true;
			 		}
			 		else {
			 			father = true;
			 		}
			 	}
			}			
		}
	}
	if(!mother || !father) {
		for(int i = 0 ;i<gene.varnodes.size(); i++) {
			for(int v = 0; v<gene.varnodes.get(i).vars.size(); v++) {
				
				entry = gene.varnodes.get(i).vars.get(v);
				for(int m = 0; m<entry.getValue().size(); m++) {
					if(!entry.getValue().get(m).inheritance) {
						continue;
					}
					if(entry.getValue().size() != 2) continue;
					entry.getValue().get(m).inheritance = false;
				}
			}
		}
		return false;
	}
	return true;
}*/

@Override
public void mouseWheelMoved(MouseWheelEvent e) {
	if ( e.getWheelRotation() < 0 ) {	
		tablescroll.getVerticalScrollBar().setValue(tablescroll.getVerticalScrollBar().getValue()-16);
		
	}
	else {			
		tablescroll.getVerticalScrollBar().setValue(tablescroll.getVerticalScrollBar().getValue()+16);
		
	}
	
}
}
