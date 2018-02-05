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


public class StatsTable extends JPanel implements MouseMotionListener, MouseListener, MouseWheelListener {

	
private static final long serialVersionUID = 1L;
	
		BufferedImage bufImage;
		ListSorter sorter = new ListSorter();
		int rowHeight = 15, geneheaderlength = 0;
		int variants = 0;
		Graphics2D buf;
		StringBuffer mutcountbuffer = new StringBuffer("");
		int width, height;
	//	ArrayList<Sample> aminoarray = new ArrayList<Sample>();
		ArrayList<Object[]> sampleArray = new ArrayList<Object[]>();
		ArrayList<Sample> samples = new ArrayList<Sample>();
		int mouseY=0, mouseX=0, pressX=0, pressY=0;
	//	VarNode hoverNode;
		final JScrollPane tablescroll;
		Object[] hoverNode, selectedNode;
		String[] posSplit, hoverString;
		int samplecount = 0;
		Enumeration<String> e;
		String base;
		
		ArrayList<String> geneheader =  new ArrayList<String>();
		String[] header =  {"Sample","Variants", "SNVs", "DELs", "INSs", "Coding", "Hetero/homo-rate", "TS/TV-rate", "T>A", "T>C", "T>G", "C>A", "C>G", "C>T", "Avg.call/cov"}; 
		
		int[][] headerlengths = new int[header.length][2];
		private int genemutcount;	
		VarNode hoverVar;
		
		int headerHover;
		boolean dot = false;		
		Polygon sortTriangle = new Polygon();
		
		ArrayList<SampleNode> vararray = new ArrayList<SampleNode>(), controlarray;
		
		VarNode varAdd;
		Map.Entry<String, ArrayList<SampleNode>>  entry;
	
		MethodLibrary.controlsorter ctrlsort = new MethodLibrary.controlsorter();
		private Color linecolor;
		private Sample sample;
		
		StatsTable(int width, int height, JScrollPane tablescroll) {			
			this.width = width;
			this.height = height;
			this.tablescroll = tablescroll;
			bufImage = MethodLibrary.toCompatibleImage(new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB));	
		//	bufImage = MethodLibrary.toCompatibleImage(new BufferedImage((int)width, (int)height, BufferedImage.TYPE_INT_ARGB));	
			buf = (Graphics2D)bufImage.getGraphics();
			this.addMouseListener(this);
			this.addMouseMotionListener(this);
		
			controlarray = new ArrayList<SampleNode>();
			headerlengths[0][0] = 0;
			headerlengths[0][1] = tablescroll.getViewport().getWidth()/header.length;
			for(int i = 1 ; i<headerlengths.length; i++) {
				headerlengths[i][0] = headerlengths[0][1]*i;
				headerlengths[i][1] = headerlengths[0][1];	
			}		
			
		}	
		
void resizeTable() {
	
	headerlengths[0][0] = 0;
	headerlengths[0][1] = tablescroll.getViewport().getWidth()/header.length;
	for(int i = 1 ; i<headerlengths.length; i++) {
		headerlengths[i][0] = headerlengths[0][1]*i;
		headerlengths[i][1] = headerlengths[0][1];	
	}
	
}

void drawScreen(Graphics g) {	

	buf.setColor(Color.black);
	
	buf.fillRect(0, 0, tablescroll.getViewport().getWidth(), tablescroll.getViewport().getHeight());	
	if(width != tablescroll.getViewport().getWidth()) {
		
		width = tablescroll.getViewport().getWidth();
		
		createPolygon();
		resizeTable();
	}


	genemutcount = 0;
	
	if(sampleArray.size() > 0) {		
		
		hoverVar = null;
		
		
		for(int i = 0; i<sampleArray.size(); i++) {
			dot = false;
						
			if((i+1+samplecount+sampleArray.size())*rowHeight < tablescroll.getVerticalScrollBar().getValue()) {
				
				continue;
			}
			
			if(i*rowHeight > tablescroll.getVerticalScrollBar().getValue() + tablescroll.getViewport().getHeight()) {
				break;
			}
			
			if(mouseY >= (rowHeight*(i+genemutcount+1)) && mouseY < (rowHeight*(i+genemutcount+2))) {
				hoverNode = sampleArray.get(i);	
			
			}			
			
			try {
				buf.setColor(Color.darkGray);
				buf.drawLine(4, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3, tablescroll.getViewport().getWidth(), (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3);	
				
				
				if(sampleArray.get(i).equals(hoverNode) || sampleArray.get(i).equals(selectedNode)) {
					linecolor = Color.yellow;
					
				}
				else {
					linecolor = Color.white;
					
				}
				buf.setColor(linecolor);
				
				for(int c=0; c<headerlengths.length; c++) {
					if(c == 0) {
						sample = (Sample)sampleArray.get(i)[c];
						buf.drawString((i+1) +".  " +sample.getName(), headerlengths[c][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
					}
					else {
						buf.drawString("" +sampleArray.get(i)[c], headerlengths[c][0]+5, (rowHeight*(i+1+genemutcount))-tablescroll.getVerticalScrollBar().getValue()+rowHeight);		
						
					}
					
					if(c < headerlengths.length-1) {
						buf.setColor(Color.black);
						buf.fillRect(headerlengths[c+1][0]+1, (rowHeight*(i+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+4, headerlengths[c][1], rowHeight-1);	
						buf.setColor(linecolor);				
					}
				}
				
			
			}
			catch(Exception e) {
				ErrorLog.addError(e.getStackTrace());
				e.printStackTrace();
			}	
			
		}
		buf.setColor(Color.darkGray);
		buf.drawLine(4, (rowHeight*(sampleArray.size()+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3, tablescroll.getViewport().getWidth(), (rowHeight*(sampleArray.size()+genemutcount+1))-tablescroll.getVerticalScrollBar().getValue()+3);	
	
	}
	headerHover = -1;
	for(int i = 0; i<header.length; i++) {
		buf.setColor(Color.darkGray);
		buf.fillRect(i*tablescroll.getViewport().getWidth()/header.length, 0, tablescroll.getViewport().getWidth()/header.length+1, rowHeight);
		
		if(mouseY-tablescroll.getVerticalScrollBar().getValue() < rowHeight && mouseX >i*tablescroll.getViewport().getWidth()/header.length && mouseX < (i+1)*tablescroll.getViewport().getWidth()/header.length) {
			headerHover = i;
			buf.setColor(Color.yellow);
		}
		else {
			
			buf.setColor(Color.white);
		}
		buf.drawString(header[i], i*tablescroll.getViewport().getWidth()/header.length+4, rowHeight-2);
		buf.setColor(Color.black);
		buf.drawLine(i*tablescroll.getViewport().getWidth()/header.length, 0, i*tablescroll.getViewport().getWidth()/header.length, rowHeight);
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
	
	
	
	if(headerHover == -1) {
		
		setCursor( Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));		
		
	}
	else {
		setCursor( Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));		
	}
	g.drawImage(bufImage, 0, tablescroll.getVerticalScrollBar().getValue(), null);
	
}	



void drawGeneheader(int y) {
	for(int i = 0; i<geneheader.size(); i++) {
		buf.setColor(Color.darkGray);
		buf.fillRect((int)(i*((tablescroll.getViewport().getWidth()-10)/(double)geneheader.size()))+10, y, (int)((tablescroll.getViewport().getWidth()-10)/(double)geneheader.size()+1), rowHeight);
		
		/*if(mouseY-tablescroll.getVerticalScrollBar().getValue() < rowHeight && mouseX >i*tablescroll.getViewport().getWidth()/geneheader.length && mouseX < (i+1)*tablescroll.getViewport().getWidth()/geneheader.length) {
			headerHover = i;
			buf.setColor(Color.yellow);
		}
		else {
			*/
			buf.setColor(Color.white);
//		}
		buf.drawString(geneheader.get(i), (int)(i*(tablescroll.getViewport().getWidth()-10)/(double)geneheader.size()+14), y+rowHeight-3);
		buf.setColor(Color.black);
		buf.drawLine((int)(i*(tablescroll.getViewport().getWidth()-10)/(double)geneheader.size())+10, y, (int)(i*(tablescroll.getViewport().getWidth()-10)/(double)geneheader.size())+10, y+rowHeight);
		
	}
}

void addEntry(String[] entry) {
	sampleArray.add(entry);	
}
/*void addEntry(AminoEntry entry) {
	aminoarray.add(entry);	
}*/
void clear() {
		samples.clear();
		sampleArray.clear();
		this.controlarray.clear();
		varAdd = null;
		hoverNode = null;
		selectedNode = null;
		entry = null;
		vararray.clear();
	//	vardraw = null;
		hoverVar = null;
		
	
}
int getTableSize() {
	return sampleArray.size();
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
		
			
		}
		case InputEvent.BUTTON3_MASK: {		
			
			selectedNode = null;
		
			repaint();
		}
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

void createPolygon() {
	if(sorter.ascending) {
		int[] x = { (sorter.index+1)*tablescroll.getViewport().getWidth()/header.length-15, (sorter.index+1)*tablescroll.getViewport().getWidth()/header.length-10, (sorter.index+1)*tablescroll.getViewport().getWidth()/header.length-5};
		int[] y = { 12, 4, 12 };
		int n = 3;
		this.sortTriangle = new Polygon(x,y,n);
	}
	else {
		int[] x = { (sorter.index+1)*tablescroll.getViewport().getWidth()/header.length-15, (sorter.index+1)*tablescroll.getViewport().getWidth()/header.length-10, (sorter.index+1)*tablescroll.getViewport().getWidth()/header.length-5};
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
			if(headerHover > -1) {
				
					if(sorter.ascending) {
						sorter.ascending = false;
					}
					else {
						sorter.ascending = true;
					}
					
					sorter.index = headerHover;
					Collections.sort(sampleArray, sorter);
					createPolygon();
					
					repaint();
				}
			
		}

		if(hoverNode != null || hoverVar != null) {
			Main.chromDraw.repaint();
		}
	}
	
}

public static class ListSorter implements Comparator<Object[]> {
	public int index;
	boolean ascending = true;
	
	
	public int compare(Object[] o1, Object[] o2) {  
		
		
		if(index == 0) {
			Sample s1 = (Sample)o1[0];
			Sample s2 = (Sample)o2[0];
			String f1 = s1.getName();
			String f2 = s2.getName();
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
		else {
						
			if(o1[index] instanceof Integer) {
				int f1 = (int)(o1[index]);
				int f2 = (int)(o2[index]);
				
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
				double f1 = (double)(o1[index]);
				double f2 = (double)(o2[index]);
				
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
			
			
		}
			
}  
}
@Override
public void mouseReleased(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}


@Override
public void mouseDragged(MouseEvent event) {
	
	// TODO Auto-generated method stub
	
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
