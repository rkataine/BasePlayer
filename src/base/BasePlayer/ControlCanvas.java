
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
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.KeyEventDispatcher;
import java.awt.KeyboardFocusManager;
import java.awt.Rectangle;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.Graphics2D;
import java.util.ArrayList;

import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.bed.BEDCodec;

import javax.swing.JPanel;

public class ControlCanvas extends JPanel implements MouseMotionListener, MouseListener, KeyListener {

private static final long serialVersionUID = 1L;
	
	BufferedImage bufImage, nodeImage;
	Graphics2D buf, nodebuf;
	int width, height, bedwidth;
	Color zoomColor = new Color(115,115,100,100);	
	Composite backupComposite;
	Rectangle remoBox = new Rectangle();
	TabixReader tabixReader = null;
	TabixReader.Iterator iterator = null;
	TabixReader.Iterator bedIterator = null;	
	Rectangle mouseRect = new Rectangle();
	BEDCodec bedcodec = null;
	private int trackstart;
	private int selectedBox = -1;
	private boolean overlapping = false;
	private boolean typing = false;
	private int cursorPosition = 0;
	private int typeTextWidth = 0;	
	FontMetrics fm;
	private int typeBox = 0;
	private int selectedPlay = 0;
	private boolean sidebar = false;	
	private int hoverIndex = -1;	
	private int removeControl = -1;
	private Rectangle sideMouseRect = new Rectangle();
	private boolean resize;
	private int pressY;
	ArrayList<Double> trackDivider = new ArrayList<Double>();
	private boolean mouseDrag;
	private int resizeDivider = 0;
	private boolean resizer;
	private int preresize = 0;
	private double preresizer = 0.0;
	private int mouseY = 0, mouseX = 0;
	boolean zoomDrag = false;	
	private ArrayList<Double> tempDivider;
	boolean lineZoomer;

	private boolean positivelock;

	private boolean negative;

	private boolean negativelock;

	private String zoomtext;

	private Rectangle2D textWidth;

	private int zoompostemp;

	private int trackheight;	
	
	ControlCanvas(int width, int height) {	
	
		this.width = width;
		this.height = height;
		bufImage = new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB);	
		buf = (Graphics2D)bufImage.getGraphics();
		buf.setRenderingHints(Draw.rh);
		nodeImage = new BufferedImage((int)Main.screenSize.getWidth(), (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB);
		nodebuf = (Graphics2D)nodeImage.getGraphics();
		nodebuf.setRenderingHints(Draw.rh);
		backupComposite = nodebuf.getComposite();
		buf.setFont(Draw.defaultFont);
		nodebuf.setFont(Draw.defaultFont);
		fm = buf.getFontMetrics();
		addMouseMotionListener(this);
		addMouseListener(this);
		addKeyListener(this);
	}
	

void drawScreen(Graphics g) {	
	
	buf.drawImage(Draw.image, Main.sidebarWidth, 0, this.getWidth(),(int)Main.screenSize.getHeight(), this);
	buf.setColor(Draw.backColor);
	buf.fillRect(Main.sidebarWidth, 0, this.getWidth(), nodeImage.getHeight());	
	buf.setColor(Color.gray);
	
	if(this.trackDivider.get(this.trackDivider.size()-1) != 1.0) {	
		 for(int i = 0 ; i<this.trackDivider.size(); i++) {
			 this.trackDivider.set(i, ((i+1)*(this.getHeight()/(double)trackDivider.size())/this.getHeight()));
		 }				
	}

	drawSidebar();
	buf.setStroke(Draw.strongStroke);
	buf.setColor(Color.black);
	buf.drawLine(Main.sidebarWidth, 0, Main.sidebarWidth, this.getHeight());
	if(resizer && !mouseDrag) {
		resizer = false;
	}
	if(Control.controlData.fileArray.size() > 1) {
		buf.setStroke(Draw.doubleStroke);
	/*	for(int i = 1; i< Control.controlData.fileArray.size(); i++) {
			buf.drawLine(0, (int)(i*(Main.controlScroll.getViewport().getHeight()/(double)Control.controlData.fileArray.size()))-2, Main.controlScroll.getViewport().getWidth(), (int)(i*(Main.controlScroll.getViewport().getHeight()/(double)Control.controlData.fileArray.size()))-2);
		}*/
		for(int i = 0 ; i<Control.controlData.fileArray.size(); i++) {
			
			if(i <Control.controlData.fileArray.size()-1) {
				buf.drawLine(0, (int)(trackDivider.get(i)*this.getHeight()), this.getWidth(), (int)(trackDivider.get(i)*this.getHeight()));
				
				if(!lineZoomer && mouseY < (int)(trackDivider.get(i)*this.getHeight())+4 && mouseY > (int)(trackDivider.get(i)*this.getHeight())-4) {
					resizer = true;
					
					if(getCursor().getType() != Cursor.N_RESIZE_CURSOR) {
						resizeDivider = i;
						setCursor(Cursor.getPredefinedCursor(Cursor.N_RESIZE_CURSOR));					
					}	
				}			
			}		
		}
	}
	
	if(!resizer && !overlapping) {		
		if(getCursor().getType() != Cursor.DEFAULT_CURSOR) {			
			setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
		}
	}
	drawZoom();
	g.drawImage(bufImage, 0, 0, null);	
}	

void drawSidebar() {
	
	buf.setColor(Draw.sidecolor);
	buf.fillRect(0, 0, Main.sidebarWidth, this.getHeight());
	buf.setColor(Draw.softColor);
	buf.fillRect(0, 0, Main.sidebarWidth, this.getHeight());
	buf.setColor(Color.black);
	buf.setStroke(Draw.doubleStroke);
	if(Control.controlData.fileArray.size() > 0) {
		
		overlapping = false;
		for(int i = 0; i<Control.controlData.fileArray.size(); i++) {
			if(i == 0) {
				trackstart = 10;
				trackheight = (int)(trackDivider.get(i)*this.getHeight());
			}
			else {
				
				trackstart = 10+(int)(trackDivider.get(i-1)*this.getHeight());
				trackheight = (int)(trackDivider.get(i)*this.getHeight())-(int)(trackDivider.get(i-1)*this.getHeight());
			}
			
			buf.drawString((i+1) +": " +Control.controlData.fileArray.get(i).getName(), 10, trackstart);
			if(trackheight > 30) {
				buf.drawString("Allele count: "+Control.controlData.fileArray.get(i).varcount, 10, trackstart+15);
			}
			if(trackheight > 60) {
				if((int)Control.controlData.fileArray.get(i).alleleBox.getBounds().getY() != trackstart+30) {
					Control.controlData.fileArray.get(i).alleleBox.setBounds(10, trackstart+30, 100, 20);
					Control.controlData.fileArray.get(i).playbox.setBounds((int)Control.controlData.fileArray.get(i).alleleBox.getMaxX()+10, trackstart+30, 20, 20);
					Control.controlData.fileArray.get(i).playTriangle.reset();
					
					Control.controlData.fileArray.get(i).playTriangle.addPoint(Control.controlData.fileArray.get(i).playbox.x+4, trackstart+34);
					Control.controlData.fileArray.get(i).playTriangle.addPoint(Control.controlData.fileArray.get(i).playbox.x+4, trackstart+46);
					Control.controlData.fileArray.get(i).playTriangle.addPoint(Control.controlData.fileArray.get(i).playbox.x+16, trackstart+40);
				}
			//	if( Control.controlData.fileArray.get(i).varcount > 2) {
					buf.setColor(Color.white);
					buf.fillRect(Control.controlData.fileArray.get(i).alleleBox.getBounds().x,Control.controlData.fileArray.get(i).alleleBox.getBounds().y,Control.controlData.fileArray.get(i).alleleBox.getBounds().width, Control.controlData.fileArray.get(i).alleleBox.getBounds().height );
					buf.setColor(Color.black);
					buf.drawString(Control.controlData.fileArray.get(i).alleletext.toString(), 12, trackstart+45);
			//	}
				
				buf.setColor(Color.lightGray);
				buf.fillRect(Control.controlData.fileArray.get(i).playbox.getBounds().x,Control.controlData.fileArray.get(i).playbox.getBounds().y,Control.controlData.fileArray.get(i).playbox.getBounds().width, Control.controlData.fileArray.get(i).playbox.getBounds().height );
				
				if(Control.controlData.fileArray.get(i).controlOn) {
					buf.setColor(Color.green);
					buf.fillRect(Control.controlData.fileArray.get(i).playTriangle.getBounds().x, Control.controlData.fileArray.get(i).playTriangle.getBounds().y, 12,12);
				}
				else {
					buf.setColor(Color.red);
					buf.fillPolygon(Control.controlData.fileArray.get(i).playTriangle);
				}				
				if(sideMouseRect.intersects(Control.controlData.fileArray.get(i).playbox)) {
					overlapping = true;
					if(getCursor().getType() != Cursor.HAND_CURSOR) {
						selectedPlay = i;
						setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
					}	
					buf.setColor(Color.white);
					buf.drawRect(Control.controlData.fileArray.get(i).playbox.getBounds().x-1,Control.controlData.fileArray.get(i).playbox.getBounds().y-1,Control.controlData.fileArray.get(i).playbox.getBounds().width+1, Control.controlData.fileArray.get(i).playbox.getBounds().height+1);
					
				}
				buf.setColor(Color.black);
				
				if(sideMouseRect.intersects(Control.controlData.fileArray.get(i).alleleBox)/* &&  Control.controlData.fileArray.get(i).varcount > 2*/) {
					overlapping = true;
					
					if(getCursor().getType() != Cursor.TEXT_CURSOR) {
						selectedBox = i;
						setCursor(Cursor.getPredefinedCursor(Cursor.TEXT_CURSOR));
					}		
					
				}
			}
		}
		if(typing && Main.controlScroll.getViewport().getHeight()/Control.controlData.fileArray.size() > 55) {
			if(typeBox > -1) {
				buf.drawLine(Control.controlData.fileArray.get(typeBox).alleleBox.getBounds().x+this.typeTextWidth+2, Control.controlData.fileArray.get(typeBox).alleleBox.getBounds().y+2, Control.controlData.fileArray.get(typeBox).alleleBox.getBounds().x+this.typeTextWidth+2, Control.controlData.fileArray.get(typeBox).alleleBox.getBounds().y+20);
			}
		}
		if(!overlapping) {			
				
			if(getCursor().getType() != Cursor.DEFAULT_CURSOR) {
				selectedBox = -1;
				selectedPlay = -1;
				setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
			}			
		}
	}
	if(sidebar ) {
		buf.setColor(Color.black);
		buf.setStroke(Draw.basicStroke);
		
		if(hoverIndex > -1 && (this.remoBox.getBounds().y != (int)(trackDivider.get(hoverIndex)*this.getHeight())-11 || this.remoBox.getBounds().x != Main.sidebarWidth-11)) {
			this.remoBox.setBounds(Main.sidebarWidth-11, (int)(trackDivider.get(hoverIndex)*this.getHeight())-11, 8, 8);
		}
	
	if(sidebar && sideMouseRect.intersects(this.remoBox)) {				
		removeControl = hoverIndex;
		buf.setColor(Color.white);
		buf.fillRect(this.remoBox.x-2, this.remoBox.y-1, this.remoBox.width+2, this.remoBox.height+2);
	}
	else {				
		removeControl = -1;
	}
	
	buf.setColor(Color.black);
	buf.drawRect(this.remoBox.x-1, this.remoBox.y, this.remoBox.width, this.remoBox.height);
	buf.drawLine(this.remoBox.x-1,this.remoBox.y,this.remoBox.x+7, this.remoBox.y +8);
	buf.drawLine(this.remoBox.x-1,this.remoBox.y+8,this.remoBox.x+7, this.remoBox.y);
	
	}
	buf.setColor(Color.black);
	
}

void drawNodes() {
	
	//nodebuf.setColor(Draw.backColor);
	nodebuf.setComposite( Main.drawCanvas.composite);		
	nodebuf.fillRect(0,0, bufImage.getWidth(),this.getHeight());	
	nodebuf.setComposite(this.backupComposite);	
	

}

void drawZoom() {	
	
	if(lineZoomer ) {
		
		buf.setColor(Color.black);
		buf.setStroke(Draw.dashed);
		
		buf.drawLine(Main.drawCanvas.pressX, pressY, mouseX, mouseY);
	
		
	}
	else if(zoomDrag) {
		
		
		buf.setStroke(Draw.dashed);
		buf.setColor(Color.black);
		buf.setFont(ChromDraw.seqFont);
		if(this.mouseX-Main.drawCanvas.pressX >= 0) {
			buf.drawRect(Main.drawCanvas.pressX, 0, this.mouseX-Main.drawCanvas.pressX, this.getHeight());	
			if(Main.drawCanvas.getDrawWidth()-this.mouseX > 200) {
				buf.drawString("" +MethodLibrary.formatNumber((int)((this.mouseX-Main.drawCanvas.pressX)/Main.drawCanvas.selectedSplit.pixel)) +"bp", Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX)+4, this.mouseY-35);
				buf.drawString("Right click to cancel zoom" , Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX)+4, this.mouseY-6);
			}
			else {
				fm = buf.getFontMetrics();
				zoomtext = ""+MethodLibrary.formatNumber((int)((this.mouseX-Main.drawCanvas.pressX)/Main.drawCanvas.selectedSplit.pixel)) +"bp";
				textWidth = fm.getStringBounds(zoomtext, buf);
				this.zoompostemp = (int)((Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX))-textWidth.getWidth());
				buf.drawString(zoomtext, this.zoompostemp, this.mouseY-35);
				
				textWidth = fm.getStringBounds("Right click to cancel zoom", buf);
				this.zoompostemp = (int)((Main.drawCanvas.pressX+(this.mouseX-Main.drawCanvas.pressX))-textWidth.getWidth())-6;
				buf.drawString("Right click to cancel zoom" , this.zoompostemp, this.mouseY-6);
	
			}
			buf.setColor(Draw.zoomColor);
			buf.fillRect(Main.drawCanvas.pressX, 0, this.mouseX-Main.drawCanvas.pressX, this.getHeight());
		
		}
		else {
			
			lineZoomer = true;
			Main.drawCanvas.lineZoomer = true;
			zoomDrag = false;
		
		
		}
		buf.setFont(Draw.defaultFont);
	}
	buf.setStroke(Draw.doubleStroke);
}

public void paint(Graphics g) {
	try {	
		
		drawScreen(g);
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

public void mouseClicked(MouseEvent event) {
	switch(event.getModifiers()) {	
	case InputEvent.BUTTON1_MASK: {	
		if(this.selectedBox > -1) {
			this.requestFocus();
			typing  = true;
			typeBox = this.selectedBox;
			if(event.getClickCount() == 2) {
				Control.controlData.fileArray.get(selectedBox).alleletext = new StringBuffer("");
				typeTextWidth = 0;
				cursorPosition = 0;
			}
			else {
				typeTextWidth = (int)fm.getStringBounds(Control.controlData.fileArray.get(typeBox).alleletext.toString(), buf).getWidth();
				if(event.getX() > 12+typeTextWidth) {
					cursorPosition = Control.controlData.fileArray.get(selectedBox).alleletext.length();
					
				}
				else {
					int pos = 0;
					while((int)fm.getStringBounds(Control.controlData.fileArray.get(typeBox).alleletext.toString().substring(0,pos), buf).getWidth()+12 <= event.getX()) {
						pos++;
					}				
					cursorPosition = pos;
					typeTextWidth = (int)fm.getStringBounds(Control.controlData.fileArray.get(typeBox).alleletext.substring(0, cursorPosition), buf).getWidth();
				}			
			}	
		}
		else {
			typing = false;
		}
		if(this.selectedPlay > -1 && Control.controlData.fileArray.get(selectedPlay).playbox.intersects(sideMouseRect)) {
			
			if(!Control.controlData.fileArray.get(selectedPlay).controlOn) {
				Control.controlData.fileArray.get(selectedPlay).alleleFreq = Double.parseDouble(Control.controlData.fileArray.get(selectedPlay).alleletext.toString());
				Control.controlData.fileArray.get(selectedPlay).controlOn = true;
				Control.controlData.controlsOn = true;
			
				Control.runner runner = new Control.runner();
				runner.execute();
			}
			else {
				Control.controlData.fileArray.get(selectedPlay).controlOn = false;
				boolean allfalse = true;
				for(int i = 0; i<Control.controlData.fileArray.size(); i++) {
					if(Control.controlData.fileArray.get(i).controlOn) {
						allfalse = false;
						break;
					}
				}
				if(allfalse) {
					Control.controlData.controlsOn = false;
				}
				if(VariantHandler.commonSlider.getValue() > 1 && VariantHandler.clusterSize > 0) {
			    	Main.drawCanvas.calcClusters(FileRead.head);
				}
			//	Control.useCheck.setSelected(false);
				Draw.updatevars = true;
				Main.drawCanvas.repaint();
			}
			
		}
		if(removeControl > -1) {			
			removeControl(removeControl);
		}
		repaint();	
		break;
	}
	case InputEvent.BUTTON3_MASK: {			
		Control.controlData.fileArray.get(hoverIndex).menu.show(this, mouseX, mouseY);
		//this.bedTrack.get(hoverIndex).getPopup().show(this, mouseX, mouseY);
	}	
	}	
}	

public void removeControl(int removeControl) {
	MethodLibrary.removeHeaderColumns(Control.controlData.fileArray.get(removeControl));
	Control.controlData.fileArray.remove(removeControl);
	hoverIndex = -1;
	trackDivider.remove(removeControl);
	
	if(Control.controlData.fileArray.size() == 0) {
		
		Main.controlScroll.setVisible(false);
		 Main.trackPane.setDividerLocation(0);		
		 Main.trackPane.setDividerSize(0);
		if(!Main.bedScroll.isVisible()) {
			Main.varpane.setDividerSize(0);
			Main.trackPane.setVisible(false);
			Main.trackPane.setDividerSize(0);
			Main.varpane.setResizeWeight(0.0);   
		}
	}
	boolean allfalse = true;
	for(int i = 0; i<Control.controlData.fileArray.size(); i++) {
		if(Control.controlData.fileArray.get(i).controlOn) {
			allfalse = false;
			break;
		}
	}
	if(allfalse) {
		Control.controlData.controlsOn = false;
	}
}
public void mouseEntered(MouseEvent arg0) {}	
public void mouseExited(MouseEvent arg0) {}	
public void mousePressed(MouseEvent event) {
	resize =false;
	pressY = event.getY();
	mouseDrag = true;
	
	this.requestFocus();
	Main.drawCanvas.pressX = event.getX();	
	Main.drawCanvas.tempDrag = Main.drawCanvas.pressX;
	
	switch(event.getModifiers()) {	
		case InputEvent.BUTTON1_MASK: {	
			if(resizer) {
				preresize = mouseY;
				preresizer = this.trackDivider.get(this.resizeDivider)*this.getHeight();
				tempDivider = (ArrayList<Double>)trackDivider.clone();			
			}		
			break;
		}
		case InputEvent.BUTTON3_MASK: {			
			this.zoomDrag = false;
		}
}
	
	
	
}	
public void mouseReleased(MouseEvent event) {
	if(!sidebar) {
		if(zoomDrag) {			
			if(mouseX-Main.drawCanvas.pressX > 0) {			
				Main.drawCanvas.gotoPos(Main.drawCanvas.selectedSplit.start+(Main.drawCanvas.pressX-Main.drawCanvas.selectedSplit.offset)/Main.drawCanvas.selectedSplit.pixel, Main.drawCanvas.selectedSplit.start+(mouseX-Main.drawCanvas.selectedSplit.offset)/Main.drawCanvas.selectedSplit.pixel);
			}
		}
	
		Main.chromDraw.updateExons = true;
		repaint();
		Main.drawCanvas.mouseDrag = false;
		Main.chromDraw.repaint();
		Draw.updatevars = true;
		Main.drawCanvas.repaint();
	}
	mouseDrag = false; 
	zoomDrag = false;	
	Main.drawCanvas.lineZoomer = false;
	lineZoomer = false;
}	
public void mouseDragged(MouseEvent event) {
	
switch(event.getModifiers()) {	
	
	case InputEvent.BUTTON1_MASK: {	
		/*if(sidebar) {
			return;
		}*/
		this.mouseX = event.getX();
		this.mouseY = event.getY();		
		
		if(resizer && this.getHeight() > this.trackDivider.size()*20) {				
				
				this.trackDivider.set(this.resizeDivider, mouseY/(double)this.getHeight());
			
				if((positivelock || negative) && this.trackDivider.get(this.resizeDivider)*this.getHeight()-preresizer >= 0) {	
					positivelock = false;
					negative = false;
					preresize = mouseY;					
					tempDivider = (ArrayList<Double>)trackDivider.clone();
					
				}
				else if ((!negative || negativelock ) && this.trackDivider.get(this.resizeDivider)*this.getHeight()-preresizer < 0) {					
					negativelock = false;
					preresize = mouseY;					
					tempDivider = (ArrayList<Double>)trackDivider.clone();
					negative = true;
				}
				
				if(negativelock) { // || (this.trackDivider.get(this.resizeDivider+1)*this.getHeight()) -(this.trackDivider.get(this.resizeDivider)*this.getHeight()) < 20) {
					this.trackDivider.set(this.resizeDivider, (this.trackDivider.get(this.resizeDivider+1)*this.getHeight()-19)/(double)this.getHeight());
				}
				if(positivelock) { // || (this.trackDivider.get(this.resizeDivider+1)*this.getHeight()) -(this.trackDivider.get(this.resizeDivider)*this.getHeight()) < 20) {
					this.trackDivider.set(this.resizeDivider, 19/(double)this.getHeight());
				}
				if(this.trackDivider.get(this.resizeDivider)*this.getHeight()-preresizer < 0) {
					
					negative = true;
					positivelock = true;
					if(this.resizeDivider > 0) {
						
						for(int i =1 ; i<this.resizeDivider+1; i++) {								
							if((this.trackDivider.get(i)*this.getHeight())-(this.trackDivider.get(i-1)*this.getHeight()) < 20) {
								
								this.trackDivider.set(i,((this.trackDivider.get(i-1)*this.getHeight())+19)/(double)this.getHeight());	
							}
							else {
								positivelock = false;
								if(i != this.resizeDivider) {
									this.trackDivider.set(i,(this.tempDivider.get(i)/preresize)*mouseY);			
								}
							}
						}						
					}
					if((this.trackDivider.get(0)*this.getHeight()) >= 20) {
						this.trackDivider.set(0,(this.tempDivider.get(0)/preresize)*mouseY);			
						positivelock = false;
					}
					else {
						this.trackDivider.set(0, 19/(double)this.getHeight());
					}
					
				}
				else {
					
					negative = false;
					negativelock = true;
					if(this.resizeDivider < this.trackDivider.size()-1) {						
						
						for(int i =this.resizeDivider ; i<this.trackDivider.size()-1; i++) {
							
							if((this.trackDivider.get(i+1)*this.getHeight())-(this.trackDivider.get(i)*this.getHeight()) < 20) {
								
								this.trackDivider.set(i,((this.trackDivider.get(i+1)*this.getHeight())-19)/(double)this.getHeight());
							}
							else {
								
								negativelock = false;
								if(i !=this.resizeDivider) {
									try {
										this.trackDivider.set(i,1-((1-this.tempDivider.get(i))/(this.getHeight()-preresize))*(this.getHeight()-mouseY));
									}
									catch(Exception e) {
									//	System.out.println(this.tempDivider.get(i));
										e.printStackTrace();
									}
								}
							}							
						}		
						if(this.getHeight()-(this.trackDivider.get(this.trackDivider.size()-2)*this.getHeight()) >= 20) {
							negativelock = false;
						}
						else {
							this.trackDivider.set(this.trackDivider.size()-2, (this.getHeight()-19)/(double)this.getHeight());
						}
						
					}
				}				
				
				preresizer = this.trackDivider.get(this.resizeDivider)*this.getHeight();
				repaint();
				
			
		}
		else if(lineZoomer) {	
			
			 if (Main.drawCanvas.selectedSplit.start > 1 || Main.drawCanvas.selectedSplit.end < Main.drawCanvas.selectedSplit.chromEnd ) {					
				 Main.drawCanvas.gotoPos(Main.drawCanvas.selectedSplit.start-(Main.drawCanvas.tempDrag-mouseX)/Main.drawCanvas.selectedSplit.pixel*2, Main.drawCanvas.selectedSplit.end+(Main.drawCanvas.tempDrag-mouseX)/Main.drawCanvas.selectedSplit.pixel*2);
			 }
			
			Main.drawCanvas.tempDrag = mouseX;
			Main.chromDraw.updateExons = true;
			repaint();
			Main.chromDraw.repaint();
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
			
		}
		else {							
			zoomDrag = true;				
			
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
public void mouseMoved(MouseEvent event) {
	mouseY = event.getY();
	mouseX = event.getX();
	this.mouseRect.setBounds(event.getX()-Main.sidebarWidth, event.getY(), 1, 1);	
	this.sideMouseRect.setBounds(event.getX(), event.getY(), 1, 1);
	
	if(!sidebar && event.getX() < Main.sidebarWidth) {
		sidebar = true;
	}
	else {
		if(sidebar && event.getX() >= Main.sidebarWidth) {
			sidebar = false;
		}
	}
	for(int i = 0; i<this.trackDivider.size(); i++) {
		if(this.trackDivider.get(i)*this.getHeight() < mouseY) {
			continue;
		}
		if(this.trackDivider.get(i)*this.getHeight() >= mouseY) {			
			if(hoverIndex != i) {
				hoverIndex = i;				
			}
			break;
		}	
	}
	repaint();
}

@Override
public void keyTyped(KeyEvent e) {
	// TODO Auto-generated method stub
	
}

@Override
public void keyPressed(KeyEvent e) {
	if(typing) {
		
		if(e.getKeyCode() == KeyEvent.VK_RIGHT) {
			if(cursorPosition < Control.controlData.fileArray.get(typeBox).alleletext.length()) {
				cursorPosition++;
			}		
		}
		else if(e.getKeyCode() == KeyEvent.VK_LEFT) {
			if(cursorPosition > 0) {
				cursorPosition--;
			}
		}
		else if(e.getKeyCode() == KeyEvent.VK_BACK_SPACE) {
			if(cursorPosition > 0) {
				Control.controlData.fileArray.get(typeBox).alleletext.deleteCharAt(cursorPosition-1);
				cursorPosition--;
			}
		}
		else if(e.getKeyCode() == KeyEvent.VK_DELETE) {
			if(cursorPosition < Control.controlData.fileArray.get(typeBox).alleletext.length()) {
				Control.controlData.fileArray.get(typeBox).alleletext.deleteCharAt(cursorPosition);
			}
		}
		else if(e.getKeyCode() == KeyEvent.VK_ENTER) {
			try {
				Control.controlData.fileArray.get(typeBox).alleleFreq = Double.parseDouble(Control.controlData.fileArray.get(typeBox).alleletext.toString());
			
					Control.controlData.fileArray.get(typeBox).controlOn = true;
					Control.controlData.controlsOn = true;
					Control.runner runner = new Control.runner();
					runner.execute();
		
			}
			catch(Exception ex) {
				
			}
		}
		else {
			Control.controlData.fileArray.get(typeBox).alleletext.insert(cursorPosition, e.getKeyChar());
			cursorPosition++;
		}			
		
		typeTextWidth = (int)fm.getStringBounds(Control.controlData.fileArray.get(typeBox).alleletext.substring(0, cursorPosition), buf).getWidth();
		repaint();
	
}
	
}


@Override
public void keyReleased(KeyEvent e) {
	// TODO Auto-generated method stub
	
}
	
	
}
