
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
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.awt.Graphics2D;

import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.bed.BEDCodec;

import javax.swing.JPanel;

public class ControlCanvas extends JPanel implements MouseMotionListener, MouseListener, KeyListener {

private static final long serialVersionUID = 1L;
	
	BufferedImage bufImage, nodeImage;
	Graphics2D buf, nodebuf;
	int width, height, bedwidth;
//	private TabixReader tabixReader;
	Color zoomColor = new Color(115,115,100,100);
	//ArrayList<BedNode> bednodes = new ArrayList<BedNode>();
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
	

	drawSidebar();
	buf.setStroke(Draw.strongStroke);
	buf.setColor(Color.black);
	buf.drawLine(Main.sidebarWidth, 0, Main.sidebarWidth, this.getHeight());
	if(Control.controlData.fileArray.size() > 1) {
		buf.setStroke(Draw.doubleStroke);
		for(int i = 1; i< Control.controlData.fileArray.size(); i++) {
			buf.drawLine(0, (int)(i*(Main.controlScroll.getViewport().getHeight()/(double)Control.controlData.fileArray.size()))-2, Main.controlScroll.getViewport().getWidth(), (int)(i*(Main.controlScroll.getViewport().getHeight()/(double)Control.controlData.fileArray.size()))-2);
		}
	}
	
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
			trackstart = 10+(i*Main.controlScroll.getViewport().getHeight()/Control.controlData.fileArray.size());
			buf.drawString((i+1) +": " +Control.controlData.fileArray.get(i).getName(), 10, trackstart);
			buf.drawString("Allele count: "+Control.controlData.fileArray.get(i).varcount, 10, trackstart+15);
			if(Main.controlScroll.getViewport().getHeight()/Control.controlData.fileArray.size() > 55) {
				if((int)Control.controlData.fileArray.get(i).alleleBox.getBounds().getY() != trackstart+30) {
					Control.controlData.fileArray.get(i).alleleBox.setBounds(10, trackstart+30, 100, 20);
					Control.controlData.fileArray.get(i).playbox.setBounds((int)Control.controlData.fileArray.get(i).alleleBox.getMaxX()+10, trackstart+30, 20, 20);
					Control.controlData.fileArray.get(i).playTriangle.reset();
					
					Control.controlData.fileArray.get(i).playTriangle.addPoint(Control.controlData.fileArray.get(i).playbox.x+4, trackstart+34);
					Control.controlData.fileArray.get(i).playTriangle.addPoint(Control.controlData.fileArray.get(i).playbox.x+4, trackstart+46);
					Control.controlData.fileArray.get(i).playTriangle.addPoint(Control.controlData.fileArray.get(i).playbox.x+16, trackstart+40);
				}
				if( Control.controlData.fileArray.get(i).varcount > 2) {
					buf.setColor(Color.white);
					buf.fillRect(Control.controlData.fileArray.get(i).alleleBox.getBounds().x,Control.controlData.fileArray.get(i).alleleBox.getBounds().y,Control.controlData.fileArray.get(i).alleleBox.getBounds().width, Control.controlData.fileArray.get(i).alleleBox.getBounds().height );
					buf.setColor(Color.black);
					buf.drawString(Control.controlData.fileArray.get(i).alleletext.toString(), 12, trackstart+45);
				}
				
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
				if(mouseRect.intersects(Control.controlData.fileArray.get(i).playbox)) {
					overlapping = true;
					if(getCursor().getType() != Cursor.HAND_CURSOR) {
						selectedPlay = i;
						setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
					}	
					buf.setColor(Color.white);
					buf.drawRect(Control.controlData.fileArray.get(i).playbox.getBounds().x-1,Control.controlData.fileArray.get(i).playbox.getBounds().y-1,Control.controlData.fileArray.get(i).playbox.getBounds().width+1, Control.controlData.fileArray.get(i).playbox.getBounds().height+1);
					
				}
				buf.setColor(Color.black);
				
				if(mouseRect.intersects(Control.controlData.fileArray.get(i).alleleBox) &&  Control.controlData.fileArray.get(i).varcount > 2) {
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
		
		if(this.remoBox.getBounds().y != ((hoverIndex+1)*Main.controlScroll.getViewport().getHeight()/Control.controlData.fileArray.size())-11) {
			this.remoBox.setBounds(Main.sidebarWidth-11, ((hoverIndex+1)*Main.controlScroll.getViewport().getHeight()/Control.controlData.fileArray.size())-11, 8, 8);
		}
	
	if(sidebar && mouseRect.intersects(this.remoBox)) {				
	//	buf.setStroke(Draw.doubleStroke);
		removeControl  = hoverIndex;
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



public void paint(Graphics g) {
	try {	
		
		drawScreen(g);
	}
	catch(Exception e) {
		e.printStackTrace();
	}
}

public void mouseClicked(MouseEvent event) {
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
	if(this.selectedPlay > -1 && Control.controlData.fileArray.get(selectedPlay).playbox.intersects(mouseRect)) {
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
		//	Control.useCheck.setSelected(false);
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
		}
		
	}
	if(removeControl > -1) {
	
		//VariantHandler.table.geneheader.remove(index)
		VariantHandler.table.geneheader.remove(VariantHandler.table.geneheaderlength +removeControl*2);
		VariantHandler.table.geneheader.remove(VariantHandler.table.geneheaderlength +removeControl*2);
		VariantHandler.table.repaint();
		Control.controlData.fileArray.remove(removeControl);
		
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
	repaint();	
}	
public void mouseEntered(MouseEvent arg0) {}	
public void mouseExited(MouseEvent arg0) {}	
public void mousePressed(MouseEvent arg0) {}	
public void mouseReleased(MouseEvent arg0) {}	
public void mouseDragged(MouseEvent arg0) {}	
public void mouseMoved(MouseEvent event) {
	this.mouseRect.setBounds(event.getX(), event.getY(), 1, 1);
	if(!sidebar && event.getX() < Main.sidebarWidth) {
		sidebar = true;
	}
	else {
		if(sidebar && event.getX() >= Main.sidebarWidth) {
			sidebar = false;
		}
	}
	if(hoverIndex != (int)(event.getY()/((double)Main.controlScroll.getViewport().getHeight()/Control.controlData.fileArray.size())) ) {
		hoverIndex = (int)(event.getY()/((double)Main.controlScroll.getViewport().getHeight()/Control.controlData.fileArray.size()));
		if(hoverIndex > Control.controlData.fileArray.size()-1) {
			hoverIndex =  Control.controlData.fileArray.size()-1;
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
