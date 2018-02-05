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
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.Serializable;

import javax.swing.JCheckBox;
import javax.swing.JPopupMenu;

public class ControlFile implements Serializable, ActionListener {
	
	private static final long serialVersionUID = 1L;
	public boolean controlled = false;
	String tabixfile;
	String sampleName;
	short index;
	double alleleFreq = 0.01;
	int varcount = 0;
	Rectangle alleleBox = new Rectangle(), playbox = new Rectangle();
	JCheckBox remOverlaps;
	private transient JPopupMenu menu;
	Polygon playTriangle = new Polygon();
	
	boolean controlOn = false;
	StringBuffer alleletext = new StringBuffer("0.01");
	
	
	public ControlFile(String sampleName, short index, String tabixfile) {
		this.tabixfile = tabixfile;
		this.sampleName = sampleName;
		this.index = index;		
		setMenu();
	}

	void setMenu() {
		if(remOverlaps == null) {
			remOverlaps = new JCheckBox("Overlap indels");
		}
		
		menu = new JPopupMenu("Options");
		menu.add(remOverlaps);
		remOverlaps.addActionListener(this);
		for(int c = 0 ; c<getPopupMenu().getComponentCount(); c++) {
			getPopupMenu().getComponent(c).setFont(Main.menuFont);
		}
	}
	public JPopupMenu getPopupMenu() {
		return this.menu;
	}
	String getTabixFile() {
		return this.tabixfile;
	}	
	short getIndex() {
		return this.index;
	}
	void setIndex(short index) {
		this.index = index;
	}	
	String getName() {
		return this.sampleName;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if(e.getSource() == remOverlaps) {
			
			Control.dismissControl(FileRead.head, this);
		}		
	}
}
