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
import java.awt.Component;
import java.awt.Dimension;



import java.awt.Font;
import java.awt.GridLayout;
import java.awt.MouseInfo;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.Vector;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.ListCellRenderer;

public class SampleDialog extends JPanel implements ActionListener, KeyListener {
	
	private static final long serialVersionUID = 1L;
	static Boolean loading = false;
	static JFrame frame = new JFrame("Sample data");
	static JLabel samplelabel = new JLabel("Sample name:");
	static JLabel motherlabel = new JLabel("Mother:");
	static JLabel fatherlabel = new JLabel("Father:");
	static JLabel sexlabel = new JLabel("Gender:");
	static JTextField sampleNameField = new JTextField("");
	static JTextField vcfpath = new JTextField("");
	static JTextField bampath = new JTextField("");
	static ComboBoxRenderer comborenderer = new ComboBoxRenderer();
	static Sample sample;
	static Vector<Sample> motherlist = new Vector<Sample>();
	static Vector<Sample> fatherlist = new Vector<Sample>();
	static Color[] colors={new Color(224,160,100),new Color(238,80,68),new Color(126,126,230),new Color(104,192,106)};
	static SampleComboBox motherdrop, fatherdrop;	
	static JComboBox<Color> colorBox;
	static JComboBox<String> sexdrop;
	static boolean fromUI = false;
	static JPanel panel = new JPanel(new GridLayout(0,1));
	static JCheckBox affected = new JCheckBox("Affected");
	static JCheckBox annotation = new JCheckBox("Annotation");
	static JCheckBox intersect = new JCheckBox("Intersect");
	private static int annotationpointer;
	ActionListener parentDropActionListener = new ActionListener() {
		public void actionPerformed(ActionEvent actionEvent) {   	
    	
    	  if (actionEvent.getActionCommand() == "comboBoxChanged") {
    		 if(fromUI) {
    			 return;
    		 }
    		 if(actionEvent.getSource() == fatherdrop) {
				 SampleComboBox box = (SampleComboBox)actionEvent.getSource();
				 if(sample.father != null) {
					 sample.father.children.remove(sample);					 
				 }
				 sample.father = (Sample)box.getSelectedItem();
				 if(sample.father != null) {
					 if(sample.father.children == null) {
						 sample.father.children = new ArrayList<Sample>();
					 }
					 sample.parents++;
					 sample.father.female = false;
					 sample.father.children.add(sample);
					 sample.father.familyColor = (Color)colorBox.getSelectedItem();
					 for(int i = 0 ; i<sample.father.children.size(); i++) {
						 sample.father.children.get(i).familyColor = (Color)colorBox.getSelectedItem();
					 }
				 }
				 else {
					 sample.parents--;
				 }
    		 }
    		 else {
    			 SampleComboBox box = (SampleComboBox)actionEvent.getSource();
    			 if(sample.mother != null) {
					 sample.mother.children.remove(sample);					 
				 }
    			 sample.mother = (Sample)box.getSelectedItem();    	
    			 if(sample.mother != null) {
    				 sample.parents++;
					 if(sample.mother.children == null) {
						 sample.mother.children = new ArrayList<Sample>();
					 }
					 sample.mother.female = true;
					 sample.mother.children.add(sample);
					 sample.mother.familyColor = (Color)colorBox.getSelectedItem();
					 for(int i = 0 ; i<sample.mother.children.size(); i++) {
						 sample.mother.children.get(i).familyColor = (Color)colorBox.getSelectedItem();
					 }
				 }
    			 else {
    				 sample.parents--;
    			 }
    		 }
    		 setPanels();
    		 checkFiles();
    		 Main.drawCanvas.repaint();
    		
		}
		}
	};
	ActionListener colorAction = new ActionListener() {
		public void actionPerformed(ActionEvent actionEvent) {    	
    	  if (actionEvent.getActionCommand() == "comboBoxChanged") {  
    		  if(fromUI) {
    			 return; 
    		  }
    		  JComboBox<Color> color = (JComboBox)actionEvent.getSource();
    		  if(sample.father != null) {
    			  sample.father.familyColor = (Color)color.getSelectedItem();
    			  for(int i = 0 ; i<sample.father.children.size(); i++) {
    				  sample.father.children.get(i).familyColor = (Color)color.getSelectedItem();
        		  }  
    		  }
    		  if(sample.mother != null) {
    			  sample.mother.familyColor = (Color)color.getSelectedItem();
    			  for(int i = 0 ; i<sample.mother.children.size(); i++) {
    				  sample.mother.children.get(i).familyColor = (Color)color.getSelectedItem();
        		  }  
    		  } 		  
    	  }
    	  Main.drawCanvas.repaint();
		}
	};
	
	
	static void checkFiles() {
		int affected = 0;
		for(int i = 0 ; i<Main.drawCanvas.sampleList.size(); i++) {
			Sample sample = Main.drawCanvas.sampleList.get(i);
			if(sample.affected) {
				affected++;
			}
			
		  /*if(sample.children != null && sample.children.size() > 0) {
			  continue;
		  }*/
		  if(sample.father != null || sample.mother != null) {
			  if(!Main.drawCanvas.sampleList.contains(sample.father)) {
				  sample.father = null;
			  }
			  if(!Main.drawCanvas.sampleList.contains(sample.mother)) {
				  sample.mother = null;
			  }
			  if(sample.mother != null && sample.father != null) {
				  sample.parents = 2;
			  }
			  else {
				  sample.parents = 1;
			  }
			  
		  }
		  else {
			  sample.parents = 0;
			  if(sample.children == null || sample.children.size() == 0) {
				  sample.familyColor = null;	
			  }
			 		
		  }
		 
		}
		FileRead.affected = affected;
	}
	
	void setDropboxes() {
		
		colorBox = new JComboBox(colors);
		colorBox.setMaximumRowCount(5);
		colorBox.setPreferredSize(new Dimension(50,20));
		colorBox.setRenderer(new MyCellRenderer());
		colorBox.addActionListener(colorAction);
		String[] sexes={"-", "Female", "Male"};
		sexdrop = new JComboBox<String>(sexes);
		sexdrop.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if (e.getActionCommand() == "comboBoxChanged") {  
		    		  if(fromUI) {
		    			
		    			 return; 
		    		  }
		    		 
		    		  JComboBox<String> box = (JComboBox)e.getSource();
		    		  if(box.getSelectedIndex() == 0) {
		    			  sample.female = null;
		    		  }
		    		  else if(box.getSelectedIndex() == 1) {
		    			  sample.female = true;  
		    		  }
		    		  else {
		    			 
		    			  sample.female = false;
		    		  }
				}
				
			}
			
		});
		if(sample.female != null) {			
			if(sample.female) {
				sexdrop.setSelectedIndex(1);
			}
			else {				
				sexdrop.setSelectedIndex(2);
			}
			if(sample.children != null && sample.children.size() > 0) {
				if(sample.equals(sample.children.get(0).mother) || sample.equals(sample.children.get(0).father)) {
					sexdrop.setEnabled(false);
				}
			}		
		}		

		ArrayList<Sample> siblings = getSiblings(sample);
		motherlist.clear();
		motherlist.add(null);
		for(int i = 0 ; i<Main.drawCanvas.sampleList.size(); i++) {
			
			if(Main.drawCanvas.sampleList.get(i).equals(sample)) {
				continue;					
			}
			if(sample.children != null && sample.children.contains(Main.drawCanvas.sampleList.get(i)) ) {
				continue;
			}
			
			if(sample.father != null && sample.father.equals(Main.drawCanvas.sampleList.get(i)) ) {
				continue;
			}
			if(Main.drawCanvas.sampleList.get(i).female != null && !Main.drawCanvas.sampleList.get(i).female) {
				continue;
			}
			if(siblings.contains(Main.drawCanvas.sampleList.get(i))) {
				continue;
			}
			if(Main.drawCanvas.sampleList.get(i).children!=null && Main.drawCanvas.sampleList.get(i).children.size() > 0) {
				if(sample.children != null && sample.children.size() > 0) {
					if(Main.drawCanvas.sampleList.get(i).children.contains(sample.children.get(0))) {
						continue;
					}
				}				
			}
			motherlist.add(Main.drawCanvas.sampleList.get(i));
		}
		motherdrop = new SampleComboBox(motherlist);
			
		motherdrop.setRenderer(comborenderer);
		if(motherdrop.getActionListeners().length == 0) {
			motherdrop.addActionListener(parentDropActionListener);			
		}
		if(sample.mother != null) {
			fromUI = true;
			motherdrop.setSelectedItem(sample.mother);
			fromUI = false;
		}
		
		fatherlist.clear();
		fatherlist.add(null);
		for(int i = 0 ; i<Main.drawCanvas.sampleList.size(); i++) {
			if(Main.drawCanvas.sampleList.get(i).equals(sample)) {
				continue;					
			}
			if(sample.children != null && sample.children.contains(Main.drawCanvas.sampleList.get(i)) ) {
				continue;
			}
			
			if(sample.mother != null && sample.mother.equals(Main.drawCanvas.sampleList.get(i)) ) {
				continue;
			}
			if(Main.drawCanvas.sampleList.get(i).female != null && Main.drawCanvas.sampleList.get(i).female) {
				continue;
			}
			if(siblings.contains(Main.drawCanvas.sampleList.get(i))) {
				continue;
			}
			if(Main.drawCanvas.sampleList.get(i).children!=null && Main.drawCanvas.sampleList.get(i).children.size() > 0) {
				if(sample.children != null && Main.drawCanvas.sampleList.get(i).children.size() > 0) {
					if(Main.drawCanvas.sampleList.get(i).children.contains(sample.children.get(0))) {
						continue;
					}
				}				
			}
			fatherlist.add(Main.drawCanvas.sampleList.get(i));
		}
		fatherdrop = new SampleComboBox(fatherlist);
			
		fatherdrop.setRenderer(comborenderer);
		if(fatherdrop.getActionListeners().length == 0) {
			fatherdrop.addActionListener(parentDropActionListener);			
		}
		
		if(sample.father != null) {
			fromUI = true;
			fatherdrop.setSelectedItem(sample.father);
			fromUI = false;
		}
	}	
	void setPanels() {
		panel.removeAll();		
		panel.add(samplelabel);				
		panel.add(sampleNameField);		
		panel.add(annotation);		
		panel.add(sexlabel);
		panel.add(sexdrop);
		panel.add(affected);	
		panel.add(motherlabel);
		panel.add(motherdrop);
		panel.add(fatherlabel);
		panel.add(fatherdrop);
		
		
		annotationpointer = 0 ;
		for(int i = 0 ; i<panel.getComponentCount(); i++) {
			if(panel.getComponent(i).equals(annotation)) {
				annotationpointer = i;
				break;
			}		
		}
		
		if(affected.getActionListeners().length == 0) {
			affected.addActionListener(this);
			annotation.addActionListener(this);
			intersect.addActionListener(this);
		}
		if(sample.children != null && sample.children.size() > 0) {
			panel.add(new JLabel("Children:"));
			for(int i=0; i<sample.children.size(); i++) {
				if(!Main.drawCanvas.sampleList.contains(sample.children.get(i))) {
					sample.children.remove(i);
					i--;
					continue;
				}
				panel.add(new JLabel(sample.children.get(i).getName()));
			}
		}
		if(annotation.isSelected()) {
			panel.add(intersect,annotationpointer+1);			
		}
		else {
			panel.remove(intersect);			
		}
		ArrayList<Sample> siblings = getSiblings(sample);
		
		
		if(siblings.size() > 0) {
			panel.add(new JLabel("Siblings:"));
			for(int i=0; i<siblings.size(); i++) {
				panel.add(new JLabel(siblings.get(i).getName()));
			}
		}
		panel.add(colorBox);
		if(sample.getTabixFile() !=null) {
			panel.add(new JLabel("VCF path:"));
			vcfpath.setPreferredSize(new Dimension(200, Main.defaultFontSize+4));
			vcfpath.setText(sample.getTabixFile());
			vcfpath.setCaretPosition(0);
			panel.add(vcfpath);
		}
		if(sample.samFile != null) {
			panel.add(new JLabel("BAM path:"));
			bampath.setPreferredSize(new Dimension(200, Main.defaultFontSize+4));
			try {
				bampath.setText(sample.samFile.getCanonicalPath());
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			bampath.setCaretPosition(0);
			panel.add(bampath);
		}
		frame.pack();
		
		 for(int i =0; i<frame.getContentPane().getComponentCount(); i++) {
			 frame.getContentPane().getComponent(i).setMinimumSize(frame.getContentPane().getComponent(i).getPreferredSize());
		 }
	}
	public SampleDialog(Sample sample) {
		super(new GridLayout(0,1));		
		try {
			SampleDialog.sample = sample;
			
			sampleNameField.setText(sample.getName());		
			sampleNameField.addKeyListener(this);
			
			if(sample.affected != null) {
				affected.setSelected(sample.affected);			
			}
			
			annotation.setSelected(sample.annotation);
			intersect.setSelected(sample.intersect);
			
			if(sample.annotation) {
				setComponentsEnabled(true); 
			}
			else {
				setComponentsEnabled(false); 
			}
			setDropboxes();			
			setFonts(Main.menuFont);
			setPanels();
			checkFiles();
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	ArrayList<Sample> getSiblings(Sample sample) {
		ArrayList<Sample> siblings = new ArrayList<Sample>();
		if(sample.father != null && sample.father.children.size() > 1) {
			
			for(int i=0; i<sample.father.children.size(); i++) {
				if(sample.father.children.get(i).equals(sample)) {
					continue;
				}
				if(!Main.drawCanvas.sampleList.contains(sample.father.children.get(i))) {
					sample.father.children.remove(i);
					i--;
					continue;
				}
				siblings.add(sample.father.children.get(i));
				
			}
		}
		if(sample.mother != null && sample.mother.children.size() > 1) {
			
			for(int i=0; i<sample.mother.children.size(); i++) {
				if(sample.mother.children.get(i).equals(sample)) {
					continue;
				}
				if(!Main.drawCanvas.sampleList.contains(sample.mother.children.get(i))) {
					sample.mother.children.remove(i);
					i--;
					continue;
				}
				if(!siblings.contains(sample.mother.children.get(i))) {
					siblings.add(sample.mother.children.get(i));
				}
				
				
			}
		}
		return siblings;
	}
	public void createAndShowGUI() {	
			
			frame.setAlwaysOnTop(true);		 	
		    //JComponent newContentPane = new SampleDialog(sample);
		    this.setOpaque(false);		   
		    frame.setContentPane(panel);
		    setFonts(Main.menuFont);
		    frame.pack();
		    int ylocation = (int)MouseInfo.getPointerInfo().getLocation().getY();
		   
		    if(ylocation + frame.getHeight() > Main.frame.getY() + Main.frame.getHeight()) {
		    	ylocation = (int)(Main.frame.getY() + Main.frame.getHeight() - frame.getHeight() -10);
		    }
		    frame.setMinimumSize(new Dimension(200, 300));
		    frame.setLocation((int)(Main.frame.getLocation().getX()+Main.sidebarWidth+10), ylocation);
		    frame.setVisible(true);
		   
		    for(int i =0; i<frame.getContentPane().getComponentCount(); i++) {
		    	frame.getContentPane().getComponent(i).setMinimumSize(frame.getContentPane().getComponent(i).getPreferredSize());
		    }		  
			
	}
	
	static void setFonts(Font menuFont) {
		for(int i = 0 ; i<SampleDialog.frame.getContentPane().getComponentCount(); i++) {
			SampleDialog.frame.getContentPane().getComponent(i).setFont(menuFont);
		}
		SampleDialog.frame.pack();
	}
	@Override
	public void actionPerformed(ActionEvent event) {		
		if(event.getSource() == affected) {
			
			
			sample.affected = affected.isSelected();
			
			checkFiles();
			Main.drawCanvas.repaint();
		}
		else if(event.getSource() == annotation) {
			
			sample.annotation = annotation.isSelected();
			if(sample.annotation) {
				Main.varsamples--;
				setComponentsEnabled(true); 
				
				panel.add(intersect,annotationpointer+1);			
			}
			else {
				Main.varsamples++;
				setComponentsEnabled(false); 
				panel.remove(intersect);		
				
			}
			panel.repaint();
			panel.revalidate();
			frame.pack();
			updateVarsamples();
			boolean annotation = false;
			boolean intersect = false;
			for(int i = Main.drawCanvas.sampleList.size()-1; i>=0; i--) {
				if(Main.drawCanvas.sampleList.get(i).annotation && Main.drawCanvas.sampleList.get(i).intersect) {
					annotation = true;
					intersect = true;
					break;
				}
				if(Main.drawCanvas.sampleList.get(i).annotation) {
					annotation = true;
					
				}
			}
			if(annotation) {
				checkAnnotation();
			}
			Main.drawCanvas.intersect = intersect;
			Main.drawCanvas.annotationOn = annotation;
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
		}
		else if(event.getSource() == intersect) {
			sample.intersect = intersect.isSelected();
			boolean annotation = false;
			boolean intersect = false;
			for(int i = Main.drawCanvas.sampleList.size()-1; i>=0; i--) {
				if(Main.drawCanvas.sampleList.get(i).annotation && Main.drawCanvas.sampleList.get(i).intersect) {
					annotation = true;
					intersect = true;
					break;
				}
				if(Main.drawCanvas.sampleList.get(i).annotation) {
					annotation = true;
					
				}
			}
			Main.drawCanvas.intersect = intersect;
			Main.drawCanvas.annotationOn = annotation;
			Draw.updatevars = true;
			Main.drawCanvas.repaint();
		}
	}
static void checkAnnotation() {
	VarNode node = FileRead.head.getNext();
	Entry<String, ArrayList<SampleNode>> entry;
	
	while(node != null) {
		boolean found = false;
		for(int v = 0 ; v<node.vars.size(); v++) {			
			entry = node.vars.get(v);			
			for(int m = 0; m<entry.getValue().size(); m++) {
				if(entry.getValue().get(m).getSample() == null) {							
					continue;				
				}
				if(!entry.getValue().get(m).getSample().annotation) {
					found = true;
					break;
				}
			}
			if(found) {
				break;
			}
		}
		
		
		if(!found) {
			
			node.annotationOnly = true;
		
		}
		else {
			node.annotationOnly = false;
		}
		node = node.getNext();
	}	
	
}
void setComponentsEnabled(boolean enabled) {
	
	if(enabled) {
		
		for(int i = annotationpointer+1 ; i<panel.getComponentCount(); i++) {
			panel.getComponent(i).setEnabled(false);
			
		}
	}
	else {
		
		for(int i = annotationpointer+1 ; i<panel.getComponentCount(); i++) {
			panel.getComponent(i).setEnabled(true);			
		}
	}
	
}
void updateVarsamples() {
	VariantHandler.commonSlider.setMaximum(Main.varsamples);
	VariantHandler.geneSlider.setMaximum(Main.varsamples);
	VariantHandler.commonSlider.setUpperValue(Main.varsamples);
	
}
@Override
public void keyPressed(KeyEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void keyReleased(KeyEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void keyTyped(KeyEvent arg0) {
	sample.setName(sampleNameField.getText() +arg0.getKeyChar());
	Main.drawCanvas.repaint();
	
}
class MyCellRenderer extends JButton implements ListCellRenderer {  
   
	private static final long serialVersionUID = 1L;
	public MyCellRenderer() {  
        setOpaque(true);
    }
    boolean b=false;
   @Override
   public void setBackground(Color bg) {
       
        if(!b)
        {
            return;
        }

       super.setBackground(bg);
   }
    public Component getListCellRendererComponent(  
        JList list,  
        Object value,  
        int index,  

        boolean isSelected,  
        boolean cellHasFocus)  
    {  
    	
        b=true;
      
        setText(" ");           
        setBackground((Color)value);        
        b=false;
        return this;  
    }  
}
}


	