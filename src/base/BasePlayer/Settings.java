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
import java.awt.Font;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.HashMap;

import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JColorChooser;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JSlider;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class Settings  extends JPanel implements ActionListener, ChangeListener, MouseListener, KeyListener {

	private static final long serialVersionUID = 1L;
	private static JSlider insertSizeSlider = new JSlider(0,10000);
	private static JLabel insertLabel = new JLabel("Maximum insert size: 1000"), readLabel = new JLabel("Read options");
	private static Dimension mindimension = new Dimension(200, 15);
	private static JLabel depthLimitLabel;
	private static JLabel coverageDistanceLabel = new JLabel("Coverage draw distance (bp):");
	
	private static JSlider depthLimitSlide = new JSlider(1,10000);
	private static JSlider rSlider = new JSlider(60,255);
	private static JSlider gSlider = new JSlider(60,255);
	private static JSlider bSlider = new JSlider(60,255);
	static SteppedComboBox fontlist, varDrawList;
    static DefaultComboBoxModel<String> fontModel, varModel;
	private static JSlider mappingQuality = new JSlider(0,60), baseQuality = new JSlider(0,60);
	private static JLabel mappingLabel = new JLabel("Mapping quality: 10"), baseLabel = new JLabel("Base quality: 10"), reloadReads = new JLabel("Click here to reload reads"); 
	
	private static JLabel rLabel = new JLabel("Red: ");
	private static JLabel gLabel = new JLabel("Green: ");
	private static JLabel bLabel = new JLabel("Blue: ");
	static JFrame frame = new JFrame("Settings");
	private static JCheckBox softclips = new JCheckBox("Show bases in softclips");
	private static JPanel readPanel = new JPanel(new GridBagLayout());
	private static JPanel varPanel = new JPanel(new GridBagLayout());
	private static JPanel generalPanel = new JPanel(new GridBagLayout());
	private static JPanel appearancePanel = new JPanel(new GridBagLayout());
	private static JTabbedPane tabPanel = new JTabbedPane(JTabbedPane.LEFT);
	private static JButton colorButton = new JButton("Select color");
	static HashMap<String, Integer> settings;
	private static JLabel windowLabel = new JLabel("Processing window size (bp):");
	
	private static JLabel bigFileLabel = new JLabel("Big file size for tracks (MB):");
	private static JTextField coverageDistanceField = new JTextField("");
	private static JTextField bigFileField = new JTextField("");
	private static JTextField windowField = new JTextField("");
	static int readDrawDistance, coverageDrawDistance,coverageAlleleFreq, windowSize, readDepthLimit, softClips, mappingQ, insertSize, baseQ;
	static Color frameColor = new Color(188,188,178);
	
	static JCheckBox bold = new JCheckBox("Bold");
	static JColorChooser colorchooser = new JColorChooser();
	static boolean constr = true;
	static int selectedVarDraw = 0;
	
	public Settings() {
		super(new GridBagLayout());		
		try {
		this.setBackground(Color.black);
		windowLabel.setToolTipText("Window size for processing large vcf and bed files.");
		bigFileLabel.setToolTipText("Maximum file size for chromosomal level drawing.");
		depthLimitLabel = new JLabel("");
		String[] sizes = {"8","10","12","14","16","18","20","22","24"};
		
		fontModel = new DefaultComboBoxModel<String>(sizes);
		fontlist = new SteppedComboBox(fontModel);		
		fontlist.addActionListener(this);				
		fontlist.setBorder(BorderFactory.createLineBorder(Color.lightGray, 1));
		fontlist.setBackground(Color.white);
		fontlist.setEditable(true);
		fontlist.setSelectedItem(""+Main.defaultFontSize);
		constr = false;
		String[] varTypes = {"Coverage","Allelic fraction","Quality","GQ","Calls"};
		
		varModel = new DefaultComboBoxModel<String>(varTypes);
		varDrawList = new SteppedComboBox(varModel);
		//varDrawList.setEditable(true);
		varDrawList.addActionListener(this);		
		
		varDrawList.setBorder(BorderFactory.createLineBorder(Color.lightGray, 1));
		
		
		reloadReads.setForeground(Color.red);
		reloadReads.setVisible(false);
		reloadReads.setMinimumSize(mindimension);
		
		
		reloadReads.addMouseListener(this);
		insertSizeSlider.setMinimumSize(mindimension);
		readPanel.add(insertSizeSlider);
		
		insertSizeSlider.setOpaque(false);
		mappingQuality.setMinimumSize(mindimension);
		
		mappingQuality.setOpaque(false);
		baseQuality.setOpaque(false);
		baseQuality.setMinimumSize(mindimension);
		insertSizeSlider.addChangeListener(this);
		mappingQuality.addChangeListener(this);
		baseQuality.addChangeListener(this);		
		depthLimitSlide.addChangeListener(this);		
		depthLimitSlide.setOpaque(false);
		softclips.addActionListener(this);
		softclips.setOpaque(false);
		
		depthLimitLabel.setMinimumSize(mindimension);
		rSlider.addChangeListener(this);	
		gSlider.addChangeListener(this);	
		bSlider.addChangeListener(this);	
		setValues();
		rSlider.setOpaque(false);
		gSlider.setOpaque(false);
		bSlider.setOpaque(false);
		colorButton.setPreferredSize(new Dimension(Main.buttonWidth, Main.buttonHeight));
		colorButton.addActionListener(this);
	//	coverageDistanceField.setOpaque(false);
		coverageDistanceField.addKeyListener(this);
		bigFileField.addKeyListener(this);
	//	windowField.setOpaque(false);
		windowField.addKeyListener(this);
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.NORTHWEST;		
		c.insets = new Insets(2,0,0,30);
		c.gridx = 0;
		c.gridy = 0;
		c.gridwidth = 1;	
		readLabel.setName("header");
		readPanel.add(readLabel,c);
				
		c.gridx = 1;
		readPanel.add(reloadReads,c);
		c.gridy++;	
		c.gridx = 0;
		
		c.gridwidth = 2;
		readPanel.add(new JSeparator(),c);
		c.gridy++;
		c.gridwidth = 1;
		readPanel.add(insertSizeSlider, c);
		c.gridx = 1;
		readPanel.add(insertLabel,c);	
		c.gridy++;
		c.gridx = 0;
		readPanel.add(mappingQuality,c);
		c.gridx = 1;
		readPanel.add(mappingLabel,c);
		c.gridx = 0;
		c.gridy++;
		readPanel.add(baseQuality,c);
		c.gridx = 1;
		readPanel.add(baseLabel,c);
		c.gridx = 0;
		c.gridy++;
		readPanel.add(depthLimitSlide,c);
		c.gridx = 1;
		readPanel.add(depthLimitLabel,c);
		c.gridx = 0;
		c.gridy++;
		readPanel.add(softclips,c);
		
		c.weightx = 1;
		c.weighty = 1;
		c.gridwidth = GridBagConstraints.REMAINDER;
		readPanel.add(new JLabel(),c);
		c.gridx = 0;
		c.gridy = 0;
		c.weightx = 0;
		c.weighty = 0;
		c.gridwidth = 3;
		JLabel appearLabel = new JLabel("Appearance");
		appearLabel.setName("header");
		appearancePanel.add(appearLabel,c);
		c.gridy++;
		
		appearancePanel.add(new JSeparator(),c);
		c.gridy++;
		c.gridwidth = 1;
		//c.gridwidth = 3;
		appearancePanel.add(new JLabel("Font size: "),c);
		c.gridx = 1;
		bold.setOpaque(false);
		bold.addActionListener(this);
		appearancePanel.add(bold,c);
		
		c.gridx = 2;
		appearancePanel.add(fontlist, c);
		c.gridwidth = 3;		
		c.gridx = 0;
		c.gridy++;
		
		appearancePanel.add(new JSeparator(),c);
		c.gridy++;
		
		JLabel colors = new JLabel("Background color");
		c.gridwidth = 1;	
		colors.setName("header");
		appearancePanel.add(colors,c);
		c.gridy++;
		appearancePanel.add(rSlider,c);
		c.gridx = 2;
		appearancePanel.add(rLabel,c);
		c.gridy++;
		c.gridx = 0;
		appearancePanel.add(gSlider,c);
		c.gridx = 2;
		appearancePanel.add(gLabel,c);
		c.gridy++;
		c.gridx = 0;
		appearancePanel.add(bSlider,c);
		c.gridx = 2;
		appearancePanel.add(bLabel,c);
		c.weightx = 1;
		c.weighty = 1;
		c.gridwidth = GridBagConstraints.REMAINDER;
		appearancePanel.add(new JLabel(""),c);
		c.gridx = 0;
		c.gridy = 0;
		c.weightx = 0;
		c.weighty = 0;
		c.gridwidth = 1;
		JLabel generalLabel = new JLabel("General options");
		generalLabel.setName("header");
		
		generalPanel.add(generalLabel,c);		
		c.gridy++;	
		
		c.gridwidth = 2;
		generalPanel.add(new JSeparator(),c);
		c.gridwidth = 1;		
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridy++;
		c.gridx = 0;
		generalPanel.add(coverageDistanceLabel,c);		
		c.gridx = 1;
		generalPanel.add(coverageDistanceField,c);
		c.gridy++;
		c.gridx = 0;
		generalPanel.add(windowLabel,c);
		
		c.gridx = 1;
		generalPanel.add(windowField,c);
		c.gridx = 0;
		c.gridy++;
		generalPanel.add(bigFileLabel,c);
		c.gridx = 1;
		generalPanel.add(bigFileField,c);
		
		c.gridx = 0;
		c.weightx = 1;
		c.weighty = 1;
		c.gridwidth = GridBagConstraints.REMAINDER;
		generalPanel.add(new JLabel(),c);
		
		c.gridx = 0;
		c.gridy = 0;
		c.weightx = 0;
		c.weighty = 0;
		c.gridwidth = 1;
		JLabel varLabel = new JLabel("Variant options");
		varLabel.setName("header");
		varPanel.add(varLabel, c);
		c.gridy++;
		c.gridwidth = 2;
		varPanel.add(new JSeparator(),c);
		c.gridwidth = 1;	
		varPanel.add(new JLabel("Variant height by: "), c);
		c.gridx = 1;
		varDrawList.setBackground(Color.white);
		varPanel.add(varDrawList, c);
		c.gridx = 0; 
		c.weightx = 1;
		c.gridy++;
		c.weighty = 1;
		c.gridwidth = GridBagConstraints.REMAINDER;
		varPanel.add(new JLabel(),c);
		
		tabPanel.add("General", generalPanel);
		tabPanel.add("Variants", varPanel);
		tabPanel.add("Reads", readPanel);
		tabPanel.add("Appearance", appearancePanel);
		
		//tabPanel.setBackground(color);
		/*generalPanel.setBackground(Draw.sidecolor);
		readPanel.setBackground(Draw.sidecolor);
		varPanel.setBackground(Draw.sidecolor);
		appearancePanel.setBackground(Draw.sidecolor);
		*///frame.getContentPane().setBackground(color);
		add(tabPanel);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		setFonts(Main.menuFont);
		
	//	setBackground(color);		
	}
	 protected void paintComponent(Graphics g) {
	        super.paintComponent(g);     
	        
	     //   g.drawImage(image, 0, 0, this.getWidth(), this.getHeight(), null);
	        g.setColor(frameColor);
	        g.fillRect(0, 0, this.getWidth(), this.getHeight());	        
	       
	    }
	static void setFonts(Font menuFont) {
		JPanel panel = null;
		tabPanel.setFont(menuFont);
		if(Main.drawCanvas == null) {
			return;
		}
		for(int i = 0 ; i<tabPanel.getComponentCount(); i++) {
			
			if(tabPanel.getComponent(i) instanceof JPanel) {
				
				panel = (JPanel)tabPanel.getComponent(i);
			//	panel.setFont(menuFont);
				for(int c = 0; c < panel.getComponentCount(); c++) {
					if(panel.getComponent(c).getName() != null) {
						panel.getComponent(c).setFont(Main.menuFontBold.deriveFont((float)(Main.defaultFontSize+2)));
						
					}
					else {
						panel.getComponent(c).setFont(menuFont);
					}
					
				}
			}
			
		}
		frame.pack();
		panel = null;
		
	}
	static void setValues() {
		if(settings == null) {
			settings = new HashMap<String, Integer>();			
			settings.put("readDrawDistance", 60000);
			settings.put("baseQ", 10);
			settings.put("readDepthLimit", 1000);
			settings.put("coverageDrawDistance", 1000000);
			settings.put("coverageAlleleFreq", 1);
			settings.put("windowSize", 1000000);
			settings.put("bigFile", 200);
			settings.put("rValue", 228);
			settings.put("gValue", 228);
			settings.put("bValue", 218);
			settings.put("insertSize", 1000);
			settings.put("mappingQuality", 10);
			settings.put("softClips", 0);		
			
		}
		
		baseQ = settings.get("baseQ");
		readDrawDistance = settings.get("readDrawDistance");
		coverageDrawDistance = settings.get("coverageDrawDistance");
		coverageAlleleFreq = settings.get("coverageAlleleFreq");
		coverageDistanceField.setText(""+settings.get("coverageDrawDistance"));
		bigFileField.setText(""+settings.get("bigFile"));
		depthLimitLabel.setText("Read depth limit: " +settings.get("readDepthLimit"));		
		windowField.setText(""+settings.get("windowSize"));
		baseQuality.setValue(settings.get("baseQ"));
	//	rSlider.setValue(settings.get("colorValue"));
		if(settings.get("rValue") == null) {
			settings.put("rValue", 228);
			settings.put("gValue", 228);
			settings.put("bValue", 218);
		}
		
		rSlider.setValue(settings.get("rValue"));
		gSlider.setValue(settings.get("gValue"));
		bSlider.setValue(settings.get("bValue"));
		
		rLabel.setText("Red: " +rSlider.getValue());
		gLabel.setText("Green: " +gSlider.getValue());
		bLabel.setText("Blue: " +bSlider.getValue());
		
		insertSizeSlider.setValue(settings.get("insertSize"));
		mappingQuality.setValue(settings.get("mappingQuality"));
		mappingQ = settings.get("mappingQuality");
		depthLimitSlide.setValue(settings.get("readDepthLimit"));
		insertLabel.setText("Maximum insert size: "+settings.get("insertSize"));
		windowSize = settings.get("windowSize");
		readDepthLimit = settings.get("readDepthLimit");
		softClips = settings.get("softClips");
		if(softClips == 1) {
			softclips.setSelected(true);
		}
		
		reloadReads.setVisible(false);
	}
	
	private static void createAndShowGUI() {	
		JFrame.setDefaultLookAndFeelDecorated(false);
		if(Main.userDir == null) {
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
			 frame.setVisible(true); 
		}
		else {
		 frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE); 
		 frame.setVisible(false);  
		}
		JFrame.setDefaultLookAndFeelDecorated(false);
		 frame.setResizable(false);    
		 frame.setAlwaysOnTop(true);
		
	    JComponent newContentPane = new Settings();
	    newContentPane.setOpaque(false); 
	   
	    frame.setContentPane(newContentPane);
	    frame.pack();
	   
	}
	public static void main(String[] args) {
		
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	          public void run() {	          
	          		createAndShowGUI();
	          	
	          }
	      });	
	}
	public void actionPerformed(ActionEvent event) {
		if(event.getSource() == bold) {
			Main.setFonts();
		}
		else if(event.getSource() == softclips) {
			if(softclips.isSelected()) {
				softClips = 1;
				settings.put("softClips", 1);
			}
			else {
				softClips = 0;
				settings.put("softClips", 0);
			}
			reloadReads.setVisible(true);
			return;
		}
		if(event.getSource() == fontlist && !constr) {
			if (event.getActionCommand().equals("comboBoxEdited")) {
				try {
					Main.defaultFontSize = Integer.parseInt(fontlist.getSelectedItem().toString().trim());
					Main.setFonts();
					Main.writeToConfig("fontSize=" +Main.defaultFontSize);
				}
				catch(Exception e) {
					
				}		    
			}
			if (event.getActionCommand().equals("comboBoxChanged")) {
				Main.defaultFontSize = Integer.parseInt(fontlist.getSelectedItem().toString());
				if(Main.drawCanvas != null) {
					Main.setFonts();
				}
				Main.writeToConfig("fontSize=" +Main.defaultFontSize);
			}
		}
		if(event.getSource() == varDrawList) {
			if (event.getActionCommand().equals("comboBoxChanged")) {
				selectedVarDraw = varDrawList.getSelectedIndex();
				switch(Settings.selectedVarDraw) {	
				
					case 0: {	
						VariantHandler.maxSlideValue = VariantHandler.maxCoverageSlider.getValue();
						break;
					}
					case 1: {
						VariantHandler.maxSlideValue = VariantHandler.callSlider.getUpperValue()/(float)100;
						break;
					}
					case 4: {
						VariantHandler.maxSlideValue = VariantHandler.maxCoverageSlider.getValue();
						break;
					}
					
				}
				Main.chromosomeDropdown.setSelectedIndex(Main.selectedChrom);
			}
		}
		
	}
	@Override
	public void stateChanged(ChangeEvent event) {
		if(event.getSource() == insertSizeSlider) {
			insertLabel.setText("Maximum insert size: " +insertSizeSlider.getValue());
			insertSize = insertSizeSlider.getValue();
			settings.put("insertSize", insertSizeSlider.getValue());
			reloadReads.setVisible(true);
			return;
		}
		else if(event.getSource() == rSlider || event.getSource() == gSlider || event.getSource() == bSlider) {
			
			rLabel.setText("Red: " +rSlider.getValue());
			gLabel.setText("Green: " +gSlider.getValue());
			bLabel.setText("Blue: " +bSlider.getValue());
			
			Draw.sidecolor = new Color(rSlider.getValue(), gSlider.getValue(), bSlider.getValue());				
			if(VariantHandler.frame != null) {
				VariantHandler.backColor = new Color(rSlider.getValue(), gSlider.getValue(), bSlider.getValue());
				VariantHandler.adder.setForeground(Draw.sidecolor);
				VariantHandler.adder2.setForeground(Draw.sidecolor);
				VariantHandler.adder3.setForeground(Draw.sidecolor);
				VariantHandler.adder4.setForeground(Draw.sidecolor);
				//VariantHandler.adder5.setForeground(Draw.sidecolor);
			}
			/*if(colorSlider.getValue() <= 120) {
				Draw.sidecolor = new Color(200, 80+colorSlider.getValue(), 120);				
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
			}
			
			else if(colorSlider.getValue() <= 240) {
				Draw.sidecolor = new Color(200-(colorSlider.getValue()-120), 200, 120);
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
				
			}
			else if(colorSlider.getValue() <= 360) {
				Draw.sidecolor = new Color(120, 200, 80+(colorSlider.getValue()-240));
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
				
			}
			else if(colorSlider.getValue() <= 480) {
				Draw.sidecolor = new Color(120, 200-(colorSlider.getValue()-360), 200);
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
				
			}
			else if(colorSlider.getValue() <= 600) {
				Draw.sidecolor = new Color(80+(colorSlider.getValue()-480),120, 200);	
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
			}
			else {
				Draw.sidecolor = new Color(200, 120, 200-(colorSlider.getValue()-600));	
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
			}
			
			*/
			
			//rLabel.setText("Red: " +settings.get("colorValue"));
			if(event.getSource() == rSlider) {
				settings.put("rValue",rSlider.getValue());
			}
			else if (event.getSource() == gSlider){
				settings.put("gValue",gSlider.getValue());
			}
			else {
				settings.put("bValue",bSlider.getValue());
			}
			
			appearancePanel.setBackground(Draw.sidecolor);
			
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
				Main.chromDraw.repaint();
				Main.bedCanvas.repaint();
				Main.panel.setBackground(new Color(Draw.sidecolor.getRed()-40, Draw.sidecolor.getGreen()-40, Draw.sidecolor.getBlue()-40));
				Main.chrompan.setBackground(Draw.sidecolor);
				Main.setbut.setBackground(Main.panel.getBackground());
				Main.panel.revalidate();
				if(VariantHandler.frame != null) {
					VariantHandler.frame.getContentPane().setBackground(new Color(Draw.sidecolor.getRed()-40, Draw.sidecolor.getGreen()-40, Draw.sidecolor.getBlue()-40));
					VariantHandler.filterpanel.setBackground(VariantHandler.backColor);
					VariantHandler.filterpanel.revalidate();
					VariantHandler.aminopanel.setBackground(VariantHandler.backColor);
					VariantHandler.aminopanel.revalidate();
					VariantHandler.comparepanel.setBackground(VariantHandler.backColor);
					VariantHandler.comparepanel.revalidate();
					VariantHandler.tabs.setBackground(VariantHandler.backColor);
					VariantHandler.tabs.revalidate();
				}
			}
		}
		
		if(event.getSource() == mappingQuality) {
			
			mappingLabel.setText("Mapping quality: " +mappingQuality.getValue());
			mappingQ = mappingQuality.getValue();
			settings.put("mappingQuality", mappingQuality.getValue());
			Draw.updateReads = true;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
			}
			return;			
		}
		if(event.getSource() == baseQuality) {
			settings.put("baseQ", baseQuality.getValue());
			baseLabel.setText("Base quality: " +baseQuality.getValue());
			baseQ = baseQuality.getValue();
			Draw.updateReads = true;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
			}
			return;
		}
		if(event.getSource() == depthLimitSlide) {
			settings.put("readDepthLimit",depthLimitSlide.getValue());
			readDepthLimit = depthLimitSlide.getValue();
			depthLimitLabel.setText("Read depth limit: " +depthLimitSlide.getValue());
			
			return;
		}
		
		/*if(event.getSource() == coverageDistanceSlide) {
			coverageDrawDistance = coverageDistanceSlide.getValue();
			coverageDistanceLabel.setText("Coverage draw distance: " +coverageDrawDistance +"bp");
			
			return;
		}*/
		
		
	}
	@Override
	public void mouseClicked(MouseEvent arg0) {
		
		
	}
	@Override
	public void mouseEntered(MouseEvent event) {
		if(event.getSource() == reloadReads) {
			if(getCursor().getType() != Cursor.HAND_CURSOR ) {
				setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));				
			}
		}
		
	}
	@Override
	public void mouseExited(MouseEvent event) {
		if(event.getSource() == reloadReads) {
			if(getCursor().getType() != Cursor.DEFAULT_CURSOR ) {
				setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));				
			}
		}
		
	}
	@Override
	public void mousePressed(MouseEvent event) {
		if(event.getSource() == reloadReads) {
			Main.drawCanvas.clearReads();
			for(int i = 0; i<Main.drawCanvas.splits.size(); i++) {
				
				Main.drawCanvas.splits.get(i).updateReads = true;
				Main.drawCanvas.drawReads(Main.drawCanvas.splits.get(i));
			}
			Main.drawCanvas.repaint();
			reloadReads.setForeground(Color.black);
		}
		
	}
	@Override
	public void mouseReleased(MouseEvent event) {
		if(event.getSource() == reloadReads) {
			reloadReads.setVisible(false);
			reloadReads.setForeground(Color.red);
		}
		
	}
	@Override
	public void keyTyped(KeyEvent e) {
		
	}
	@Override
	public void keyPressed(KeyEvent e) {
		if(e.getSource() == coverageDistanceField) {
			try {
				settings.put("coverageDrawDistance", Integer.parseInt(coverageDistanceField.getText()));
				coverageDrawDistance = Integer.parseInt(coverageDistanceField.getText());
			}
			catch(Exception ex) {
				settings.put("coverageDrawDistance", Integer.MAX_VALUE);
				coverageDrawDistance= Integer.MAX_VALUE;
				
			}
		} 
		else if(e.getSource() == windowField) {
			try {
				settings.put("windowSize", Integer.parseInt(windowField.getText()));
				windowSize = Integer.parseInt(windowField.getText());
			}
			catch(Exception ex) {
				settings.put("windowSize", Integer.MAX_VALUE);
				windowSize = Integer.MAX_VALUE;
				
			}
		} 
		else if(e.getSource() == bigFileField) {
			try {
				settings.put("bigFile", Integer.parseInt(bigFileField.getText()));
				for(int i = 0 ; i< Main.bedCanvas.bedTrack.size(); i++) {
					if(Main.bedCanvas.bedTrack.get(i).file.length() / 1048576 < Settings.settings.get("bigFile")) {
						Main.bedCanvas.bedTrack.get(i).small = true;			    	      	 
				    }
					else {
						Main.bedCanvas.bedTrack.get(i).small = false;	
					}
				}				
			}
			catch(Exception ex) {
				settings.put("bigFile", 200);				
			}
		}
	}
	@Override
	public void keyReleased(KeyEvent e) {
		
	}
	
	
}
