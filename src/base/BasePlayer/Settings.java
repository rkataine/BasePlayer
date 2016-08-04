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
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class Settings  extends JPanel implements ActionListener, ChangeListener, MouseListener, KeyListener{

	private static final long serialVersionUID = 1L;
	static JSlider insertSize = new JSlider(0,2000);
	static JLabel insertLabel = new JLabel("Maximum insert size: 1000"), readLabel = new JLabel("Read filters:");
	static JLabel depthLimitLabel;
	static JLabel coverageDistanceLabel = new JLabel("Coverage draw distance");
	static JTextField coverageDistanceSlide = new JTextField("1000000");
	static JSlider depthLimitSlide = new JSlider(1,10000);
	static JSlider colorSlider = new JSlider(0,720);
	static JSlider mappingQuality = new JSlider(0,60), baseQuality = new JSlider(0,60);
	static JLabel mappingLabel = new JLabel("Mapping quality: 10"), baseLabel = new JLabel("Base quality: 10"), reloadReads = new JLabel("Click here to reload reads"); 
	static JLabel colorLabel = new JLabel("Color slider");
	static JFrame frame = new JFrame("Settings");
	JPanel readPanel = new JPanel(new GridLayout(9,2));
	static int readDrawDistance = 60000, readDepthLimit = 1000, coverageDrawDistance = 1000000, coverageAlleleFreq = 1;
	public static boolean softClips = true;
	
	static int getMaxInsertSize() {
		return insertSize.getValue();
	}
	static int getMappingQuality() {
		return mappingQuality.getValue();
	}
	static int getBaseQuality() {
		return baseQuality.getValue();
	}
	public Settings() {
		super(new GridLayout());		
		this.setBackground(Color.black);
		depthLimitLabel = new JLabel("Read depth limit: " +readDepthLimit);
		frame.setPreferredSize(new Dimension(400,200));
		readPanel.setPreferredSize(new Dimension(400,200));
		readPanel.add(readLabel);
		readPanel.add(new JLabel());
		readPanel.add(new JSeparator());
		readPanel.add(new JSeparator());
		reloadReads.setForeground(Color.red);
		reloadReads.setVisible(false);
		readPanel.add(new JLabel());
		readPanel.add(reloadReads);
		reloadReads.addMouseListener(this);
		readPanel.add(insertSize);
		insertSize.setValue(1000);
		insertSize.setOpaque(false);
		mappingQuality.setValue(10);
		mappingQuality.setOpaque(false);
		baseQuality.setOpaque(false);
		insertSize.addChangeListener(this);
		mappingQuality.addChangeListener(this);
		baseQuality.addChangeListener(this);
		baseQuality.setValue(10);
		depthLimitSlide.addChangeListener(this);
		depthLimitSlide.setValue(readDepthLimit);
		depthLimitSlide.setOpaque(false);
		colorSlider.addChangeListener(this);
		colorSlider.setValue(138);
		coverageDistanceSlide.setOpaque(false);
		coverageDistanceSlide.addKeyListener(this);
		readPanel.add(insertLabel);		
		readPanel.add(mappingQuality);
		readPanel.add(mappingLabel);
		readPanel.add(baseQuality);
		readPanel.add(baseLabel);
		readPanel.add(depthLimitSlide);
		readPanel.add(depthLimitLabel);
		readPanel.add(colorSlider);
		readPanel.add(colorLabel);
		colorSlider.setOpaque(false);
		readPanel.add(coverageDistanceSlide);
		readPanel.add(coverageDistanceLabel);
		add(readPanel);
		readPanel.setBackground(Color.white);
		setBackground(Color.white);
	}
	
	private static void createAndShowGUI() {	
	
		if(Main.userDir == null) {
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
			 frame.setVisible(true); 
		}
		else {
		 frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE); 
		 frame.setVisible(false);  
		}
		 frame.setResizable(true);    
		
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
	public void actionPerformed(ActionEvent arg0) {		
		
	}
	@Override
	public void stateChanged(ChangeEvent event) {
		if(event.getSource() == insertSize) {
			insertLabel.setText("Maximum insert size: " +insertSize.getValue());
			reloadReads.setVisible(true);
			return;
		}
		if(event.getSource() == colorSlider) {
			if(colorSlider.getValue() <= 120) {
				Draw.sidecolor = new Color(200, 80+colorSlider.getValue(), 80);				
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
			}
			
			else if(colorSlider.getValue() <= 240) {
				Draw.sidecolor = new Color(200-(colorSlider.getValue()-120), 200, 80);
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
				
			}
			else if(colorSlider.getValue() <= 360) {
				Draw.sidecolor = new Color(80, 200, 80+(colorSlider.getValue()-240));
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
				
			}
			else if(colorSlider.getValue() <= 480) {
				Draw.sidecolor = new Color(80, 200-(colorSlider.getValue()-360), 200);
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
				
			}
			else if(colorSlider.getValue() <= 600) {
				Draw.sidecolor = new Color(80+(colorSlider.getValue()-480), 80, 200);	
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
			}
			else {
				Draw.sidecolor = new Color(200, 80, 200-(colorSlider.getValue()-600));	
				VariantHandler.backColor = new Color(Draw.sidecolor.getRed(), Draw.sidecolor.getGreen(), Draw.sidecolor.getBlue(), 200);
			}
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
				Main.chromDraw.repaint();
				Main.bedCanvas.repaint();
				Main.panel.setBackground(Draw.sidecolor);
				Main.panel.revalidate();
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
		if(event.getSource() == mappingQuality) {
			
			mappingLabel.setText("Mapping quality: " +mappingQuality.getValue());
			Draw.updateReads = true;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
			}
			return;			
		}
		if(event.getSource() == baseQuality) {
			
			baseLabel.setText("Base quality: " +baseQuality.getValue());
			Draw.updateReads = true;
			if(Main.drawCanvas != null) {
				Main.drawCanvas.repaint();
			}
			return;
		}
		if(event.getSource() == depthLimitSlide) {
			readDepthLimit = depthLimitSlide.getValue();
			depthLimitLabel.setText("Read depth limit: " +readDepthLimit);
			
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
		if(e.getSource() == coverageDistanceSlide) {
			try {
				Settings.coverageDrawDistance = Integer.parseInt(coverageDistanceSlide.getText());
			}
			catch(Exception ex) {
				Settings.coverageDrawDistance = Integer.MAX_VALUE;
			}
		} 
		
	}
	@Override
	public void keyReleased(KeyEvent e) {
		
	}
	
	
}
