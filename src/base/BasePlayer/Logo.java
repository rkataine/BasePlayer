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
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.Image;
import java.awt.Toolkit;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class Logo extends JPanel {

	static JFrame frame = new JFrame("Launcher");    
	static JLabel image = new JLabel();
	private static final long serialVersionUID = 1L;
	static Image iconimage;
	
	public Logo(){
		
		    iconimage=Toolkit.getDefaultToolkit().getImage(getClass().getResource("Logo.png"));
			ImageIcon icon = new ImageIcon(iconimage); 
			image.setIcon(icon);
			frame.setIconImage(iconimage);
			frame.add(image);
			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));	
		
	}
	
	public static void main(String[] args) {
		
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	         

			public void run() {
	      //  	   timer = System.currentTimeMillis();
	           		createAndShowGUI();
	           		
	           	
	           }
	       });    
	   
	}
	private static void createAndShowGUI() {	
		   try {
		//	   image.addMouseMotionListener(this);
				
				
		JFrame.setDefaultLookAndFeelDecorated(false);
		//Logo logo = new Logo();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
	   frame.setAlwaysOnTop(true);
	   int width = Toolkit.getDefaultToolkit().getScreenSize().width;
	   int height = Toolkit.getDefaultToolkit().getScreenSize().height;
	   try {
		GraphicsDevice gd = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice();
		width = gd.getDisplayMode().getWidth();
		height = gd.getDisplayMode().getHeight();
	   }
	   catch(Exception e) {
		   e.printStackTrace();
	   }
		 frame.setPreferredSize(new Dimension( iconimage.getWidth(null),iconimage.getHeight(null)));
		 frame.setLocation((int)width/2-300, (int)height/2-300);
		 frame.setSize(iconimage.getWidth(null),iconimage.getHeight(null));
		 frame.setResizable(false);    
	   frame.setUndecorated(true);
	   
	 //   JComponent newContentPane = new Logo();
	  //  newContentPane.setOpaque(false);			
	    image.setOpaque(false);
	   
//	    frame.setContentPane(newContentPane);
	    try {
	    	frame.setBackground(new Color(1.0f,1.0f,1.0f,0.0f));
	    }
	    catch(Exception e) {
	    	frame.setBackground(new Color(1.0f,1.0f,1.0f));
	    }
	    frame.setUndecorated(true);
		
	    frame.setVisible(true); 
	    //frame.pack();
		   }
		   catch(Exception e) {
			   e.printStackTrace();
		   }
	   
	}
}
