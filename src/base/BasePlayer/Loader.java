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
import java.awt.Cursor;
import java.awt.Graphics;

import javax.swing.JOptionPane;
import javax.swing.SwingWorker;



public class Loader extends SwingWorker<String, Object>  {
	String text = "Loading";
	int loadBar = 0;
	Graphics g = Main.frame.getGlassPane().getGraphics();
	static boolean memory = false;
	boolean updater = false;
	Loader(String text) {
		
		if(!text.equals("note")) {
			Main.frame.setCursor( Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		}
		memory = false;
		Main.drawCanvas.loadingtext = text;
		Main.drawCanvas.loading = true;
		
		Draw.setGlasspane(true);
		
	}

	@Override
	protected String doInBackground() throws Exception {
		try {
				
				while(Main.drawCanvas.loading) {	
				
					if(memory) {
						String loadtext = Main.drawCanvas.loadingtext;
						Main.cancel(1);
						Main.drawCanvas.variantsStart = 0;
						Main.drawCanvas.variantsEnd = 1;
						memory = false;
						
						if(loadtext.contains("Annotating")) {
							JOptionPane.showMessageDialog(Main.chromDraw, "Memory allocation exceeded in " +loadtext +"\nTip1: Goto \"Tools > Settings > General\" and change 'Processing window size (bp)' to 100,000.\nTip2: try to increase memory allocation for BasePlayer\n", "Information", JOptionPane.INFORMATION_MESSAGE);
							
						}
						else {
							JOptionPane.showMessageDialog(Main.chromDraw, "Try to search more specific region to visualize all data\nTip: try to increase memory allocation for BasePlayer (config.txt in your BasePlayer folder)\nNote: Variant annotation is still possible through Variant handler.", "Information", JOptionPane.INFORMATION_MESSAGE);
						}
						
						System.gc();
						Main.chromDraw.repaint();
					}
				      
					Main.glassPane.repaint();				
					Thread.sleep(100);
				
			}
			
			}
			catch(Exception e) {
				ErrorLog.addError(e.getStackTrace());
				e.printStackTrace();
			}
		Main.glassPane.repaint();	
		Draw.setGlasspane(false);
		return null;
	}

}
