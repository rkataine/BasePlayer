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
import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
//import java.io.Serializable;
import java.util.ArrayList;

public class Serializer {
	
	
	public void serialize(File outfile) {
	
		try {
			FileOutputStream fout = new FileOutputStream(outfile);
			ObjectOutputStream oos = new ObjectOutputStream(fout);   
			ArrayList<Sample> sampleTemp = new ArrayList<Sample>();
			for(int i = 0; i<Main.drawCanvas.sampleList.size(); i++) {
				if(!Main.drawCanvas.sampleList.get(i).multipart) {
					sampleTemp.add(Main.drawCanvas.sampleList.get(i));
				}
			}
			oos.writeObject(sampleTemp);
			sampleTemp =null;
			oos.writeObject(Main.drawCanvas.splits);
			Main.drawCanvas.drawVariables.scrollbarpos = Main.drawScroll.getVerticalScrollBar().getValue();
			Main.drawCanvas.drawVariables.projectFile = outfile;
			Main.drawCanvas.drawVariables.projectName = outfile.getName().substring(0, outfile.getName().indexOf(".ses"));
			Main.frame.setTitle("BasePlayer - Project: " +Main.drawCanvas.drawVariables.projectName);
			
			oos.writeObject(Main.drawCanvas.drawVariables);
			oos.writeObject(Control.controlData);
			oos.writeObject(Main.bedCanvas.bedTrack);
			
			oos.writeObject(Settings.settings);
			VariantHandler.saveValues();
			oos.writeObject(VariantHandler.variantSettings);
			
			oos.flush();
			oos.close();
				
			
		}
		catch(Exception e) { 
			e.printStackTrace();
		}
	}
}
