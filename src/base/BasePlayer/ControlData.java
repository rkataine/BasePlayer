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
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;


public class ControlData implements Serializable{

	
	
	private static final long serialVersionUID = 1L;
//	java.util.List<Object[]> data = Collections.synchronizedList(new ArrayList<Object[]>());   
	int sampleCount = 0;
//	double alleleThreshold = 0.01;
//	java.util.List<String> sampleArray = Collections.synchronizedList(new ArrayList<String>());	
	java.util.List<ControlFile> fileArray = Collections.synchronizedList(new ArrayList<ControlFile>());
	int total = 0;
	boolean controlsOn = false;
	
	
	
	
}
