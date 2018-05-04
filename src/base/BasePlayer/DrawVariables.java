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
import java.io.Serializable;
import java.util.ArrayList;


public class DrawVariables implements Serializable {	
	private static final long serialVersionUID = 1L;
	double sampleHeight = Draw.defaultSampleHeight;
	Boolean somatic = false; 
	short visiblesamples=1,visiblestart=0;
	int scrollbarpos = 0;
	String projectName = "Untitled";
	File projectFile;
	
	ArrayList<QualEntry> advQDraw = null;
	ArrayList<QualEntry> advQDrawIndel = null;
	
}