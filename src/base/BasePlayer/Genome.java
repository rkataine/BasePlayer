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

import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.ArrayList;

public class Genome implements Serializable {
	
	private static final long serialVersionUID = 1L;
	RandomAccessFile referenceFile;
	ArrayList<String> annotations = new ArrayList<String>();
	String referenceName = null;
	String cytobands = null;	
	int defaultAnnotation = 0;
	
	public Genome(String ref) {
		try {
			this.referenceFile = new RandomAccessFile(ref,"r");
		}
		catch(Exception e) {
			e.printStackTrace();
			Main.showError(e.getMessage(), "Error");
		}
	}
	public void removeAnnotation(String annotation) {
		this.annotations.remove(annotation);
		if(this.defaultAnnotation > this.annotations.size()-1) {
			this.defaultAnnotation = 0;
		}
	}
	public void addAnnotation(String annotation) {
		if(!annotations.contains(annotation)) {
			annotations.add(annotation);
		}		
		this.defaultAnnotation = this.annotations.indexOf(annotation);
	}
	public void setCyto(String cyto) {
		this.cytobands = cyto;
	}
	public void setReferenceFile(String ref) {
		try {
			this.referenceFile = new RandomAccessFile(ref,"r");
		}
		catch(Exception e) {
			e.printStackTrace();
			Main.showError(e.getMessage(), "Error");
		}
	}
	public void setName(String name) {
		this.referenceName = name;
	}
}
