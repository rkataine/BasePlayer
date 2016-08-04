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

public class AminoEntry {

	private VarNode node;
	private String[] row;
	
	public AminoEntry(String[] row, VarNode node) {
		this.node = node;
		this.row = row;
	}
	
	public VarNode getNode() {
		return this.node;
	}
	public String[] getRow() {
		return this.row;
	}
	
}
