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
public class Level implements Comparable<Level> {

	int level, pos, id;
	
	public Level(int level, int pos) {
		this.level = level;
		this.pos = pos;
	
	}
	
	public String toString() {
		return this.pos +"\t" +this.level;
	}
	
	public int compareTo(Level lev) {
		if (this.pos < lev.pos) return -1;
		if (this.pos > lev.pos) return 1;
		if (this.level < lev.level) return -1;
		if (this.level > lev.level) return 1;
	//	if (this.id < lev.id) return -1;
	//	if (this.id > lev.id) return 1;
		return 0;
	}

}
