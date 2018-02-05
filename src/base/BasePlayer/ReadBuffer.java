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
import htsjdk.samtools.SAMRecord;

public class ReadBuffer {
	SAMRecord record;
	ReadBuffer prev, next;
	
	public ReadBuffer(SAMRecord record, ReadBuffer prev) {
		this.record = record;
		this.prev = prev;
		
	}
}
