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
public class SampleNode {
//	private String variation;

	private final Short quality,coverage, calls;
	private final Sample sample;
	private final ControlFile controlSample;
	private final boolean genotype;	
	Integer alleles, allelenumber;
	
	public SampleNode(int alleles, int allelenumber, ControlFile sample) {
		this.alleles = alleles;
		this.allelenumber = allelenumber;
		this.sample = null;
		this.controlSample = sample;
		this.quality = null;
		this.coverage = null;
		this.calls = null;
		this.genotype = false;
	}
	public SampleNode(short coverage, short calls, boolean genotype, short quality, Sample sample) {
		this.sample = sample;
		this.controlSample = null;
		if(coverage > Short.MAX_VALUE) {
			this.calls = (short)(calls/(double)coverage * Short.MAX_VALUE);
			this.coverage = Short.MAX_VALUE;
			
		}
		else {
			this.coverage = coverage;
			this.calls = calls;
		}
		
		this.quality = quality;
		this.genotype = genotype;
	}
	/*public String getVariation() {
		return " ";
	//	return this.variation;
	}*/
	
	public void addAlleles(int add) {
		this.alleles += add;
	}
	public Short getCoverage() {
		return this.coverage;
	}
	public Short getCalls() {
		return this.calls;
	}
	public Short getQuality() {
		return this.quality;
	}
	public Sample getSample() {
		return this.sample;
	}
	public ControlFile getControlSample() {
		return this.controlSample;
	}
	public Double getAlleleFraction() {
		return this.calls/(double)(this.coverage);
	}
	public boolean isHomozygous() {
		return this.genotype;
	}
	
}
