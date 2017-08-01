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
import java.util.HashMap;

public class SampleNode {
//	private String variation;

	private static final long serialVersionUID = 1L;
	private final Integer coverage, calls;
	private final Float quality, gq;
	float heightValue = 0;
	private final Sample sample;
	private final ControlFile controlSample;
	private final boolean genotype;	
	Integer alleles, allelenumber;
	//HashMap<String, Float> advQualities;
	
	public SampleNode(int alleles, int allelenumber, ControlFile sample) {
		this.alleles = alleles;
		this.allelenumber = allelenumber;
		this.sample = null;
		this.controlSample = sample;
		this.quality = null;
		this.gq = null;
		this.coverage = null;
		this.calls = null;
		this.genotype = false;
	}
	public SampleNode(int coverage, int calls, boolean genotype, Float quality, Float gq, HashMap<String, Float> advquals, Sample sample) {
		this.sample = sample;
		//this.advQualities = advquals;
		this.controlSample = null;
		/*if(coverage > Short.MAX_VALUE) {
			this.calls = (Float)(calls/(double)coverage * Short.MAX_VALUE);
			this.coverage = Short.MAX_VALUE;
			
		}
		else {*/
			this.coverage = coverage;
			this.calls = calls;
	//	}
		this.gq = gq;
		this.quality = quality;
		this.genotype = genotype;
		switch(Settings.selectedVarDraw) {	
		
			case 0: {	
				if(sample.getMaxCoverage() < coverage) {
					
					sample.setMaxCoverage((float)(coverage));
					
				}	
				this.heightValue = coverage;
				break;
			}
			case 1: {	
				if(sample.getMaxCoverage() < calls/(float)coverage) {
					
					sample.setMaxCoverage(calls/(float)coverage);
					
				}	
				this.heightValue = calls/(float)coverage;
				break;
			}
			case 2: {	
				if(quality != null) {
					if(sample.getMaxCoverage() < quality) {
						
						sample.setMaxCoverage(quality);
						if(VariantHandler.maxSlideValue < quality) {
							VariantHandler.maxSlideValue = quality;
						}
					}	
					this.heightValue = quality;
				}
				break;
			}
			case 3: {	
				if(gq != null) {
					if(sample.getMaxCoverage() < gq) {
						
						sample.setMaxCoverage(gq);
						if(VariantHandler.maxSlideValue < gq) {
							VariantHandler.maxSlideValue = gq;
						}
					}	
					this.heightValue = gq;
				}
				
				break;
			}
			case 4: {	
				if(sample.getMaxCoverage() < calls) {
					
					sample.setMaxCoverage((float)(calls));
					
				}	
				this.heightValue = calls;
				break;
			}
		}
	}
	/*public String getVariation() {
		return " ";
	//	return this.variation;
	}*/
	
	public void addAlleles(int add) {
		this.alleles += add;
	}
	public int getCoverage() {
		return this.coverage;
	}
	public int getCalls() {
		return this.calls;
	}
	public Float getQuality() {
		
		return this.quality;
	}
	public Float getGQ() {		
		return this.gq;
	}
	public String getGQString() {
		if(this.gq == null) {
			return ".";
		}
		else {
			return ""+this.gq;
		}
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
