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
import java.awt.Composite;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.Serializable;
import java.util.ArrayList;

public class SplitClass implements Serializable {
	
	
	private static final long serialVersionUID = 1L;
	int chromEnd = 0, offset, chromOffset;
	public double pixel = 0.0, start=1, end=0, viewLength=0;
	String chrom;	
	//private transient int minReadStart = Integer.MAX_VALUE, maxReadEnd = 0;
//	private transient ArrayList<Transcript> transcripts = new ArrayList<Transcript>();
	private transient ArrayList<Gene> genes = new ArrayList<Gene>();
	//StringBuffer sequence = null;
	//int seqStartPos = 0;	
	boolean updateReads = false, removed = false, clearedReads = true, clearedCoverages = true;	
	//public StringBuffer readSequence = null;
	private transient ReferenceSeq drawReference, readReference;
	private transient Double divider = 4.0;
	public boolean splitRead = false;
	private transient SplitDraw splitDraw; 
	private transient ArrayList<String[]> chromBands = new ArrayList<String[]>();
	int transStart = 0;
	
	public SplitClass() {
		
		splitDraw = new SplitDraw();		
	}		
	void resetSplits() {		
		splitDraw = new SplitDraw();
	}
	
	/*ArrayList<Transcript> getTranscripts() {
		return this.transcripts;
	}
	void setTranscripts(ArrayList<Transcript> trans) {
		this.transcripts = trans;
	}*/
	public void clearGenes() {
		this.genes.clear();
	}
	public void setChromBands(ArrayList<String[]> bands) {
		this.chromBands = bands;
	}
	public ArrayList<String[]> getChromBands() {
		return this.chromBands;
	}
	ReferenceSeq getReference() {
		return this.drawReference;
	}
	void setReference(ReferenceSeq ref) {
		this.drawReference = ref;
	}
	ReferenceSeq getReadReference() {
		return this.readReference;
	}
	void setReadReference(ReferenceSeq ref) {
		this.readReference = ref;
	}
	void nullRef() {
		this.readReference = null;
		this.drawReference = null;
	}
	
	/*int getMinReadStart() {
		return this.minReadStart;
	}
	void setMinReadStart(int value) {
		this.minReadStart = value;
	}
	int getMaxReadEnd() {
		return this.maxReadEnd;
	}
	void setMaxReadEnd(int value) {
		this.maxReadEnd = value;
	}*/
	ArrayList<Gene> getGenes() {
		return this.genes;
	}
	void setGenes(ArrayList<Gene> genes) {
		this.genes = genes;
	}
	void setDivider(Double divider) {
		this.divider = divider;
	}
	SplitDraw getSplitDraw() {
		return this.splitDraw;
	}
	Double getDivider() {
		return this.divider;
	}
	Graphics2D getExonImageBuffer() {
		return splitDraw.exonImageBuffer;
	}
	Graphics2D getReadBuffer() {
		return splitDraw.readBuffer;
	}
	Graphics2D getSelectbuf() {
		return splitDraw.selectbuf;
	}
	BufferedImage getExonImage() {
		if(splitDraw == null) {
			return null;
		}
		else {
			return splitDraw.exonImage;
		}
	}
	BufferedImage getReadImage() {
		if(splitDraw == null) {
			return null;
		}
		return splitDraw.readImage;
	}
	BufferedImage getSelectbuffer() {
		return splitDraw.selectbuffer;
	}
	Composite getBackupe() {
		return splitDraw.backupe;
	}
	Composite getBackupr() {
		return splitDraw.backupr;
	}
	Composite getBackups() {
		return splitDraw.backups;
	}
	BufferedImage getCytoImage() {
		return splitDraw.cytoImage;
	}
	void setCytoImage(BufferedImage image) {
		splitDraw.cytoImage = image;
	}
	void removeSplitDraw() {
		this.splitDraw = null;
	}
	public class SplitDraw {
		Graphics2D exonImageBuffer;
		BufferedImage exonImage;
		Graphics2D readBuffer;
		BufferedImage readImage;
		Graphics2D selectbuf;
		BufferedImage selectbuffer;
		Composite backupe, backupr, backups;
		BufferedImage cytoImage;
		
		public SplitDraw() {
			resizeImages((int)Main.screenSize.getWidth());		
		//	exonImageBuffer.setRenderingHints(Draw.rh);			
		//	readBuffer.setRenderingHints(Draw.rh);			
		//	selectbuf.setRenderingHints(Draw.rh);
			selectbuf.setStroke(Draw.strongStroke);
			
		}
		public void resizeImages(int width) {
		
			exonImage = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
			exonImageBuffer = (Graphics2D)exonImage.getGraphics();
			exonImageBuffer.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB);
			
			readImage = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
			readBuffer = (Graphics2D)readImage.getGraphics();
			readBuffer.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB);
			backupr = readBuffer.getComposite();
			backupe = exonImageBuffer.getComposite();		
			selectbuffer = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
			selectbuf = (Graphics2D)selectbuffer.getGraphics();
			selectbuf.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB);
			backups = selectbuf.getComposite();
			
		}
		
		
	}
}
