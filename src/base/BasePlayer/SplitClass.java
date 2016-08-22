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
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.Serializable;
import java.util.ArrayList;

public class SplitClass implements Serializable {
	
	
	private static final long serialVersionUID = 1L;
	int chromEnd = 0, offset;
	double pixel = 0.0, start=1, end=0, viewLength=0;
	String chrom;	
//	private transient ArrayList<Transcript> transcripts = new ArrayList<Transcript>();
	private transient ArrayList<Gene> genes = new ArrayList<Gene>();
	//StringBuffer sequence = null;
	//int seqStartPos = 0;	
	boolean updateReads = false, removed = false, clearedReads = true, clearedCoverages = true;	
	//public StringBuffer readSequence = null;
	ReferenceSeq drawReference;
	private transient Double divider = 4.0;
	public boolean splitRead = false;
	private transient SplitDraw splitDraw; 
	int transStart = 0;
	
	public SplitClass() {
		
		splitDraw = new SplitDraw();		
	}		
	void resetSplits() {
		
		splitDraw = new SplitDraw();
		FileRead reader = new FileRead();
	//	setTranscripts(reader.getExons(chrom));
		setGenes(reader.getExons(chrom));
	}
	
	/*ArrayList<Transcript> getTranscripts() {
		return this.transcripts;
	}
	void setTranscripts(ArrayList<Transcript> trans) {
		this.transcripts = trans;
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
		return splitDraw.exonImage;
	}
	BufferedImage getReadImage() {
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
			exonImageBuffer.setRenderingHints(Draw.rh);			
			readBuffer.setRenderingHints(Draw.rh);			
			selectbuf.setRenderingHints(Draw.rh);
			selectbuf.setStroke(Draw.strongStroke);
		}
		public void resizeImages(int width) {
			exonImage = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
			exonImageBuffer = (Graphics2D)exonImage.getGraphics();
			
			readImage = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
			readBuffer = (Graphics2D)readImage.getGraphics();
			
			backupr = readBuffer.getComposite();
			backupe = exonImageBuffer.getComposite();		
			selectbuffer = MethodLibrary.toCompatibleImage(new BufferedImage(width, (int)Main.screenSize.getHeight(), BufferedImage.TYPE_INT_ARGB));	
			selectbuf = (Graphics2D)selectbuffer.getGraphics();
			
			backups = selectbuf.getComposite();
			
		}
		
		
	}
}
