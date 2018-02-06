package base.BasePlayer;

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

import htsjdk.samtools.SAMRecord;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Iterator;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingWorker;

import base.BBfile.BigWigIterator;
import base.BBfile.WigItem;

public class PeakCaller extends JPanel implements ActionListener, ComponentListener {
	
	private static final long serialVersionUID = 1L;
	static Boolean loading = false;
	static JFrame frame = new JFrame("Peak Caller");
	JLabel testing = new JLabel("Read Peak Caller test");
	static JTextField minreads = new JTextField("10");
	JLabel minreadlabel = new JLabel("Minimum coverage for calling:");
	static JTextField minwidth = new JTextField("10");
	JLabel minwidthlabel = new JLabel("Minimum width of region:");
	static String savepath = "";
	static JLabel info = new JLabel("");
	//JCheckBox onlySel = new JCheckBox("Calc only for selected");
	static JCheckBox allchroms = new JCheckBox("All chromosomes");
	JButton execute = new JButton("Execute");
	public PeakCaller(boolean caller) {
		
	}
	public PeakCaller() {
		super(new GridLayout(8,1));		
		try {
			frame.addComponentListener(this);
			boolean bams = false;
			/*for(int i = 0 ;i<Main.samples; i++) {
				if(Main.drawCanvas.sampleList.get(i).samFile != null) {
					bams = true;
					break;
				}
			}*/
			for(int i = 0 ;i<Main.bedCanvas.bedTrack.size(); i++) {
				if(Main.bedCanvas.bedTrack.get(i).file.getName().toLowerCase().endsWith(".bw") || Main.bedCanvas.bedTrack.get(i).file.getName().toLowerCase().endsWith(".bigwig")) {
					bams = true;
					break;
				}
			}
			if(!bams) {
				//info.setText("Open BAM/CRAM files to call the read peaks.");
				info.setText("Open BigWig files to call the read peaks.");
				execute.setEnabled(false);
			}
			
			PeakCaller.savepath = savepath +"/";
			testing.setForeground(Color.red);
			execute.addActionListener(this);
			allchroms.addActionListener(this);
		
			
			add(testing);			
			add(info);			
			add(minreadlabel);
			add(minreads);
			add(minwidthlabel);
			add(minwidth);
			add(allchroms);		
			add(execute);		
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void createAndShowGUI() {	
	
			if(Main.frame != null) {
				frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE); 
				frame.setVisible(false);		
			}
			else {
				frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
				frame.setVisible(true);		
			}
			frame.setAlwaysOnTop(true);
		 	
		    JComponent newContentPane = new PeakCaller();
		    newContentPane.setOpaque(false);		   
		    frame.setContentPane(newContentPane);
		    setFonts(Main.menuFont);
		    frame.pack();
		   
		    for(int i =0; i<frame.getContentPane().getComponentCount(); i++) {
		    	frame.getContentPane().getComponent(i).setMinimumSize(frame.getContentPane().getComponent(i).getPreferredSize());
		    }
		    minreads.setMinimumSize(new Dimension(minreads.getWidth(),minreads.getHeight()));	
		    minwidth.setMinimumSize(new Dimension(minwidth.getWidth(),minwidth.getHeight()));	
			
	}
	public static void main(String[] args) {
		
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	          public void run() {	          
	          		createAndShowGUI();
	          }
	      });	
	}
	
	static void setFonts(Font menuFont) {
		for(int i = 0 ; i<PeakCaller.frame.getContentPane().getComponentCount(); i++) {
			PeakCaller.frame.getContentPane().getComponent(i).setFont(menuFont);
		}
		PeakCaller.frame.pack();
	}
	@Override
	public void actionPerformed(ActionEvent event) {		
		if(event.getSource() == execute) {
			try {				
				 JFileChooser chooser = new JFileChooser();
		    	  chooser.setAcceptAllFileFilterUsed(false);		
		    	  chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		    	  chooser.setDialogTitle("Save peak results in...");
		    	  
		          int returnVal = chooser.showSaveDialog((Component)this.getParent());	           	  
			       
			      if(returnVal == JFileChooser.APPROVE_OPTION) {   			
			    	  
			    	  PeakRunner caller = new PeakRunner(chooser.getSelectedFile().getCanonicalPath());
			    	  caller.execute();											
			      }
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
	
	}

public class PeakRunner extends SwingWorker<String, Object> {
		int minwidth, minquality, minreadquality, mincoverage; 
		double readlevel;
		boolean onlysel;
		String savepath = "";
		public PeakRunner(int readlevel, int minquality, int minreadquality, int mincoverage, boolean onlysel) {
			this.readlevel = readlevel;
			this.minquality = minquality;
			this.minreadquality = minreadquality;
			this.mincoverage = mincoverage;
			this.onlysel = onlysel;
		}
		public PeakRunner(String savepath) {
			this.savepath = savepath +"/";
			int minwidth; 
			double readlevel;
			try {
				readlevel = Double.parseDouble(PeakCaller.minreads.getText());
				minwidth = Integer.parseInt(PeakCaller.minwidth.getText());
			}
			catch(Exception e) {
				PeakCaller.minreads.setForeground(Color.red);
				return;
			}			
		
			this.readlevel = readlevel;			
			this.minwidth = minwidth;
			//this.onlysel = onlySel.isSelected();
		}
		void callWigPeaks() {
			
			PeakCaller.loading = true;
			
			if(!Main.drawCanvas.loading) {
				Main.drawCanvas.loading("Calling peaks");
			}
			else {
				Main.drawCanvas.loadingtext = "Calling peaks";
			}
			int chromindex = 0;
			ArrayList<BedTrack> samples = new ArrayList<BedTrack>();
			ArrayList<BufferedWriter> samplewriter = new ArrayList<BufferedWriter>();
			String samplename;
			try {
				for(int i = 0 ; i<Main.bedCanvas.bedTrack.size(); i++) {
					if(Main.bedCanvas.bedTrack.get(i).file != null) {
						
						samples.add(Main.bedCanvas.bedTrack.get(i));
						samplename = Main.bedCanvas.bedTrack.get(i).file.getName().toLowerCase().replace(".bw", "").replace(".bigwig",  "") +"_peaks.bed";						
						BufferedWriter writer = new BufferedWriter(new FileWriter(savepath +samplename));
						samplewriter.add(writer);
						writer.write("#chrom\tstart\tend\tOverlappingGenes\tHighestPeak\tRegionWidth\n");						
					}
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			
			if(allchroms.isSelected()) {
				Main.nothread = true;
				Main.chromosomeDropdown.setSelectedIndex(0);
			}
			BigWigIterator wigiter = null;
			WigItem features = null;
			int loadbarsamplevalue = 0;
			while(chromindex < 25) {
				if(!Main.drawCanvas.loading) {
					break;
				}
				for(int s = 0; s<samples.size();s++) {	
					try {	
						if(!Main.drawCanvas.loading) {
							break;
						}
						loadbarsamplevalue = (int)(s/(double)samples.size()*100);  
						Main.drawCanvas.loadBarSample = loadbarsamplevalue;						
						BedTrack track = samples.get(s);															
							
						int startpos = (int)Main.drawCanvas.splits.get(0).start;
						int endpos = (int)Main.drawCanvas.splits.get(0).end;						
						int coveragestart = 0, coverageend = 0, prevend = 0;
												
						boolean start = false;
						wigiter = track.getBBfileReader().getBigWigIterator(track.chr+Main.chromosomeDropdown.getSelectedItem().toString() , startpos, track.chr+Main.chromosomeDropdown.getSelectedItem().toString(), endpos, false);
						track.setZoomlevel(1);
						int regionstart = 0, regionend = 0;
						Float maxvalue = 0F;
						while(wigiter.hasNext()) {
							   if(!Main.drawCanvas.loading) {							
								   break;
							   }
						    	try {								
						    		features = wigiter.next();
						    	}
						    	catch(Exception ex) {
						    		ex.printStackTrace();
						    	}
						    	if(features.getWigValue() >= readlevel && !start) {
						    		start = true;
						    		regionstart = features.getStartBase();
						    		maxvalue = features.getWigValue();
						    		continue;
						    	}
						    	if(features.getWigValue() < readlevel && start) {
						    		regionend = features.getEndBase();
						    		String genes = MethodLibrary.getOverlappingGenes(regionstart, regionend, Main.drawCanvas.splits.get(0));
									Main.drawCanvas.loadbarAll  = (int)(regionstart/(double)Main.drawCanvas.splits.get(0).chromEnd*100);
									Main.drawCanvas.loadBarSample = loadbarsamplevalue + ((int)(regionstart/(double)Main.drawCanvas.splits.get(0).chromEnd*100)/samples.size());
									
						    		if(regionend-regionstart >= minwidth) {
						    			try {
						    				
						    				samplewriter.get(s).write(Main.chromosomeDropdown.getSelectedItem() +"\t" +(regionstart+1) +"\t" +regionend +"\t" +genes +"\t" +maxvalue +"\t" +(regionend-regionstart) +"\n");
						    			}
						    			catch(Exception e) {
						    				e.printStackTrace();
						    			}
						    			//System.out.println(features.getChromosome().replace("chr", "") +"\t" +(regionstart+1) +"\t" +regionend +"\t" +genes +"\t" +maxvalue +"\t" +(regionend-regionstart));
						    		}
						    		
							    	start = false;
						    		continue;
						    	}
						    	if(maxvalue < features.getWigValue()) {
						    		maxvalue = features.getWigValue();
						    	}
						    	
						}
															
					}
					catch(Exception e) {
						e.printStackTrace();
						break;
					}						
				}			
				if(!allchroms.isSelected()) {
					break;
				}
				chromindex++;
				if(chromindex == 25) {
					break;
				}
				if(allchroms.isSelected()) {
					Main.nothread = true;
					Main.chromosomeDropdown.setSelectedIndex(chromindex);
				}
			
			}
			try {
				for(int i = 0 ; i<samplewriter.size(); i++) {
					samplewriter.get(i).close();					
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			Main.drawCanvas.ready("Calling peaks");			
					
			
			Main.drawCanvas.repaint();
			PeakCaller.loading = false;
			
		}
		void callPeaks() {

			if(onlysel && (Main.drawCanvas.selectedSample == null || Main.drawCanvas.selectedSample.getTabixFile() != null || Main.drawCanvas.selectedSample.samFile == null)) {
				PeakCaller.loading = false;
				return;			
			}		
			
			FileRead reader = new FileRead();
			
			PeakCaller.loading = true;
			
			if(!Main.drawCanvas.loading) {
				Main.drawCanvas.loading("Calling peaks");
			}
			else {
				Main.drawCanvas.loadingtext = "Calling peaks";
			}
			int chromindex = 0;
			ArrayList<Sample> samples = new ArrayList<Sample>();
			ArrayList<BufferedWriter> samplewriter = new ArrayList<BufferedWriter>();
			String samplename;
			try {
				for(int i = 0 ; i<Main.samples; i++) {
					if(Main.drawCanvas.sampleList.get(i).samFile != null) {
						samples.add(Main.drawCanvas.sampleList.get(i));
						samplename = Main.drawCanvas.sampleList.get(i).getName().replace(".bam", "") +"_peaks.bed";
						
						BufferedWriter writer = new BufferedWriter(new FileWriter(savepath +samplename));
						samplewriter.add(writer);
						writer.write("#chrom\tstart\tend\tOverlappingGenes\tHighestPeak\n");
						
					}
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			
			if(allchroms.isSelected()) {
				Main.nothread = true;
				Main.chromosomeDropdown.setSelectedIndex(0);
			}
			int loadbarsamplevalue = 0;
			while(chromindex < 25) {
				if(!Main.drawCanvas.loading) {
					break;
				}
				for(int s = 0; s<samples.size();s++) {	
					try {	
						
						if(!Main.drawCanvas.loading) {
							break;
						}
						loadbarsamplevalue = (int)(s/(double)Main.samples*100);  
						Main.drawCanvas.loadBarSample = loadbarsamplevalue;
						
						Sample sample = samples.get(s);
						sample.basequalsum = 0;
						sample.basequals = 0;
						if(onlysel && !sample.equals(Main.drawCanvas.selectedSample)) {
							continue;
						}
						
						Reads reads = null;
						try {
							reads = (Reads)sample.getreadHash().get(Main.drawCanvas.splits.get(0)).clone();
							
						}
						catch(Exception e) {
							e.printStackTrace();
						}				
																	
							if(!Main.drawCanvas.loading) {
								Main.drawCanvas.ready("all");
								
								break;
							}
							int startpos = (int)Main.drawCanvas.splits.get(0).start;
							int endpos = (int)Main.drawCanvas.splits.get(0).end;
							
							Iterator<SAMRecord> bamIterator = reader.getBamIterator(reads,reads.sample.chr + Main.drawCanvas.splits.get(0).chrom,startpos, endpos);	
							SAMRecord samRecord = null;
							int coveragestart = 0, coverageend = 0, prevend = 0;
							ArrayList<Integer> coverages = new ArrayList<Integer>();
							for(int i = 0 ; i<2000; i++) {
								coverages.add(0);
							}
							boolean start = false, printed = false;
							
							int regionstart = 0, regionend = 0, maxvalue = 0;
							while(bamIterator != null && bamIterator.hasNext()) {	
								if(!Main.drawCanvas.loading) {
									Main.drawCanvas.ready("all");
									break;
								}
								
								try {
									samRecord = bamIterator.next(); 	
									
								}
								catch(htsjdk.samtools.SAMFormatException ex) {
									ex.printStackTrace();		
									continue;
								}
								/*if(samRecord.getMappingQuality() < minReadQuality) {
									continue;
								}*/
								if(samRecord.getReadUnmappedFlag()) {					
									continue;
								}						
								
								if(samRecord.getUnclippedEnd() < startpos) { //this.readSeqStart+1) {
									
									continue;
								}
								
								if(samRecord.getUnclippedStart() >= endpos) {					
									break;
								}		
								if(coveragestart != 0) {
									printed = false;
									maxvalue = 0;
									if(prevend < samRecord.getAlignmentStart()) {
										for(int i = 0;i<coverages.size(); i++) {										
											if(coverages.get(i) == 0) {
												break;
											}
											if(start && coverages.get(i) < readlevel) {
												start = false;
												regionend = coveragestart + i;
												String genes = MethodLibrary.getOverlappingGenes(regionstart, regionend, Main.drawCanvas.splits.get(0));
												Main.drawCanvas.loadbarAll  = (int)(regionstart/(double)Main.drawCanvas.splits.get(0).chromEnd*100);
												Main.drawCanvas.loadBarSample = loadbarsamplevalue + ((int)(regionstart/(double)Main.drawCanvas.splits.get(0).chromEnd*100)/samples.size());
												samplewriter.get(s).write(Main.drawCanvas.splits.get(0).chrom +"\t" +regionstart +"\t" +regionend +"\t" +genes +"\t" +maxvalue +"\n");
												//System.out.println(Main.drawCanvas.splits.get(0).chrom +"\t" +regionstart +"\t" +regionend +"\t" +genes +"\t" +maxvalue);
											}
											if(start == false && coverages.get(i) >= readlevel) {
												start = true;
												regionstart = coveragestart + i;
											}
											if(start) {
												if (maxvalue < coverages.get(i)) {
													maxvalue = coverages.get(i);
												}
											//	System.out.print(coverages.get(i) +" ");
											}
											
											coverages.set(i, 0);
										}
										
										
										coverageend = 0;
									}								
									
								}
								if(coverageend == 0) {
									coveragestart = samRecord.getAlignmentStart();
									
								}
								coverageend = samRecord.getAlignmentEnd();
								prevend = coverageend;
								
								for(int i = 0; i<samRecord.getReadLength(); i++) {
									if(samRecord.getAlignmentStart()-coveragestart+i > coverages.size()-1) {
										for(int j = 0 ; j<2000; j++) {
											coverages.add(0);
										}
									}
									coverages.set(samRecord.getAlignmentStart()-coveragestart+i, coverages.get(samRecord.getAlignmentStart()-coveragestart+i)+1);
								}
								
							}
						reads = null;					
					}
					catch(Exception e) {
						e.printStackTrace();
						break;
					}						
				}			
				if(!allchroms.isSelected()) {
					break;
				}
				chromindex++;
				if(chromindex == 25) {
					break;
				}
				if(allchroms.isSelected()) {
					Main.nothread = true;
					Main.chromosomeDropdown.setSelectedIndex(chromindex);
				}
			
			}
			try {
				for(int i = 0 ; i<samplewriter.size(); i++) {
					samplewriter.get(i).close();					
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			Main.drawCanvas.ready("Calling peaks");			
					
			
			Main.drawCanvas.repaint();
			PeakCaller.loading = false;
			
		}
		@Override
		protected String doInBackground() throws Exception {
			callWigPeaks();
			return null;
		}
	}

@Override
public void componentResized(ComponentEvent e) {
	// TODO Auto-generated method stub
	boolean bams = false;
	/*for(int i = 0 ;i<Main.samples; i++) {
		if(Main.drawCanvas.sampleList.get(i).samFile != null) {
			bams = true;
			break;
		}
	}*/
	for(int i = 0 ;i<Main.bedCanvas.bedTrack.size(); i++) {
		if(Main.bedCanvas.bedTrack.get(i).file.getName().toLowerCase().endsWith(".bw") || Main.bedCanvas.bedTrack.get(i).file.getName().toLowerCase().endsWith(".bigwig")) {
			bams = true;
			break;
		}
	}
	
	if(!bams) {
		//info.setText("Open BAM/CRAM files to call the read peaks.");
		info.setText("Open BigWig files to call the read peaks.");
		execute.setEnabled(false);
	}
	else {
		info.setText("");
		execute.setEnabled(true);
	}
}
@Override
public void componentMoved(ComponentEvent e) {
	// TODO Auto-generated method stub
	
}
@Override
public void componentShown(ComponentEvent e) {
	// TODO Auto-generated method stub
	
}
@Override
public void componentHidden(ComponentEvent e) {
	// TODO Auto-generated method stub
	
}
}

