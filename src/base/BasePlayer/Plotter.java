package base.BasePlayer;

import java.awt.Color;
import java.awt.Graphics;

import javax.swing.JPanel;

public class Plotter extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	int[][] array = null;
	int width, barwidth;
	Color lightgray = new Color(220,220,220);
	int samplewidth = 0, scale = 0, valuecount = 0;
	int bottom = 310, highest = 0;
	public Plotter(int[][] array, int width) {
		this.width = width;
		samplewidth = this.width; ///Main.samples;
		this.valuecount = array[0].length;
		barwidth = samplewidth/valuecount;
		this.array = array;
		int sum = 0, countsum = 0;
		for(int i = 0; i<array.length; i++) {			
			sum = 0;
			countsum = 0;
			for(int j = 0; j<array[0].length; j++) {
				if(array[i][j] > highest) {
					highest = array[i][j];					
				}
				sum+=array[i][j]*(j+1);
				countsum += (j+1);
				System.out.print(array[i][j] +" ");
			}
			System.out.println();
			System.out.println(sum/(double)countsum);
		}		
	}
	public void paint(Graphics g) {
		try {					
			
			if(array != null) {
				for(int i = 0; i<array.length; i++) {
					if(i%2 == 0) {
						g.setColor(lightgray);
						g.fillRect(samplewidth*i, 10, samplewidth, bottom);
					}
				}
				for(int i = 0;i<11;i++) {
					g.setColor(Color.LIGHT_GRAY);
					g.fillRect(0,bottom-i*30,width+30,1);
					g.setColor(Color.black);
					g.fillRect(width+20,bottom-i*30,10,1);					
					g.drawString(""+MethodLibrary.map(i*30, 0,300,0,highest), width,bottom-i*30+10 );
				}
				/*for(int i = 0; i<array.length; i++) {					
					g.setColor(Color.black);
					for(int j = 0; j<array[0].length; j++) {
						g.fillRect(samplewidth*i+(j-VariantHandler.callSlider.getValue())*barwidth, bottom-MethodLibrary.map(array[i][j],0,highest,0,300), barwidth, MethodLibrary.map(array[i][j],0,highest,0,300));
						g.drawString(""+j,samplewidth*i+(j-VariantHandler.callSlider.getValue())*barwidth+2,bottom+10);
					
					}
					g.drawString(Main.drawCanvas.sampleList.get(i).getName(),samplewidth*i, bottom+24);
			
				}*/
								
					g.setColor(Color.black);
					for(int j = 0; j<array[0].length; j++) {
						g.fillRect((j-VariantHandler.callSlider.getValue())*barwidth, bottom-MethodLibrary.map(array[Main.drawCanvas.selectedIndex][j],0,highest,0,300), barwidth, MethodLibrary.map(array[Main.drawCanvas.selectedIndex][j],0,highest,0,300));
						g.drawString(""+j,(j-VariantHandler.callSlider.getValue())*barwidth+2,bottom+10);
					
					}
					//g.drawString(Main.drawCanvas.sampleList.get(i).getName(),samplewidth*i, bottom+24);
			
				
				
			}
			
		}
		catch(Exception e) {
			ErrorLog.addError(e.getStackTrace());
			e.printStackTrace();
		}
	}
	
	
	
}
