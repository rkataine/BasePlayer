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
import java.awt.Dimension;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.JButton;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;


public class ErrorLog extends JPanel implements ActionListener {
	
	private static final long serialVersionUID = 1L;
	static JFrame frame = new JFrame("BasePlayer log");    
	static ArrayList<String> errors = new ArrayList<String>();

	JButton save = new JButton("Save");
	
	static JTextArea errorField  = new JTextArea();	
	static boolean found = false;
	
		
	static void addError(StackTraceElement[] elements) {
		found = false;
		for(int i = 0; i<elements.length; i++) {
			if(elements[i].getLineNumber() < 0) {
				continue;
			}
			if(!errors.contains("Error in " +elements[i].getClassName() +":" +elements[i].getMethodName() +" at line " +elements[i].getLineNumber())) {
				errors.add("Error in " +elements[i].getClassName() +":" +elements[i].getMethodName() +" at line " +elements[i].getLineNumber());
				errorField.append(errors.get(errors.size()-1) +"\n");
				found = true;
			}
			
		}
		if(found) {
			errorField.setText(errorField.getText()  +"--------------------------\n");
			if(errors.size() > 100) {
				while(errors.size() > 100) {
					errors.remove(0);
				}
			}			
		}	
	}
	
	static void addError(String line) {
		
	
		errorField.append(line +"\n");			
		errorField.setText(errorField.getText()  +"--------------------------\n");
		
	}
	public static void main(String[] args) {
		
		 javax.swing.SwingUtilities.invokeLater(new Runnable() {
	           public void run() {
	           
	           		createAndShowGUI();
	           	
	           }
	       });
		
	}
	private static void createAndShowGUI() {	
		   
		JFrame.setDefaultLookAndFeelDecorated(true);
		if(Main.userDir == null) {
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
			 frame.setVisible(true); 
		}
		else {
		 frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE); 
		 frame.setVisible(false);  
		}
		 errorField.setLineWrap(true);
		 errorField.setText("Error log:\n");
		 frame.getContentPane().add(new JScrollPane(errorField));
		 frame.setResizable(true);    
		
		 frame.setPreferredSize(new Dimension(400,400));	
	     frame.pack();
	   
	}
	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		
	}
	
}
