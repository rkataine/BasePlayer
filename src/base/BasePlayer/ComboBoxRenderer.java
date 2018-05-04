package base.BasePlayer;

import java.awt.Component;

import javax.swing.DefaultListCellRenderer;
import javax.swing.JList;

public class ComboBoxRenderer extends DefaultListCellRenderer {

   
	private static final long serialVersionUID = 1L;

	public Component getListCellRendererComponent(
                                   JList list,
                                   Object value,
                                   int index,
                                   boolean isSelected,
                                   boolean cellHasFocus) {
		if (value instanceof Sample) {        	
        	value = ((Sample)value).getName();
        }
		else {
			value = "-";
		}
        
        super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
        return this;
    }
}