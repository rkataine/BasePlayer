package base.BasePlayer;


import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class ProxySettings extends JPanel {
	
	private static final long serialVersionUID = 1L;
	static JCheckBox useProxy;
	static String host, port, username, password;
	static JTextField hostField, portField, userField, passField;
	static String[] types = {"HTTP", "SOCKS", "DIRECT"};
	static JComboBox<String> proxytypes;	
	static JButton save;	
	
	public ProxySettings() {
		super(new GridBagLayout());
		useProxy = new JCheckBox("Use proxy");
		proxytypes = new JComboBox<String>(types);
		hostField = new JTextField("Proxy host");
		portField = new JTextField("Port");
		save = new JButton("Save");
		//userField = new JTextField("Username");
		//passField = new JTextField("Pass");
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.NORTHWEST;		
		c.insets = new Insets(2,0,0,30);
		c.gridx = 0;
		c.weightx = 1;
		
		c.gridy = 0;
		c.gridwidth = 1;	
		add(useProxy,c);
		c.gridy++;
		add(proxytypes,c);
		c.gridy++;
		add(hostField,c);
		c.gridy++;
		add(portField,c);
		c.gridy++;
		c.fill = GridBagConstraints.NONE;
		
		add(save,c);
		save.setPreferredSize(Main.buttonDimension);
		c.gridy++;
		c.weighty = 1;
		add(new JLabel(),c);
		//add(userField);
		//add(passField);
		save.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				if(useProxy.isSelected()) {
					Main.writeToConfig("isProxy=true");
				}
				else {
					Main.writeToConfig("isProxy=false");
				}
				Main.writeToConfig("proxyHost=" +hostField.getText());
				Main.writeToConfig("proxyPort=" +portField.getText());
				Main.writeToConfig("proxyType=" +proxytypes.getSelectedItem().toString());
				
			}
			
		});
	}

	
	
}
