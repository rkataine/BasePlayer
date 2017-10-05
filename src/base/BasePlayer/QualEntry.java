package base.BasePlayer;

import java.io.Serializable;

public class QualEntry implements Serializable {

	private static final long serialVersionUID = 1L;
	String key;
	float value;
	String format;
	public QualEntry(String key, Float value, String format) {
		this.key = key;
		this.value = value;
		this.format = format;
	}
}