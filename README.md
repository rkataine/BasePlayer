# BasePlayer
BasePlayer for genetic data analysis

Go to https://baseplayer.fi for more information about the software.
Download BasePlayer package with executable jar file or operating system specific installers from https://baseplayer.fi/downloads.html.

Build instructions for developers:

<b>Build instructions using Eclipse Oxygen IDE</b> (https://www.eclipse.org/oxygen/):

1. Clone or download/unzip the BasePlayer git (project https://github.com/rkataine/BasePlayer.git)
to your computer to any folder (in this example we use /BasePlayer/).

2. Open Eclipse and create and launch new Workspace as prompted.

3. Create new Java project (File -> New -> Java Project).
 
4. Set project name e.g. BasePlayer and uncheck "Use default location". Set your BasePlayer folder (/BasePlayer/ in this case) as a project location. Press "Finish".

5. In "Package Explorer", expand BasePlayer -> src -> base.BasePlayer and double click Main.java.

6. Go to Run -> Run as -> Java Application (BasePlayer should start).

7. To create executable JAR file, go to File -> Export -> Java, select "Runnable JAR file" and press Next.

8. Use "Main - BasePlayer" as a Launch configuration and your BasePlayer folder (/BasePlayer/ in this case) as an Export destination.

9. Press "Finish" and executable JAR file will be created to your BasePlayer folder.

10. Double-click BasePlayer.jar to start the program or type java -jar BasePlayer.jar in the /baseplayer/ folder

<b>Build instructions for command line (linux):</b>
 
1. Download and install Java SE Development kit 8 or newer at http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
and add it to your environment path


2. Clone or download/unzip the BasePlayer git (project https://github.com/rkataine/BasePlayer.git)
to your computer to any folder (in this example we use /BasePlayer/).

3. Goto /BasePlayer/src directory and unpack all JAR files from /BasePlayer/jars to this folder:
	jar -xvf ../jars/commons-compress-1.17.jar
	jar -xvf ../jars/commons-io-2.4.jar
	jar -xvf ../jars/htsjdk.jar
	jar -xvf ../jars/WigReader.jar
	jar -xvf ../jars/commons-net-3.5.jar

4. Stay in /BasePlayer/src folder and compile java-files and copy external JAR folders using the following commands (you can ignore possible warnings or notes):
	javac -cp ./base/BasePlayer/*.java -d ../build -classpath .
 	cp -R ./base/BBfile/ ../build/base/
	mv htsjdk/ org/ ../build/

5. Go back to /BasePlayer/ folder and copy resource folders to build folder:
	cp -R ./src/base/BasePlayer/SELEX/ ./src/base/BasePlayer/icons/ ./build/base/BasePlayer/
	
6. Create executable BasePlayer.jar by using following command in /BasePlayer/ folder:
	jar cvfm BasePlayer.jar ./build/META-INF/MANIFEST.MF -C ./build .
	
7. Double-click BasePlayer.jar to start the program or type java -jar BasePlayer.jar in the /baseplayer/ folder
