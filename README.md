# BasePlayer
BasePlayer for genetic data analysis

Go to https://baseplayer.fi for more information about the software.
Download BasePlayer package with executable jar file or operating system specific installers from https://baseplayer.fi/downloads.html.

Build instructions for developers:

Build instructions using Eclipse Oxygen IDE (https://www.eclipse.org/oxygen/):

1. Clone or download/unzip the BasePlayer git (project https://github.com/rkataine/BasePlayer.git)
to your computer to any folder (in this example we use /BasePlayer/).

2. Open Eclipse and create and launch new Workspace as prompted.

3. Create new Java project (File -> New -> Java Project).
 
4. Set project name e.g. BasePlayer and uncheck "Use default location". Set your BasePlayer folder (/BasePlayer/ in this case) as a project location. Press "Finish".

5. In "Package Explorer", expand BasePlayer -> src -> base.BasePlayer and double click Main.java.

6. Go to Run -> Run as -> Java Application and BasePlayer should start.

7. To create executable JAR file, go to File -> Export -> Java, select "Runnable JAR file" and press Next.

8. Use "Main - BasePlayer" as a Launch configuration and your BasePlayer folder (/BasePlayer/ in this case) as an Export destination.

9. Press "Finish" and executable JAR file will be created to your BasePlayer folder.
