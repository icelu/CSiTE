/*****************************************************************************
 * Copyright (c) 2017 Minh Duc Cao. All rights reserved.                     *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * File: Deploy.java
 * 15/11/2013 - Minh Duc Cao: Created
 * 01/02/2018 - Bingxin Lu: Revised
 *
 ****************************************************************************/

package japsa.util.deploy;


import japsa.sim.SimulateCaptureCmd;

import japsa.util.CommandLine;
import japsa.util.StringSeparator;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Date;
import java.util.Scanner;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;


/**
 * This class is used to deploy tools: create a makefile to generate scripts
 *
 * @author Minh Duc Cao (http://www.caominhduc.org/)
 */
public class Deploy {
	public static ArrayList<Object> tools = new ArrayList<Object>();
	public static String VERSION = "1.7-05b";
	public static final String FIELD_SEP = "\t";

	public static boolean DEBUG = true;

	//private static String AUTHORS = "Minh Duc Cao";
	static{
			//jsa.sim
		tools.add("Alignment with Finite State Machines");
		tools.add(new SimulateCaptureCmd());
	}

	/**
	 * Method to set up Japsa directory
	 * @param japsaJar
	 * @return
	 * @throws IOException
	 */
	static private String setupJapsaDir(String japsaJar) throws IOException{
		if (japsaPath.startsWith("~/")) {
			japsaPath = System.getProperty("user.home") + japsaPath.substring(1);
		}

		File japsaLib = new File(japsaPath + File.separator + "lib" + File.separator + "japsa");
		File japsaBin = new File(japsaPath + File.separator + "bin");

		japsaBin.mkdirs();
		if (!japsaBin.exists()){
			System.err.println("Folder `" +japsaBin.getCanonicalPath() + "' cannot be created");
			return null;
		}

		japsaLib.mkdirs();
		if (!japsaLib.exists()){
			System.err.println("Folder `" +japsaLib.getCanonicalPath() + "' cannot be created");
			return null;
		}

		//Get the classpath now that the directory of installation is known
        File from = new File(japsaJar);

		File to = new File (japsaLib.getCanonicalPath()  + File.separator + japsaJar);
		try{
			Files.copy(from.toPath(),to.toPath(),REPLACE_EXISTING);
		}catch (IOException e){
			System.err.println(e.getMessage());
			return null;
		}

		String cp = to.getCanonicalPath();
		StringSeparator ss = new StringSeparator(libs, ':');
		while (ss.hasNext()) {
			String l = ss.next();
			if (l.length() > 0){
				from = new File("libs" + File.separator + l);
				to = new File (japsaLib.getCanonicalPath()  + File.separator + l);
				try{
                    Files.copy(from.toPath(),to.toPath(),REPLACE_EXISTING);
					//Files.copy(from,to);
				}catch (IOException e){
					System.err.println(e.getMessage());
					return null;
				}
				cp = cp + File.pathSeparator + to.getCanonicalPath();
			}
		}
		System.out.println("  Directory " + japsaLib.getCanonicalPath() + " set up");
		return cp;
	}

	public static int server = -1;

	public static String maxMem = null;
	public static String japsaPath = null;
	// public static String jlp = null;
	public static String libs = null;
	public static String japsaJar = null;
	public static String compiler = null;

	private static String classPath = null;
	private static String javaCommand = null;

	private static void guessJapsaPath(){
		if (japsaPath == null){
			boolean isWindows = System.getProperty("os.name").toLowerCase().indexOf("win") >= 0;
			japsaPath = isWindows?
					"c:\\Japsa"
					: System.getProperty("user.home") + "/.usr/local";
		}
	}
	/**
	 * Prepare the directory to copy libraries and scripts for instalation.
	 * This method also set up classpath, java command, library path, and
	 * (default) memory allocation in the process
	 *
	 * Set up classPath and javaCommand
	 * @throws IOException
	 */
	public static void setUpDirectory() throws IOException{

		boolean isWindows = System.getProperty("os.name").toLowerCase().indexOf("win") >= 0;
		boolean isMac = System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0;

		classPath = japsaJar;
		Scanner scanner = new Scanner(System.in);
		String line = null;

		System.out.println("Setting up Japsa Directory and copying libraries");
		////////////////////////////////////////////////////////////////////////////
		if (japsaPath == null){
			guessJapsaPath();
			while (true){
				System.out.print("Directory to install japsa: [" + japsaPath + "]");
				line = scanner.nextLine();
				line = line.trim();

				if (line.length() > 0){
					japsaPath = line;
				}

				classPath =  setupJapsaDir(japsaJar);
				if (classPath != null)
					break;//while
			}//while
		}else{
			classPath =  setupJapsaDir(japsaJar);
			if (classPath == null)
				System.exit(-1);
		}
		////////////////////////////////////////////////////////////////////////////
		//Get default memory if not previously set
		if (maxMem == null){
			maxMem = isWindows? "1000m":"7000m";
			System.out.print("Default memory allocated to jvm: [" + maxMem + "]");
			line = scanner.nextLine();
			line = line.trim();
			if (line.length() > 0){
				maxMem = line;
			}
		}

		javaCommand = isWindows?
				"java -Xmx%JSA_MEM% -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8"
				:"java -Xmx${JSA_MEM} -ea -Djava.awt.headless=true -Dfile.encoding=UTF-8";


		//Get server mode or client mode
		////////////////////////////////////////////////////////////////////////////
		if (server < 0){
			System.out.print("Enforce your jvm to run on server mode: [n]?(y/n)");
			line = scanner.nextLine();
			line = line.trim().toLowerCase();
			server = (line.equals("y") || line.equals("yes"))?1:0;
		}

		if (server > 0)
			javaCommand += " -server";

		scanner.close();
		System.out.println("Finished copying libraries\n");
	}
	/**
	 *
	 * @param masterScript
	 * @throws IOException
	 */
	public static void setUpScripts(ArrayList<Object> toolList, String masterScript)
			throws IOException{
		System.out.println("Set upting scripts in " + masterScript + ":");
		boolean isWindows = System.getProperty("os.name").toLowerCase().indexOf("win") >= 0;
		//Set up differences between windows and the rest
		String suffixStr = isWindows?".bat":"";
		String echoStr = isWindows?"echo ":"";

		File outJsa = new File(japsaPath + File.separator +  "bin" + File.separator + masterScript + suffixStr);
		PrintStream outJsaMain = new PrintStream(new FileOutputStream(outJsa));
		if (!isWindows){
			outJsaMain.println("#!/bin/sh\n\ncat << EOF");
		}else{
			outJsaMain.println("@echo off");
		}

		outJsaMain.print(echoStr + "Japsa: A Java Package for Statistical Sequence Analysis\n"
				+ echoStr + "Version " + VERSION + ", Built on " + (new Date()));

		if (compiler != null){
			outJsaMain.println(" with " + compiler);
		}
		else
			outJsaMain.println();

		outJsaMain.println(echoStr);
		outJsaMain.println(echoStr + " List of tools:");

		for (Object obj : toolList) {

			//A string separated
			if ((obj instanceof String)){
				outJsaMain.println(echoStr);
				outJsaMain.printf(echoStr + "%s\n",obj);
				continue;
			}
			if (obj instanceof CommandLine){
				Class<?> tool = obj.getClass();

				Deployable annotation = (Deployable) tool.getAnnotation(Deployable.class);
				File file = new File(japsaPath + File.separator +  "bin" + File.separator + annotation.scriptName() + suffixStr);

				if(isWindows){
					PrintStream out = new PrintStream(new FileOutputStream(file));
					out.println("@echo off");

					out.println("if \"%JSA_MEM%\"==\"\" (set JSA_MEM=" + maxMem+")");

					out.println("set JSA_CP=" + classPath);
					out.println();
					out.println(javaCommand +" -classpath %JSA_CP% " + tool.getCanonicalName() + " %*");
					out.close();
				}else{
					PrintStream out = new PrintStream(new FileOutputStream(file));
					out.println("#!/bin/sh\n");
					out.println("case $JSA_MEM in\n  '')JSA_MEM="+maxMem +";;\n  *);;\nesac\n\n");
					out.println("case $JSA_CP in\n  '')JSA_CP="
							+ classPath
							+ ";;\n  *)echo \"[INFO] Use ${JSA_CP} as path \" 1>&2;;\nesac\n\n");

					//out.println("JSA_CMD=\"`basename $0` $@\"\n");

					out.println(javaCommand + " -classpath ${JSA_CP} "
							+ tool.getCanonicalName() + " \"$@\"");
					out.close();

					Runtime.getRuntime().exec(
							"chmod a+x " + file.getCanonicalPath());
				}
				System.out.println(" " + file.getCanonicalPath() + " created");
				outJsaMain.printf(echoStr + "  %-23s  %s\n", annotation.scriptName(),	annotation.scriptDesc());

				//CommandLine cmdTool = (CommandLine) obj;
				//if (cmdTool.galaxy()){
				//	String wrapper = cmdTool.generateGalaxyWrapper();
				//	System.out.println(wrapper);
				//}

			}//if
			else{
				System.err.println("Cannot generate script for " + obj + "  " + obj.getClass());
			}

		}
		//Done
		if (!isWindows){
			outJsaMain.println("\nEOF");

		}
		outJsaMain.close();
		if (!isWindows){
			Runtime.getRuntime().exec(
					"chmod a+x " + outJsa.getCanonicalPath());
		}
		System.out.println("Done " + masterScript + "\n");
	}


	public static void setUpGalaxyScripts(ArrayList<Object> toolList)
			throws IOException{
		System.out.println("Set galaxy wrapper :");
		PrintStream masterFile = new PrintStream(new File("galaxy"  + File.separator + "japsa.xml"));
		masterFile.println("  <section id=\"japsa_tools\" name=\"Japsa Tools\">");

 		for (Object obj : toolList) {
			//A string separated
			if ((obj instanceof String)){
				continue;
			}
			if (obj instanceof CommandLine){

				CommandLine cmdTool = (CommandLine) obj;
				if (cmdTool.galaxy()){

					Class<?> tool = obj.getClass();
					Deployable annotation = (Deployable) tool.getAnnotation(Deployable.class);

					String fileName = "japsa" + File.separator + annotation.scriptName().replace(".", "_") + ".xml";
					File file = new File("galaxy"  + File.separator  + fileName);
					String wrapper = cmdTool.generateGalaxyWrapper();

					PrintStream ps = new PrintStream(file);
					ps.println(wrapper);
					ps.close();

					masterFile.println("    <tool file=\"" + fileName + "\" />");

					masterFile.println("  </section>");


					System.out.println(" " + file.getCanonicalPath() + " created");

				}
			}
		}
		masterFile.println("  </section>");
		masterFile.close();

		System.out.println("Done galaxy \n");
	}




	public static boolean uninstallLibraries() throws IOException{
		guessJapsaPath();

		if (japsaPath.startsWith("~/")) {
			japsaPath = System.getProperty("user.home") + japsaPath.substring(1);
		}

		File japsaLib = new File(japsaPath + File.separator + "lib" + File.separator + "japsa");
		File japsaMaster = new File(japsaPath + File.separator + "bin" + File.separator + "jsa");
		if (!japsaLib.exists() || !japsaMaster.exists() ){
			System.err.println("No instalation of japsa found at " + japsaPath);
			return false;
		}

		File to = new File (japsaLib.getCanonicalPath()  + File.separator + japsaJar);
		to.delete();

		StringSeparator ss = new StringSeparator(libs, ':');
		while (ss.hasNext()) {
			String l = ss.next();
			if (l.length() > 0){
				to = new File (japsaLib.getCanonicalPath()  + File.separator + l);
				to.delete();
			}
		}
		return true;

	}


	public static void uninstallScripts(ArrayList<Object> toolList, String masterScript) throws IOException{
		boolean isWindows = System.getProperty("os.name").toLowerCase().indexOf("win") >= 0;
		//Set up differences between windows and the rest
		String suffixStr = isWindows?".bat":"";

		// Delete all the scripts
		for (Object obj : toolList) {
			if (!(obj instanceof CommandLine)){
				continue;
			}

			Class<?> tool = obj.getClass();
			Deployable annotation = tool.getAnnotation(Deployable.class);

			File file = new File(japsaPath + File.separator +  "bin" + File.separator + annotation.scriptName());
			System.out.println("rm " + file.getCanonicalPath());
			file.delete();
		}
		File jsa = new File(japsaPath + File.separator +  "bin" + File.separator + masterScript + suffixStr);
		System.out.println("rm " + jsa.getCanonicalPath());
		jsa.delete();
	}


	public static void main(String[] args) throws NoSuchFieldException,
	SecurityException, IOException {
		CommandLine cmdLine = new CommandLine();
		cmdLine.addString("mode", "install", "install or uinstall");
		cmdLine.addString("libs", "", "list of extenal libraries");
		cmdLine.addString("installDir", null, "the directory to install");
		// cmdLine.addString("jlp", null, "Directories to libhdf5");
		cmdLine.addString("xmx", null, "Set default maximum memory");
		cmdLine.addString("compiler", null, "Compiler version");
		cmdLine.addBoolean("version", false, "Get version and exit");
		cmdLine.addString("server", "na", "Run on server mode: yes/true for yes; no/false for no");

		cmdLine.addStdHelp();
		args = cmdLine.stdParseLine(args);
		/**********************************************************************/
		if (cmdLine.getBooleanVal("version")){
			System.out.println(VERSION);
			System.exit(0);
		}

		///Get command lines option
		String mode = cmdLine.getStringVal("mode");

		japsaPath = cmdLine.getStringVal("installDir");
		compiler  = cmdLine.getStringVal("compiler");
		// jlp       = cmdLine.getStringVal("jlp");
		libs      = cmdLine.getStringVal("libs");
		maxMem    = cmdLine.getStringVal("xmx");
		String serverOpt = cmdLine.getStringVal("server").toLowerCase();
		if(serverOpt.equals("yes") || serverOpt.equals("true"))
			server = 1;
		if(serverOpt.equals("no") || serverOpt.equals("false"))
			server = 0;


		japsaJar = "japsa.jar";
		/**********************************************************************/


		if ("install".equals(mode)) {
			setUpDirectory();
			setUpScripts(tools, "jsa");
			//System.out.println("Japsa installtion complete, please set your path to " + japsaPath + File.separator+"bin");
			System.out.println("Japsa installtion complete\nFor your convenience, please add the following directory your PATH: " + Deploy.japsaPath + File.separator+"bin\n");
		} else if ("uninstall".equals(mode)) {
			//japsaPath must have been set
			if (uninstallLibraries())
				uninstallScripts(tools, "jsa");
		}else if ("galaxy".equals(mode)) {
			setUpGalaxyScripts(tools);
		}
		else {
			System.err.println("Mode " + mode + " not recognised");
			System.err.println(cmdLine.errors() + "\n" + "Usage: " + cmdLine.usage() + "\nOptions:\n" + cmdLine.options());
			System.exit(-1);
		}
	}
}
