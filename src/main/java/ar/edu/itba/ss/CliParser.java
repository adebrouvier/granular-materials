package ar.edu.itba.ss;

import org.apache.commons.cli.*;

import static java.lang.System.exit;

public class CliParser {

    static double width = 1;
    static double height = 2;
    static double opening = 0.15;
    static double time = 10;
    static double dt2;

    private static Options createOptions(){
        Options options = new Options();
        options.addOption("h", "help", false, "Shows this screen.");
        options.addOption("w", "width", true, "Width of the silo.");
        options.addOption("H", "height", true, "Height of the silo.");
        options.addOption("d", "opening", true, "Size of opening.");
        options.addOption("t", "time", true, "Total time of the simulation.");
        options.addOption("dt2", "deltatime2", true, "Time step for the animation.");
        return options;
    }

    public static void parseOptions(String[] args){
        Options options = createOptions();
        CommandLineParser parser = new BasicParser();

        try{
            CommandLine cmd = parser.parse(options, args);
            if(cmd.hasOption("h")){
                help(options);
            }

            if (cmd.hasOption("H")) {
                height = Double.parseDouble(cmd.getOptionValue("H"));
            }

            if (cmd.hasOption("w")) {
                width = Double.parseDouble(cmd.getOptionValue("w"));
            }

            if (cmd.hasOption("d")) {
                opening = Double.parseDouble(cmd.getOptionValue("d"));
            }

            if (cmd.hasOption("t")) {
                time = Double.parseDouble(cmd.getOptionValue("t"));
            }

            if (cmd.hasOption("dt2")) {
                dt2 = Double.parseDouble(cmd.getOptionValue("dt2"));
            }

            if (!(height > width) && !(width > opening)){
                System.err.println("d < W < h");
            }
        }catch (Exception e){
            System.out.println("Command not recognized.");
            help(options);
        }
    }

    private static void help(Options options){
        HelpFormatter helpFormatter = new HelpFormatter();
        helpFormatter.printHelp("granular-materials", options);
        exit(0);
    }
}
