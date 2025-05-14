package utilities;

import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

public class HBILogger {
    private static Logger logger;

    static {
        try {
            // Create or get the logger
            logger = Logger.getLogger(HBILogger.class.getName());
            logger.setUseParentHandlers(false); // Disable default console handler

            // Console handler
            ConsoleHandler consoleHandler = new ConsoleHandler();
            consoleHandler.setLevel(Level.INFO);
            logger.addHandler(consoleHandler);

            FileHandler fileHandler = new FileHandler("HBILogger.log", true); // true = append mode
            fileHandler.setLevel(Level.ALL);
            fileHandler.setFormatter(new SimpleFormatter());
            logger.addHandler(fileHandler);

            // Set global logging level
            logger.setLevel(Level.ALL);

        } catch (Exception e) {
            System.err.println("Failed to initialize logger: " + e.getMessage());
        }
    }

    public static void info(String msg) {
        logger.info(msg);
    }

    public static void warning(String msg) {
        logger.warning(msg);
    }

    public static void error(String msg) {
        logger.severe(msg);
    }

    public static void debug(String msg) {
        logger.fine(msg);
    }

    public static void trace(String msg) {
        logger.finest(msg);
    }

}

