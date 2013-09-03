package edu.mcmaster.daostorm;

import imagej.ImageJ;
import imagej.data.Dataset;
import imagej.ui.UIService;

import javax.swing.*;
import java.io.File;


public class LoadAndDisplayDatasetJ2 {

    public static void main(final String... args) throws Exception {
        // create the ImageJ application context with
        // all available services
        final ImageJ ij = new ImageJ();

        // ask the user for a file to open
        final JFileChooser chooser = new JFileChooser();
        final int ret = chooser.showOpenDialog(null);
        if (ret != JFileChooser.APPROVE_OPTION)
            return;
        final File file = chooser.getSelectedFile();

        // load the dataset
        final Dataset dataset = ij.dataset().open(file.getAbsolutePath());

        UIService ui = ij.ui();

        // display the dataset
        ui.show(dataset);
    }

}
