package edu.mcmaster.daostorm;

import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.io.Opener;

import java.io.File;
import java.util.Random;

public class LoadAndDisplayDatasetIJ1 {

    public static void main(String...  args) {
        ImageJ ij = new ImageJ();
        File file = new  File("/home/jake/STORM_DATA/comp.tif");
        ImagePlus imp = new Opener().openImage(
                            file.getAbsolutePath());

        Random rand = new Random(123);
        Overlay overlay = new Overlay();
        int width = imp.getWidth();
        int height = imp.getHeight();
        for (int i=0; i < 1000; i++) {
            double x = ((double) width) * rand.nextDouble();
            double y = ((double) height) * rand.nextDouble();
            PointRoi pt = new PointRoi(x, y);
            overlay.add(pt);
        }
        imp.setOverlay(overlay);
        imp.show();
        /*
        imp.show();
        Img imp2 = ImagePlusAdapter.wrap(imp);
        ImageJFunctions.show(imp2);
        */
    }
}
