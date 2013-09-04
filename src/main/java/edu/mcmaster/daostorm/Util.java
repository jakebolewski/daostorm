package edu.mcmaster.daostorm;

import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.io.Opener;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.Collections;

public class Util {

    public static void printStats(ArrayList<Peak> peaks) {
        int[] stats = fitStats(peaks);
        double total = (double) stats[3];
        System.out.format("%d good:%t%d bad:%t%d unconverged:%t%d total:%n",
                          stats[0] / total, stats[1] / total, stats[2] / total, (int) total);
    }

    public static double[] calcSxSy(double[] wxParams,
                                    double[] wyParams,
                                    double z)
    {
        assert(wxParams != null && wyParams != null);
        assert((wxParams.length == 5) &&
               (wxParams.length == wyParams.length));
        double zx = (z - wxParams[1]) / wxParams[2];
        double sx = 0.5 * wxParams[0] * Math.sqrt(1.0 + zx*zx + wxParams[3]*zx*zx*zx + wxParams[4]*zx*zx*zx*zx);
        double zy = (z - wyParams[1]) / wyParams[2];
        double sy = 0.5 * wyParams[0] * Math.sqrt(1.0 + zy*zy + wyParams[3]*zy*zy*zy + wyParams[4]*zy*zy*zy*zy);
        return new double[]{sx, sy};
    }

    public static ArrayList<Peak> copyPeakList(ArrayList<Peak> peaks) {
        ArrayList<Peak> copy = new ArrayList<Peak>(peaks.size());
        for (Peak p : peaks) {
            copy.add(p.copyPeak());
        }
        return copy;
    }

    public static ArrayList<Peak> getConvergedPeaks(ArrayList<Peak> peaks,
                                                    double minHeight,
                                                    double minWidth) {
        assert(peaks != null);
        assert(minHeight >= 0.0 && minWidth >= 0.0);

        if (peaks.size() == 0)
            return peaks;

        ArrayList<Peak> convergedPeaks = new ArrayList<Peak>();
        minWidth *= 0.5;
        for (Peak p : peaks) {
            if (p.hasStatus(PeakStatus.CONVERGED) &&
                p.getHeight() > minHeight &&
                p.getXWidth() > minWidth &&
                p.getYWidth() > minWidth)
                convergedPeaks.add(p);
        }
        return convergedPeaks;
    }

    public static ArrayList<Peak> getGoodPeaks(ArrayList<Peak> peaks,
                                               double minHeight,
                                               double minWidth) {
        assert(peaks != null);
        assert(minHeight >= 0.0 && minWidth >= 0.0);

        if (peaks.size() == 0)
            return peaks;
        ArrayList<Peak> goodPeaks = new ArrayList<Peak>();
        minWidth *= 0.5;
        for (Peak p : peaks) {
            if (!p.hasStatus(PeakStatus.ERROR) &&
                 p.getHeight() > minHeight &&
                 p.getXWidth() > minWidth &&
                 p.getYWidth() > minWidth)
                goodPeaks.add(p);
        }
        return goodPeaks;
    }

    public static int[] fitStats(ArrayList<Peak> peaks) {
        assert(peaks != null);
        int total = peaks.size();
        int numBad = 0;
        int numConverged = 0;
        int numUnconverged = 0;
        for (Peak p : peaks) {
            PeakStatus status = p.getStatus();
            if (status == PeakStatus.BADPEAK)
                numBad++;
            else if (status == PeakStatus.CONVERGED)
                numConverged++;
            else if (status == PeakStatus.RUNNING)
                numUnconverged++;
        }
        return new int[]{numConverged, numBad, numUnconverged, total};
    }

    public static FloatProcessor subtractImageProcessors(FloatProcessor ip1,
                                                         FloatProcessor ip2) {
        assert(ip1 != null && ip2 != null);
        assert(ip1.getWidth() == ip2.getWidth() &&
               ip1.getHeight() == ip2.getHeight());
        FloatProcessor out = (FloatProcessor) ip1.duplicate();
        for (int y=0; y < ip1.getHeight(); y++) {
            for (int x=0; x < ip2.getWidth(); x++) {
                out.setf(x, y, ip1.getf(x, y) - ip2.getf(x, y));
            }
        }
        return out;
    }

    public static ImageProcessor subtractImageProcessors(ImageProcessor ip1,
                                                         ImageProcessor ip2) {
        assert(ip1 != null && ip2 != null);
        assert(ip1.getWidth() == ip2.getWidth() &&
               ip1.getHeight() == ip2.getHeight());
        ImageProcessor out = ip1.duplicate();
        for (int y=0; y < ip1.getHeight(); y++) {
            for (int x=0; x < ip2.getWidth(); x++) {
                out.set(x, y, ip1.get(x, y) - ip2.get(x, y));
            }
        }
        return out;
    }

    public static FloatProcessor addImageProcessors(FloatProcessor ip1,
                                                    FloatProcessor ip2) {
        assert(ip1 != null && ip2 != null);
        assert(ip1.getWidth() == ip2.getWidth() &&
               ip1.getHeight() == ip2.getHeight());
        FloatProcessor out = (FloatProcessor) ip1.duplicate();
        for (int y=0; y < ip1.getHeight(); y++) {
            for (int x=0; x < ip2.getWidth(); x++) {
                out.setf(x, y, ip1.getf(x, y) + ip2.getf(x, y));
            }
        }
        return out;
    }

    public static ImageProcessor addImageProcessors(ImageProcessor ip1,
                                                    ImageProcessor ip2) {
        assert(ip1 != null && ip2 != null);
        assert(ip1.getWidth() == ip2.getWidth() &&
               ip1.getHeight() == ip2.getHeight());
        ImageProcessor out = ip1.duplicate();
        for (int y=0; y < ip1.getHeight(); y++) {
            for (int x=0; x < ip2.getWidth(); x++) {
                out.set(x, y, ip1.get(x, y) + ip2.get(x, y));
            }
        }
        return out;
    }

    public static void addPeakROIS(ImagePlus imp, ArrayList<Peak> peaks) {
        assert(imp != null && peaks != null);
        Overlay overlay = new Overlay();
        for (Peak peak : peaks) {
            PointRoi pt = new PointRoi(peak.getXCenter() + 0.5,
                                       peak.getYCenter() + 0.5);
            overlay.add(pt);
        }
        imp.setOverlay(overlay);
    }

    public static void main(String... args) {
        ImageJ ij = new ImageJ();
        ImagePlus imp = new Opener().openImage("/home/jake/STORM_DATA", "comp.tif");
        ImageProcessor improc = imp.getProcessor();
        float[][] imgArrayFloat = improc.getFloatArray();

        int margin = 10;
        int imageWidth = improc.getWidth();
        int imageHeight = improc.getHeight();
        double[] imgArrayDouble = new double[imageWidth * imageHeight];

        // ImageJ defines 0,0 as the bottom LEFT hand corner...
        for (int i=0; i < improc.getHeight(); i++) {
            for (int j=0; j < improc.getWidth(); j++) {
                int idx = i * imageWidth + j;
                // need to flip axes to align with imagej's layout
                imgArrayDouble[idx] = (double) imgArrayFloat[j][i];
            }
        }

        /* Fast copy (not transposed)
        imgWidth = improc.getWidth();
        for (int i=0; i < improc.getHeight(); i++) {
            System.arrayCopy(imgArrayFloat[i], 0, imgArrayDouble, i * imgWidth, imgWidth);
        }
        */

        int[] taken = new int[imgArrayDouble.length];
        double threshold = 300.0;
        double radius = 3.0;
        double background = 100.0;
        double sigma = 1.2;

        ArrayList<Peak> peaks = findLocalMaxima(imgArrayDouble,
                                                taken,
                                                threshold,
                                                radius,
                                                background,
                                                sigma,
                                                imageWidth,
                                                imageHeight,
                                                margin);
        System.out.format("Found %d peaks....%n", peaks.size());
        addPeakROIS(imp, peaks);
        imp.show();
    }

    /*
     * findLocalMaxima(image, x, y, h, threshold, background, sigma, image_size, margin, peak_size)
     *
     * Finds the locations of all the local maxima in an image with
     * intensity greater than threshold. Adds them to the list if
     * that location has not already been used.
     *
     * image - the image to analyze, assumed square.
     * taken - spots in the image where peaks have already been
     * peak - pre-allocated storage for peak data.
     * threshold - minumum peak intensity.
     * radius - circle in which the peak is maximal.
     * background - initial value for peak background.
     * sigma - initial value for peak width.
     * image_size_x - size of the image in x (fast axis).
     * image_size_y - size of the image in y (slow axis).
     * margin - number of pixels around the edge to ignore.
     * peak_size - size of the peaks array.
     *
     * Returns the number of peaks found.
     */
    public static ArrayList<Peak> findLocalMaxima(double[] image,
                                                  int[] taken,
                                                  double threshold,
                                                  double radius,
                                                  double background,
                                                  double sigma,
                                                  int imageSizeX,
                                                  int imageSizeY,
                                                  int margin)
    {
        ArrayList<Peak> newPeaks = new ArrayList<Peak>();

        int r = (int) Math.ceil(radius);
        double radiusSq = radius * radius;
        for (int i=margin; i < (imageSizeY - margin); i++) {
            for (int j=margin; j < (imageSizeX - margin); j++) {
                int idx = i * imageSizeX + j;
                double pixel = image[idx];
                boolean isMax = false;
                if (pixel > threshold) {
                    isMax = true;
                    int k = -r;
                    while ((k <= r) && isMax) {
                        int idx2 = idx + k * imageSizeX;
                        for (int l=-r; l <= r; l++) {
                            if ((k*k + l*l) < radiusSq) {
                                if ((k <= 0) && (l <= 0)) {
                                    if (pixel < image[idx2 + l]) {
                                        isMax = false;
                                    }
                                } else {
                                    if (pixel <= image[idx2 + l]) {
                                        isMax = false;
                                    }
                                }
                            }
                        }
                        k++;
                    }
                    if (isMax) {
                        if (taken[idx] < 2) {
                            Peak newPeak = new Peak(
                                                pixel - background, // height
                                                background,         // background
                                                j,                  // x center
                                                i,                  // y center
                                                sigma,              // width x
                                                sigma);             // width y
                            newPeaks.add(newPeak);
                            taken[idx] += 1;
                        }
                    }
                }
            }
        }
        return newPeaks;
    }

    public static double[] peakToPeakDist(ArrayList<Peak> peaks1,
                                          ArrayList<Peak> peaks2) {
        double[] bestDists = new double[peaks1.size()];
        int idx = 0;
        for (Peak p1 : peaks1) {
            double bestDist = Double.MAX_VALUE;
            for (Peak p2 : peaks2) {
                double dist = p1.sqDist2D(p2);
                if (dist < bestDist) {
                    bestDist = dist;
                }
            }
            bestDists[idx] = Math.sqrt(bestDist);
            idx++;
        }
        return bestDists;
    }

    public static int[] peakToPeakIndex(ArrayList<Peak> peaks1,
                                        ArrayList<Peak> peaks2) {
        int[] index = new int[peaks1.size()];
        for (Peak p1 : peaks1) {
            int idx = 0;
            int bestIdx = 0;
            double bestDist = p1.sqDist2D(peaks2.get(0));
            for (Peak p2 : peaks2) {
                double dist = p1.sqDist2D(p2);
                if (dist < bestDist) {
                    bestDist = dist;
                    bestIdx = idx;
                }
                idx++;
            }
            index[idx] = bestIdx;
        }
        return index;
    }

    public static ArrayList<Peak> mergeNewPeaks(ArrayList<Peak> inPeaks,
                                                ArrayList<Peak> newPeaks,
                                                double radius,
                                                double neighborhood)
    {
        int numInPeaks = inPeaks.size();
        ArrayList<Peak> outPeaks = new ArrayList<Peak>(numInPeaks);
        outPeaks.addAll(inPeaks);

        // check new peaks and add if they are ok
        double radiusSq = radius * radius;
        double neighborhoodSq = neighborhood * neighborhood;
        for (int i=0; i < numInPeaks; i++) {
            Peak newPeak = newPeaks.get(i);
            boolean bad = false;
            int j = 0;
            while ((j < numInPeaks) && (!bad)) {
                Peak inPeak = inPeaks.get(j);
                double distSq = newPeak.sqDist2D(inPeak);
                if (distSq < radiusSq) {
                    bad = true;
                } else if (distSq < neighborhoodSq) {
                    // Fixme: this could mark as running peaks that are close to a bad peak which do not get added
                    Peak outPeak = outPeaks.get(j);
                    outPeak.setZCenter(0.0);
                    outPeak.setStatus(PeakStatus.RUNNING);
                }
                j++;
            }
            if (!bad) outPeaks.add(newPeak);
        }
        return outPeaks;
    }


    // TODO: Spatial data structures could really speed this up
    public static ArrayList<Peak> removeClosePeaks(ArrayList<Peak> inPeaks,
                                                   double radius,
                                                   double neighborhood)
    {
        ArrayList<Peak> outPeaks = new ArrayList<Peak>();
        int numInPeaks = inPeaks.size();
        double radiusSq = radius * radius;
        double neighborhoodSq = neighborhood * neighborhood;

        // 1. flag peaks to be removed
        for (int i=0; i < numInPeaks; i++) {
            Peak p1 = inPeaks.get(i);
            boolean bad = false;
            int j = 0;
            while ((j < numInPeaks) && (!bad)) {
                if (j != i) {
                    Peak p2 = inPeaks.get(j);
                    if ((p1.sqDist2D(p2) < radiusSq) &&
                        (p2.getHeight() > p1.getHeight()))
                        bad = true;
                }
            }
            if (bad) p1.setStatus(PeakStatus.BADPEAK);
        }

        // 2. flag non-bad neighbors of bad peaks as running
        for (Peak p1 : inPeaks) {
            if (p1.getStatus() == PeakStatus.BADPEAK) {
                for (Peak p2 : inPeaks) {
                    if (!p1.equals(p2)) {
                        if (p1.dist2D(p2) < neighborhoodSq &&
                            p2.getStatus() != PeakStatus.BADPEAK) {
                            p2.setStatus(PeakStatus.RUNNING);
                        }
                    }
                }
            }
        }

        // 3. create a new list with the bad peaks removed
        for (Peak p : inPeaks)
            if (p.getStatus() != PeakStatus.BADPEAK)
                outPeaks.add(p);

        return outPeaks;
    }


    public static ArrayList<Peak> removeNeighbors(ArrayList<Peak> inPeaks,
                                                  double radius)
    {
        ArrayList<Peak> outPeaks = new ArrayList<Peak>();
        double radius2 = radius * radius;
        int numPeaks = inPeaks.size();
        for (int i=0; i < numPeaks; i++) {
            Peak p1 = inPeaks.get(i);
            boolean bad = false;
            int j = 0;
            while ((j < numPeaks) && (!bad)) {
                 if (j != i) {
                     Peak p2 = inPeaks.get(j);
                     if (p1.sqDist2D(p2) < radius2)
                         bad = true;
                 }
                 j++;
            }
            if (!bad) outPeaks.add(p1);
        }
        return outPeaks;
    }
}
