package edu.mcmaster.daostorm;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;
import java.util.Arrays;

public class MultiFit {

    private final double HYSTERESIS = 0.6;
    private final int MARGIN = 10;

    // fits tolerance
    private double tolerance;

    // imagedata
    private double[] imageData;
    private int imageSizeX;
    private int imageSizeY;

    // background data
    private double[] foregroundData;

    // fits (foreground) data
    private double[] backgroundData;

    // number of peaks covering a particular pixel
    private int[] backgroundCounts;

    // Z fitting
    private double[] wxZParams;
    private double[] wyZParams;
    private double minZ = -0.5; //um
    private double maxZ = +0.5; //um

    private ArrayList<FitPeak> fits;

    //-------------------------------------------
    // Constructor
    //-------------------------------------------
    public MultiFit(double[] image, ArrayList<Peak> params, double tolerance,
                    int imageSizeX, int imageSizeY, boolean zFitting) {
        assert(tolerance < 1.0e-1);
        assert(image != null && image.length == (imageSizeX * imageSizeY));

        this.tolerance = tolerance;

        this.imageData = image;
        this.imageSizeX = imageSizeX;
        this.imageSizeY = imageSizeY;

        this.foregroundData = new double[image.length];
        this.backgroundData = new double[image.length];
        this.backgroundCounts = new int[image.length];

        this.fits = new ArrayList<FitPeak>(params.size());
        for (Peak p : params) {
            FitPeak newFit = new FitPeak(p, 10);
            if (newFit.peak.getStatus() == PeakStatus.RUNNING) {
                newFit.error = 0.0;
                newFit.errorOld = 0.0;
            } else {
                newFit.error = newFit.peak.getIError();
                newFit.errorOld = newFit.error;
            }
            if (zFitting) {
                calcWidthsFromZ(newFit);
            } else {
                newFit.peak.setXWidth(
                    1.0 / (2.0 * (Math.pow(newFit.peak.getXWidth(), 2))));
                newFit.peak.setYWidth(
                    1.0 / (2.0 * (Math.pow(newFit.peak.getYWidth(), 2))));
            }

            newFit.xc = (int) newFit.peak.getXCenter();
            newFit.yc = (int) newFit.peak.getYCenter();
            // possible bug casting fromm double to int in original version
            newFit.wx = calcWidth(newFit.peak.getXWidth(), -10);
            newFit.wy = calcWidth(newFit.peak.getYWidth(), -10);

            // todo these are annoying constants
            newFit.setClampHeight(1000.0);
            newFit.setClampBackground(100.0);
            newFit.setClampXCenter(1.0);
            newFit.setClampYCenter(1.0);
            newFit.setClampXWidth(0.3);
            newFit.setClampYWidth(0.3);
            newFit.setClampZCenter(0.1);

            this.fits.add(newFit);
        }
        calcFit();
        calcError();
    }

    // -------------------------------------------------
    // Public Interface
    //--------------------------------------------------

    // initialize Z parameters
    public void initZParameters(double[] wxVsZ,
                                double[] wyVsZ,
                                double zMin,
                                double zMax) {
        assert(wxVsZ.length == 5 && wyVsZ.length == 5);
        assert(zMin < zMax);

        this.wxZParams = Arrays.copyOf(wxVsZ, 5);
        this.wyZParams = Arrays.copyOf(wyVsZ, 5);
        this.wxZParams[0] = Math.pow(this.wxZParams[0], 2);
        this.wyZParams[0] = Math.pow(this.wyZParams[0], 2);
        this.minZ = zMin;
        this.maxZ = zMax;
    }

    public void iter2DFixed() {
        update2DFixed();
        calcError();
    }

    // getError
    public double getTotalError() {
        double error = 0.0;
        for (FitPeak fit : this.fits) {
            error += fit.error;
        }
        return error;
    }

    public int getNumUnconverged() {
        int count = 0;
        for (FitPeak cur : this.fits) {
            if (cur.peak.getStatus() == PeakStatus.RUNNING) {
                count++;
            }
        }
        return count;
    }

    public void getResidual(double[] residual) {
        assert(residual.length == this.imageData.length);
        calcFit();
        for (int i=0; i < residual.length; i++) {
            residual[i] = this.imageData[i] - this.foregroundData[i];
        }
    }

    public ArrayList<Peak> getResults() {
        ArrayList<Peak> results = new ArrayList<Peak>(this.fits.size());
        for (FitPeak fit : this.fits) {
            if (fit.peak.getStatus() != PeakStatus.ERROR) {
                fit.peak.setXWidth(Math.sqrt(1.0 / (2.0 * fit.peak.getXWidth())));
                fit.peak.setYWidth(Math.sqrt(1.0 / (2.0 * fit.peak.getYWidth())));
            } else  {
                fit.peak.setXWidth(1.0);
                fit.peak.setYWidth(1.0);
            }
            fit.peak.setIError(fit.error);
            results.add(fit.peak);
        }
        return results;
    }

    // -----------------------------------------------
    // Private Interface
    // -----------------------------------------------
    private int calcWidth(double peakWidth, int oldWidth) {
        if (peakWidth < 0.0) {
            return 1;
        } else {
            int newWidth = oldWidth;
            double tmp = 4.0 * Math.sqrt(1.0 / (2.0 * peakWidth));
            if (Math.abs(tmp - ((double) oldWidth) - 0.0) > this.HYSTERESIS) {
                newWidth = (int) tmp;
            }
            newWidth = Math.min(newWidth, this.MARGIN);
            return newWidth;
        }
    }

    // TODO: make this a pure function
    private void calcWidthsFromZ(FitPeak fit) {

        double z0, z1, z2, z3, tmp;

        // wx
        z0 = (fit.peak.getZCenter() - this.wxZParams[1]) / this.wxZParams[2];
        z1 = z0 * z0;
        z2 = z1 * z0;
        z3 = z2 * z0;
        tmp = 1.0 + z1 + this.wxZParams[3]*z2 + this.wxZParams[4]*z3;
        fit.wxTerm = tmp * tmp;
        fit.peak.setXWidth(2.0 / (this.wxZParams[0] * tmp));

        // wy
        z0 = (fit.peak.getZCenter() - this.wyZParams[1]) / this.wyZParams[2];
        z1 = z0 * z0;
        z2 = z1 * z0;
        z3 = z2 * z0;
        tmp = 1.0 + z1 + this.wyZParams[3] * z2 + this.wyZParams[4] * z3;
        fit.wyTerm = tmp * tmp;
        fit.peak.setYWidth(2.0 / (this.wyZParams[0] * tmp));
    }


    private void fitDataUpdate(FitPeak fit, double[] deltas) {

        // update sign and clamp solution if appears to be oscillating
        for (int i=0; i < deltas.length; i++) {
            if (fit.sign[i] != 0) {
                if ((fit.sign[i] == 1) && (deltas[i] < 0.0)) {
                    fit.clamp[i] *= 0.5;
                } else if ((fit.sign[i] == -1) && (deltas[i] > 0.0)) {
                    fit.clamp[i] *= 0.5;
                }
            }
            if (deltas[i] > 0.0) {
                fit.sign[i] = 1;
            } else {
                fit.sign[i] = -1;
            }
            // update values based on delta and clamp
            if (deltas[i] != 0.0) {
                double update = deltas[i] /
                        (1.0 + Math.abs(deltas[i]) / fit.clamp[i]);
                fit.peak.addToParameter(i, -update);
            }
        }

        // update peak center with HYSTERESIS
        if (Math.abs(fit.peak.getXCenter() - ((double) fit.xc) - 0.5) > this.HYSTERESIS) {
            fit.xc = (int) fit.peak.getXCenter();
        }
        if (Math.abs(fit.peak.getYCenter() - ((double) fit.yc) - 0.5) > this.HYSTERESIS) {
            fit.yc = (int) fit.peak.getYCenter();
        }

        // check that the peak has not moved too close of the
        // edge of the image, flag peak as bad if it has
        int xc = fit.xc;
        int yc = fit.yc;

        if ((xc < this.MARGIN) ||
            (xc >= (imageSizeX - this.MARGIN)) ||
            (yc < this.MARGIN) ||
            (yc >= (imageSizeY - this.MARGIN))) {
            fit.peak.setStatus(PeakStatus.BADPEAK);
        }

        // check for negative background or height
        if ((fit.peak.getBackground() < 0.0) ||
            (fit.peak.getHeight() < 0.0)) {
            fit.peak.setStatus(PeakStatus.ERROR);
        }

        // check for negative widths
        if ((fit.peak.getXWidth() < 0.0) ||
            (fit.peak.getYWidth() < 0.0)) {
            fit.peak.setStatus(PeakStatus.ERROR);
        }
        // Opt 1 : peak errors if z is out of range
        /*
         if ((cur.peak.zCenter < this.minZ) ||
             (cur.peak.zCenter > this.maxZ)) {
             cur.status = PeakStatus.Error;
         }
        */
        // Opt 2: Clamp Z value Range
        if (fit.peak.getZCenter() < this.minZ)
            fit.peak.setZCenter(this.minZ);
        if (fit.peak.getZCenter() > this.maxZ)
            fit.peak.setZCenter(this.maxZ);
    }

    private void calcError() {
        for (FitPeak cur : this.fits) {
            if (cur.peak.getStatus() == PeakStatus.RUNNING) {
                int offset = cur.offset;
                int wx = cur.wx;
                int wy = cur.wy;
                double error = 0.0;

                for(int i=-wy; i <= wy; i++) {
                    for (int j=-wx; j <= wx; j++) {
                        int idx = (i * this.imageSizeX) + (j + offset);
                        double fi = this.foregroundData[idx] +
                                    (this.backgroundData[idx] / ((double) this.backgroundCounts[idx]));
                        double xi = this.imageData[idx];
                        error += 2 * (fi - xi) - 2 * xi * Math.log(fi / xi);
                    }
                }

                cur.errorOld = cur.error;
                cur.error = error;

                if (Math.abs(error - cur.errorOld) / error < this.tolerance) {
                    cur.peak.setStatus(PeakStatus.CONVERGED);
                }
            }
        }
    }

    private void calcFit() {
        Arrays.fill(this.foregroundData, 1.0);
        Arrays.fill(this.backgroundData, 0.0);
        Arrays.fill(this.backgroundCounts, 0);

        for (FitPeak fit : this.fits) {
            if (fit.peak.getStatus() != PeakStatus.ERROR) {
                addPeak(fit);
            }
        }
    }

    private void addPeak(FitPeak fit) {
        int xc = fit.xc;
        int yc = fit.yc;
        fit.offset = yc * this.imageSizeX + xc;

        int wx = fit.wx;
        int wy = fit.wy;

        double xCenter = fit.peak.getXCenter();
        double xWidth = fit.peak.getXWidth();
        for (int i=(xc - wx); i < (xc + wx); i++) {
            double xt = ((double) i) - xCenter;
            int n = i - xc + wx;
            fit.xt[n] = xt;
            fit.ext[n] = Math.exp(-xt * xt * xWidth);
        }

        double yCenter = fit.peak.getYCenter();
        double yWidth = fit.peak.getYWidth();
        for (int i=(yc - wy); i < (yc + wy); i++) {
            double yt = ((double) i) - yCenter;
            int n = i - yc + wy;
            fit.yt[n] = yt;
            fit.eyt[n] = Math.exp(-yt * yt * yWidth);
        }

        // gaussian function
        int offset = fit.offset;
        double background = fit.peak.getBackground();
        double height = fit.peak.getHeight();
        for (int i=-wy; i <= wy; i++) {
            double eyt = fit.eyt[i + wy];
            for (int j=-wx; j <= wx; j++) {
                double ext = fit.ext[j + wx];
                int idx = (i * this.imageSizeX) + (j + offset);
                this.foregroundData[idx] += height * eyt * ext;
                this.backgroundData[idx] += background;
                this.backgroundCounts[idx] += 1;
            }
        }
    }

    private void subtractPeak(FitPeak fit) {
        int wx = fit.wx;
        int wy = fit.wy;

        // gaussian function
        int offset = fit.offset;
        double background = fit.peak.getBackground();
        double height = fit.peak.getHeight();

        for (int i=-wy; i <= wy; i++) {
            double eyt = fit.eyt[i + wy];
            for (int j=-wx; j <= wx; j++) {
                double ext = fit.ext[j + wx];
                int idx = (i * this.imageSizeX) + (j + offset);
                this.foregroundData[idx] -= height * eyt * ext;
                this.backgroundData[idx] -= background;
                this.backgroundCounts[idx] -= 1;
            }
        }
    }

    public void updateZ() {

        double[] delta = new double[Peak.NFITPARAMS];
        double[] jt = new double[5];
        double[] jacobian = new double[5];
        double[][] hessian = new double[5][5];

        for (FitPeak fit : this.fits) {
            if (fit.peak.hasStatus(PeakStatus.RUNNING)) {

                for (int i=0; i < 5; i++) {
                    jacobian[i] = 0;
                    hessian[i][0] = 0.0;
                    hessian[i][1] = 0.0;
                    hessian[i][2] = 0.0;
                    hessian[i][3] = 0.0;
                    hessian[i][4] = 0.0;
                }

                int offset = fit.offset;
                int wx = fit.wx;
                int wy = fit.wy;

                double height = fit.peak.getHeight();
                double xwidth = fit.peak.getXWidth();
                double ywidth = fit.peak.getYWidth();

                // calculate dwx vs z
                double z0 = (fit.peak.getZCenter() - this.wxZParams[1]) /
                             this.wxZParams[2];
                double z1 = z0 * z0;
                double z2 = z1 * z0;
                double zt = (2.0 * z0) +
                            (3.0 * this.wxZParams[3] * z1) +
                            (4.0 * this.wxZParams[4] * z2);
                double gx = -2.0 * zt / (this.wxZParams[0] * fit.wxTerm);

                // calculate dwy vs z
                z0 = (fit.peak.getZCenter() - this.wxZParams[1]) /
                             this.wxZParams[2];
                z1 = z0 * z0;
                z2 = z1 * z0;
                zt = (2.0 * z0) +
                     (3.0 * this.wyZParams[3] * z1) +
                     (4.0 * this.wyZParams[4] * z2);
                double gy = -2.0 * zt / (this.wyZParams[0] * fit.wxTerm);

               for (int i=-wy; i < +wy; i++) {
                   double yt = fit.yt[i + wy];
                   double eyt = fit.eyt[i + wy];

                   for(int j=-wx; j < +wx; j++) {
                       double xt = fit.xt[j + wx];
                       double ext = fit.ext[j + wx];

                       int idx = (i * this.imageSizeX) + (j + offset);
                       double xi = this.imageData[idx];
                       double fi = this.foregroundData[idx] /
                               (this.backgroundData[idx] / ((double) this.backgroundCounts[idx]));

                       // first derivatives
                       double et = ext * eyt;
                       jt[0] = et;
                       jt[1] = 2.0 * height * xwidth * xt * et;
                       jt[2] = 2.0 * height * ywidth * yt * et;
                       jt[3] = (-height * xt * xt * gx * et) -
                               (height * yt * yt * gy * et);
                       jt[4] = 1.0;

                       // calculate jacobian
                       double t1 = 2.0 * (1.0 - xi/fi);
                       jacobian[0] += t1*jt[0];
                       jacobian[1] += t1*jt[1];
                       jacobian[2] += t1*jt[2];
                       jacobian[3] += t1*jt[3];
                       jacobian[4] += t1*jt[4];

                       // calculate hessian
                       double t2 = 2.0 * xi / (fi * fi);

                       // calculate hessian without second derivative terms.
                       hessian[0][0] += t2*jt[0]*jt[0];
                       hessian[0][1] += t2*jt[0]*jt[1];
                       hessian[0][2] += t2*jt[0]*jt[2];
                       hessian[0][3] += t2*jt[0]*jt[3];
                       hessian[0][4] += t2*jt[0]*jt[4];

                       hessian[1][1] += t2*jt[1]*jt[1];
                       hessian[1][2] += t2*jt[1]*jt[2];
                       hessian[1][3] += t2*jt[1]*jt[3];
                       hessian[1][4] += t2*jt[1]*jt[4];

                       hessian[2][2] += t2*jt[2]*jt[2];
                       hessian[2][3] += t2*jt[2]*jt[3];
                       hessian[2][4] += t2*jt[2]*jt[4];

                       hessian[3][3] += t2*jt[3]*jt[3];
                       hessian[3][4] += t2*jt[3]*jt[4];

                       hessian[4][4] += t2*jt[4]*jt[4];
                   }
                }

                // subtract the old peak out of the foreground and background arrays
                subtractPeak(fit);

                // use lapack to solve Ax=B to calculate update vector;
                boolean error = false;
                RealMatrix hessianMatrix = MatrixUtils.createRealMatrix(hessian);
                RealVector jacobianVector = MatrixUtils.createRealVector(jacobian);
                try {
                    MatrixUtils.solveUpperTriangularSystem(hessianMatrix, jacobianVector);
                } catch(Exception ex) {
                    fit.peak.setStatus(PeakStatus.ERROR);
                    error = true;
                }
                if (!error) {
                    // update parameters
                    delta[Peak.HEIGHT] = jacobianVector.getEntry(0);
                    delta[Peak.XCENTER] = jacobianVector.getEntry(1);
                    delta[Peak.YCENTER] = jacobianVector.getEntry(2);
                    delta[Peak.ZCENTER] = jacobianVector.getEntry(3);
                    delta[Peak.BACKGROUND] = jacobianVector.getEntry(4);

                    fitDataUpdate(fit, delta);

                    if (!fit.peak.hasStatus(PeakStatus.ERROR)) {
                        // calculate new x, y, width and update fit area
                        calcWidthsFromZ(fit);
                        fit.wx = calcWidth(fit.peak.getXWidth(), fit.wx);
                        fit.wy = calcWidth(fit.peak.getYWidth(), fit.wy);
                        addPeak(fit);
                    }
                }
            }
        }
    }


    public void update3D() {
        double[] delta = new double[Peak.NFITPARAMS];
        double[] jt = new double[6];
        double[] jacobian = new double[6];
        double[][] hessian = new double[6][6];

        for (FitPeak fit : this.fits) {
            if (fit.peak.hasStatus(PeakStatus.RUNNING)) {

                for (int i=0; i < 6; i++) {
                    jacobian[i] = 0.0;
                    hessian[i][0] = 0.0;
                    hessian[i][1] = 0.0;
                    hessian[i][2] = 0.0;
                    hessian[i][3] = 0.0;
                    hessian[i][4] = 0.0;
                    hessian[i][5] = 0.0;
                }

                int offset = fit.offset;
                int wx = fit.wx;
                int wy = fit.wy;
                double height = fit.peak.getHeight();
                double xwidth = fit.peak.getXWidth();
                double ywidth = fit.peak.getYWidth();

                for (int i=-wy; i < +wy; i++) {
                    double yt = fit.yt[i + wy];
                    double eyt = fit.yt[i + wy];

                    for (int j=-wx; j < + wx; j++) {
                        double xt = fit.xt[j + wx];
                        double ext = fit.ext[j + wx];

                        int idx = (i * this.imageSizeX) + (j + offset);

                        double xi = this.imageData[idx];
                        double fi = this.foregroundData[idx] +
                                (this.backgroundData[idx] / ((double) this.backgroundCounts[idx]));

                        double et = ext * eyt;
                        jt[0] = et;
                        jt[1] = 2.0 * height * xwidth * xt * et;
                        jt[2] = -height * xt * xt * et;
                        jt[3] = 2.0 * height * ywidth * yt * et;
                        jt[4] = -height * yt * yt * et;
                        jt[5] = 1.0;

                        // calculate jacobian
                        double t1 = 2.0*(1.0 - xi/fi);
                        jacobian[0] += t1*jt[0];
                        jacobian[1] += t1*jt[1];
                        jacobian[2] += t1*jt[2];
                        jacobian[3] += t1*jt[3];
                        jacobian[4] += t1*jt[4];
                        jacobian[5] += t1*jt[5];

                        // calculate hessian
                        double t2 = 2.0*xi/(fi*fi);
                                                                        // hessian without second derivative terms.
                        hessian[0][0] += t2*jt[0]*jt[0];
                        hessian[0][1] += t2*jt[0]*jt[1];
                        hessian[0][2] += t2*jt[0]*jt[2];
                        hessian[0][3] += t2*jt[0]*jt[3];
                        hessian[0][4] += t2*jt[0]*jt[4];
                        hessian[0][5] += t2*jt[0]*jt[5];

                        hessian[1][1] += t2*jt[1]*jt[1];
                        hessian[1][2] += t2*jt[1]*jt[2];
                        hessian[1][3] += t2*jt[1]*jt[3];
                        hessian[1][4] += t2*jt[1]*jt[4];
                        hessian[1][5] += t2*jt[1]*jt[5];

                        hessian[2][2] += t2*jt[2]*jt[2];
                        hessian[2][3] += t2*jt[2]*jt[3];
                        hessian[2][4] += t2*jt[2]*jt[4];
                        hessian[2][5] += t2*jt[2]*jt[5];

                        hessian[3][3] += t2*jt[3]*jt[3];
                        hessian[3][4] += t2*jt[3]*jt[4];
                        hessian[3][5] += t2*jt[3]*jt[5];

                        hessian[4][4] += t2*jt[4]*jt[4];
                        hessian[4][5] += t2*jt[4]*jt[5];

                        hessian[5][5] += t2*jt[5]*jt[5];

                    }
                }

                // subtract the old peak out of the foreground and background
                // arrays
                subtractPeak(fit);

                // use lapack to solve Ax=B to calculate update vector;
                boolean error = false;
                RealMatrix hessianMatrix = MatrixUtils.createRealMatrix(hessian);
                RealVector jacobianVector = MatrixUtils.createRealVector(jacobian);
                try {
                    MatrixUtils.solveUpperTriangularSystem(hessianMatrix, jacobianVector);
                } catch(Exception ex) {
                    fit.peak.setStatus(PeakStatus.ERROR);
                    error = true;
                }
                if (!error) {
                    delta[Peak.HEIGHT] = jacobianVector.getEntry(0);
                    delta[Peak.XCENTER] = jacobianVector.getEntry(1);
                    delta[Peak.XWIDTH] = jacobianVector.getEntry(2);
                    delta[Peak.YCENTER] = jacobianVector.getEntry(3);
                    delta[Peak.YWIDTH] = jacobianVector.getEntry(4);
                    delta[Peak.BACKGROUND] = jacobianVector.getEntry(5);

                    fitDataUpdate(fit, delta);

                    // add the new peak to foreground and background arrays
                    fit.wx = calcWidth(fit.peak.getXWidth(), fit.wx);
                    fit.wy = calcWidth(fit.peak.getYWidth(), fit.wy);
                    addPeak(fit);
                }
            }
        }
    }
    public void update2D() {
        double[] delta = new double[Peak.NFITPARAMS];
        double[] jt = new double[5];
        double[] jacobian = new double[5];
        double[][] hessian = new double[5][5];

        for (FitPeak fit : this.fits) {
            if (fit.peak.hasStatus(PeakStatus.RUNNING)) {

                for (int i=0; i < 5; i++) {
                    jacobian[i] = 0.0;
                    hessian[i][0] = 0.0;
                    hessian[i][1] = 0.0;
                    hessian[i][2] = 0.0;
                    hessian[i][3] = 0.0;
                    hessian[i][4] = 0.0;
                }

                int offset = fit.offset;
                int wx = fit.wx;
                int wy = fit.wy;
                double height = fit.peak.getHeight();
                double width = fit.peak.getXWidth();

                for (int i=-wy; i <= +wy; i++) {
                    double yt = fit.yt[i + wy];
                    double eyt = fit.eyt[i + wy];

                    for (int j=-wx; j <= +wx; j++) {
                        double xt = fit.xt[j + wx];
                        double ext = fit.ext[j + wx];

                        int idx = (i * this.imageSizeX) + (j + offset);
                        double xi = this.imageData[idx];
                        double fi = this.foregroundData[idx] +
                                (this.backgroundData[idx] / ((double) this.backgroundCounts[idx]));

                        double et = ext * eyt;
                        jt[0] = et;
                        jt[1] = 2.0 * height * width * xt * et;
                        jt[2] = 2.0 * height * width * yt * et;
                        jt[3] = (-height * xt * xt * et) - (height * yt * yt * et);
                        jt[4] = 1.0;

                        // calculate jacobian
                        double t1 = 2.0 * (1.0 - xi / fi);
                        jacobian[0] += t1*jt[0];
                        jacobian[1] += t1*jt[1];
                        jacobian[2] += t1*jt[2];
                        jacobian[3] += t1*jt[3];
                        jacobian[4] += t1*jt[4];

                        // calculate hessian
                        double t2 = 2.0 * xi / (fi*fi);

                        // calculate hessian without second derivative terms.
                        // (symmetric upper triangular)
                        hessian[0][0] += t2*jt[0]*jt[0];
                        hessian[0][1] += t2*jt[0]*jt[1];
                        hessian[0][2] += t2*jt[0]*jt[2];
                        hessian[0][3] += t2*jt[0]*jt[3];
                        hessian[0][4] += t2*jt[0]*jt[4];

                        hessian[1][1] += t2*jt[1]*jt[1];
                        hessian[1][2] += t2*jt[1]*jt[2];
                        hessian[1][3] += t2*jt[1]*jt[3];
                        hessian[1][4] += t2*jt[1]*jt[4];

                        hessian[2][2] += t2*jt[2]*jt[2];
                        hessian[2][3] += t2*jt[2]*jt[3];
                        hessian[2][4] += t2*jt[2]*jt[4];

                        hessian[3][3] += t2*jt[3]*jt[3];
                        hessian[3][4] += t2*jt[3]*jt[4];

                        hessian[4][4] += t2*jt[4]*jt[4];
                    }
                }

                // subtract the old peak out of the foreground and background arrays
                subtractPeak(fit);

                // use lapack to solve Ax=B to calculate update vector;
                boolean error = false;
                RealMatrix hessianMatrix = MatrixUtils.createRealMatrix(hessian);
                RealVector jacobianVector = MatrixUtils.createRealVector(jacobian);

                try {
                    MatrixUtils.solveUpperTriangularSystem(hessianMatrix, jacobianVector);
                } catch(Exception ex) {
                    fit.peak.setStatus(PeakStatus.ERROR);
                    error = true;
                }
                if (!error) {
                    // TODO: Rearranged the jacobian vector entries
                    delta[Peak.HEIGHT] = jacobianVector.getEntry(0);        // height
                    delta[Peak.XCENTER] = jacobianVector.getEntry(1);       // x center
                    delta[Peak.YCENTER] = jacobianVector.getEntry(2);       // y center
                    delta[Peak.XWIDTH] = jacobianVector.getEntry(3);        // width
                    delta[Peak.YWIDTH] = jacobianVector.getEntry(3);        // width
                    delta[Peak.BACKGROUND] = jacobianVector.getEntry(4);    // background

                    // update fits data
                    fitDataUpdate(fit, delta);

                    // add the peak to the foreground and background arrays
                    // recalculate peak fit area as the peak width may have changed
                    if (!fit.peak.hasStatus(PeakStatus.ERROR)) {
                        fit.wx = calcWidth(fit.peak.getXWidth(), fit.wx);
                        fit.wy = fit.wx;
                        addPeak(fit);
                    }
                }
            }
        }
    }


    public void update2DFixed() {
        double[] delta = new double[Peak.NFITPARAMS];
        double[] jt = new double[4];
        double[] jacobian = new double[4];
        double[][] hessian = new double[4][4];

        for (FitPeak fit : this.fits) {
            if (fit.peak.getStatus() == PeakStatus.RUNNING) {

                for (int i=0; i < 4; i++) {
                    jacobian[i] = 0.0;
                    hessian[i][0] = 0.0;
                    hessian[i][1] = 0.0;
                    hessian[i][2] = 0.0;
                    hessian[i][3] = 0.0;
                }

                int offset = fit.offset;
                int wx = fit.wx;
                int wy = fit.wy;
                double height = fit.peak.getHeight();
                double width = fit.peak.getXWidth();

                for (int i=-wy; i <= +wy; i++) {
                    double yt = fit.yt[i + wy];
                    double eyt = fit.eyt[i + wy];
                    for (int j=-wx; j <= +wx; j++) {
                        int idx = (i * this.imageSizeX) + (j + offset);

                        // est = foreground + average background
                        double fi = this.foregroundData[idx] +
                            (this.backgroundData[idx] / ((double) this.backgroundCounts[idx]));

                        double xi = this.imageData[idx];
                        double xt = fit.xt[j + wx];
                        double ext = fit.ext[j + wx];
                        double et = ext * eyt;

                        // TODO: Remember to change order to match peak array in peak class
                        jt[0] = et;
                        jt[1] = 2.0 * height * width * xt * et;
                        jt[2] = 2.0 * height * width * yt * et;
                        // TODO: this should be # 2
                        jt[3] = 1.0; // this is background

                        // calculate jacobian;
                        double t1 = 2.0 * (1.0 - xi / fi);
                        jacobian[0] += t1 * jt[0];
                        jacobian[1] += t1 * jt[1];
                        jacobian[2] += t1 * jt[2];
                        jacobian[3] += t1 * jt[3];

                        // calculate hessian (without second deriv terms)
                        double t2 = 2.0 * xi / (fi * fi);
                        hessian[0][0] += t2*jt[0]*jt[0];
                        hessian[0][1] += t2*jt[0]*jt[1];
                        hessian[0][2] += t2*jt[0]*jt[2];
                        hessian[0][3] += t2*jt[0]*jt[3];

                        // hessian[4]
                        hessian[1][1] += t2*jt[1]*jt[1];
                        hessian[1][2] += t2*jt[1]*jt[2];
                        hessian[1][3] += t2*jt[1]*jt[3];

                        // hessian[8]
                        // hessian[9]
                        hessian[2][2] += t2*jt[2]*jt[2];
                        hessian[2][3] += t2*jt[2]*jt[3];

                        // hessian[12]
                        // hessian[13]
                        // hessian[14]
                        hessian[3][3] += t2*jt[3]*jt[3];
                    }
                }

                // subtract old peak out of foreground and background arrays
                subtractPeak(fit);

                // use lapack to solve Ax=B to calculate update vector;
                boolean error = false;
                RealMatrix hessianMatrix = MatrixUtils.createRealMatrix(hessian);
                RealVector jacobianVector = MatrixUtils.createRealVector(jacobian);

                try {
                    MatrixUtils.solveUpperTriangularSystem(hessianMatrix, jacobianVector);
                } catch(Exception ex) {
                    fit.peak.setStatus(PeakStatus.ERROR);
                    error = true;
                }
                if (!error) {
                    // update parameters
                    // height
                    // TODO: Rearranged the jacobian vector entries
                    delta[Peak.HEIGHT] = jacobianVector.getEntry(0); // height
                    delta[Peak.BACKGROUND] = jacobianVector.getEntry(3); // background
                    delta[Peak.XCENTER] = jacobianVector.getEntry(1); // x center
                    delta[Peak.YCENTER] = jacobianVector.getEntry(2); // y center

                    // update fits data
                    fitDataUpdate(fit, delta);

                    // add the new peak to the foreground and background arrays
                    if (fit.peak.getStatus() != PeakStatus.ERROR) {
                        addPeak(fit);
                    }
                }
            }
        }
    }
}