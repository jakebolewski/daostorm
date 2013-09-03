package edu.mcmaster.daostorm;

public class FitPeak {

    public static final int NPEAKPAR = 7;

    public byte[] sign;
    public double[] clamp;

    public int offset;
    public int wx;
    public int wy;
    public int xc;
    public int yc;
    public double error;
    public double errorOld;

    // Z Fitting parameters
    public double wxTerm;
    public double wyTerm;
    public double[] xt;
    public double[] ext;
    public double[] yt;
    public double[] eyt;

    public Peak peak;

    public FitPeak(Peak p, int margin) {
        this.peak = p;

        this.sign = new byte[Peak.NPARAMS];
        this.clamp = new double[Peak.NPARAMS];

        this.xt = new double[2 * margin + 1];
        this.yt = new double[2 * margin + 1];
        this.ext = new double[2 * margin + 1];
        this.eyt = new double[2 * margin + 1];
    }

    public FitPeak(int margin) {
        assert(margin >= 5);
        this.peak = new Peak();
        this.sign = new byte[Peak.NPARAMS];
        this.clamp = new double[Peak.NPARAMS];

        this.xt = new double[2 * margin + 1];
        this.yt = new double[2 * margin + 1];
        this.ext = new double[2 * margin + 1];
        this.eyt = new double[2 * margin + 1];
    }

    public void setClampXCenter(double c) {
        this.clamp[Peak.XCENTER] = c;
    }

    public double getClampXCenter() {
        return this.clamp[Peak.XCENTER];
    }

    public void setClampYCenter(double c) {
        this.clamp[Peak.YCENTER] = c;
    }

    public double getClampYCenter() {
        return this.clamp[Peak.YCENTER];
    }

    public void setClampZCenter(double c) {
        this.clamp[Peak.ZCENTER] = c;
    }

    public double getClampZCenter() {
        return this.clamp[Peak.ZCENTER];
    }

    public void setClampHeight(double c) {
        this.clamp[Peak.HEIGHT] = c;
    }

    public double getClampHeight() {
        return this.clamp[Peak.HEIGHT];
    }

    public void setClampBackground(double c) {
        this.clamp[Peak.BACKGROUND] = c;
    }

    public double getClampBackground() {
        return this.clamp[Peak.BACKGROUND];
    }

    public void setClampXWidth(double c) {
        this.clamp[Peak.XWIDTH] = c;
    }

    public double getClampXWidth() {
        return this.clamp[Peak.XWIDTH];
    }

    public void setClampYWidth(double c) {
        this.clamp[Peak.YWIDTH] = c;
    }

    public double getClampYWidth() {
        return this.clamp[Peak.YWIDTH];
    }

}
