package edu.mcmaster.daostorm;


// TODO: This class should be immutable
public class Peak {
    public static final int NPARAMS= 8;
    public static final int NFITPARAMS = 7;

    public static final int HEIGHT = 0;
    public static final int BACKGROUND = 1;
    public static final int XCENTER = 2;
    public static final int YCENTER = 3;
    public static final int XWIDTH = 4;
    public static final int YWIDTH = 5;
    public static final int ZCENTER = 6;
    public static final int IERROR = 7;

    private PeakStatus status;
    private double paramArray[];


    /*public double height;
    public double background;
    public double xCenter;
    public double yCenter;
    public double xWidth;
    public double yWidth;
    public double zCenter;
    public double integratedError;
    */

    // Default Constructor
    Peak() {
        this.paramArray = new double[NPARAMS];
        this.status = PeakStatus.NEW;
    }

    // 2D Constructor
    Peak(double height,
         double background,
         double xCenter,
         double yCenter,
         double xWidth,
         double yWidth)
    {
        this.paramArray = new double[NPARAMS];
        this.status = PeakStatus.RUNNING;

        this.paramArray[HEIGHT] = height;
        this.paramArray[BACKGROUND] = background;
        this.paramArray[XCENTER] = xCenter;
        this.paramArray[YCENTER] = yCenter;
        this.paramArray[XWIDTH] = xWidth;
        this.paramArray[YWIDTH] = yWidth;
    }

    // 3D Constructor
    Peak(double height,
         double background,
         double xCenter,
         double yCenter,
         double xWidth,
         double yWidth,
         double zCenter)
    {
        this.paramArray = new double[NPARAMS];
        this.status = PeakStatus.RUNNING;

        this.paramArray[HEIGHT] = height;
        this.paramArray[BACKGROUND] = background;
        this.paramArray[XCENTER] = xCenter;
        this.paramArray[YCENTER] = yCenter;
        this.paramArray[XWIDTH] = xWidth;
        this.paramArray[YWIDTH] = yWidth;
        this.paramArray[ZCENTER] = zCenter;
    }

    /* ---------- Private Interface ----------------*/
    private void setParamArray(double[] params) {
        assert(params.length == 5);
        this.paramArray = params;
    }

    /*----------- Public Interface ------------------*/
    public Peak copyPeak() {
        Peak copy = new Peak();
        copy.setParamArray(this.paramArray.clone());
        copy.setStatus(this.status);
        return copy;
    }

    public double[] getParameterArray() {
        return this.paramArray;
    }

    public boolean hasStatus(PeakStatus s) {
        return this.status == s;
    }

    public void addToParameter(int idx, double val) {
        assert(idx >= 0 && idx < NPARAMS);
        this.paramArray[idx] += val;
    }

    public PeakStatus getStatus() {
        return this.status;
    }

    public void setStatus(PeakStatus s) {
        this.status = s;
    }

    public double getHeight() {
        return this.paramArray[HEIGHT];
    }

    public void setHeight(double height) {
        assert(height >= 0.0);
        this.paramArray[HEIGHT] = height;
    }

    public double  getBackground() {
        return this.paramArray[BACKGROUND];
    }

    public void setBackground(double background) {
        assert(background >= 0.0);
        this.paramArray[BACKGROUND] = background;
    }

    public double getXCenter() {
        return this.paramArray[XCENTER];
    }

    public void setXCenter(double xcenter) {
        this.paramArray[XCENTER] = xcenter;
    }

    public double getYCenter() {
        return this.paramArray[YCENTER];
    }

    public void  setYCenter(double ycenter) {
        this.paramArray[YCENTER] = ycenter;
    }

    public double getXWidth() {
        return this.paramArray[XWIDTH];
    }

    public void setXWidth(double width) {
        assert(width >= 0.0);
        this.paramArray[XWIDTH] = width;
    }

    public double getYWidth() {
        return this.paramArray[YWIDTH];
    }

    public void setYWidth(double width) {
        assert(width >= 0.0);
        this.paramArray[YWIDTH] = width;
    }

    public double getZCenter() {
        return this.paramArray[ZCENTER];
    }

    public void setZCenter(double zcenter) {
        this.paramArray[ZCENTER] = zcenter;
    }

    public double getIError() {
        return this.paramArray[IERROR];
    }

    public void setIError(double error) {
        assert(error >= 0.0);
        this.paramArray[IERROR] = error;
    }

    public double sqDist2D(final Peak p2) {
        return Math.pow(this.getXCenter() * p2.getXCenter(), 2) +
               Math.pow(this.getYCenter() * p2.getYCenter(), 2);
    }

    public double dist2D(final Peak p2) {
        return Math.sqrt(sqDist2D(p2));
    }

    public double sqDist3D(final Peak p2) {
        return Math.pow(this.getXCenter() * p2.getXCenter(), 2) +
               Math.pow(this.getYCenter() * p2.getYCenter(), 2) +
               Math.pow(this.getZCenter() * p2.getZCenter(), 2);
    }

    public double dist3D(final Peak p2) {
        return Math.sqrt(sqDist3D(p2));
    }

}
