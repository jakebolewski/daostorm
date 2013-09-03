package edu.mcmaster.daostorm;

public class ClampParameters {

    private final double height;
    private final double background;
    private final double xcenter;
    private final double ycenter;
    private final double xwidth;
    private final double ywidth;
    private final double zcenter;

    // 2D Constructor
    ClampParameters(final double height,
                    final double background,
                    final double xcenter,
                    final double ycenter,
                    final double xwidth,
                    final double ywidth)
    {
        assert(height >= 0 &&
               background >= 0 &&
               xwidth >= 0 &&
               ywidth >= 0);
        this.height = height;
        this.background = background;
        this.xcenter = xcenter;
        this.ycenter = ycenter;
        this.xwidth = xwidth;
        this.ywidth = ywidth;
        this.zcenter = 0.0;
    }

    // 3D Constructor
    ClampParameters(final double height,
                    final double background,
                    final double xcenter,
                    final double ycenter,
                    final double xwidth,
                    final double ywidth,
                    final double zcenter)
    {
        assert(height >= 0 &&
               background >= 0 &&
               xwidth >= 0 &&
               ywidth >= 0);

        this.height = height;
        this.background = background;
        this.xcenter = xcenter;
        this.ycenter = ycenter;
        this.xwidth = xwidth;
        this.ywidth = ywidth;
        this.zcenter = zcenter;
    }

    double height() {
        return this.height;
    }

    double background() {
        return this.background;
    }

    double xcenter() {
        return this.xcenter;
    }

    double ycenter() {
        return this.ycenter;
    }

    double xwidth() {
        return this.xwidth;
    }

    double ywidth() {
        return this.ywidth;
    }

    double zcenter() {
        return this.zcenter;
    }
};
