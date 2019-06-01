package fft;

public class RealDoubleFFT extends RealDoubleFFT_Mixed {
    public double norm_factor;
    private double[] wavetable;
    private int ndim;

    public RealDoubleFFT(int n) {
        this.ndim = n;
        this.norm_factor = (double)n;
        if (this.wavetable == null || this.wavetable.length != 2 * this.ndim + 15) {
            this.wavetable = new double[2 * this.ndim + 15];
        }

        this.rffti(this.ndim, this.wavetable);
    }

    public void ft(double[] x) {
        if (x.length != this.ndim) {
            throw new IllegalArgumentException("The length of data can not match that of the wavetable");
        } else {
            this.rfftf(this.ndim, x, this.wavetable);
        }
    }

    public void ft(double[] x, Complex1D y) {
        if (x.length != this.ndim) {
            throw new IllegalArgumentException("The length of data can not match that of the wavetable");
        } else {
            this.rfftf(this.ndim, x, this.wavetable);
            if (this.ndim % 2 == 0) {
                y.x = new double[this.ndim / 2 + 1];
                y.y = new double[this.ndim / 2 + 1];
            } else {
                y.x = new double[(this.ndim + 1) / 2];
                y.y = new double[(this.ndim + 1) / 2];
            }

            y.x[0] = x[0];
            y.y[0] = 0.0D;

            for(int i = 1; i < (this.ndim + 1) / 2; ++i) {
                y.x[i] = x[2 * i - 1];
                y.y[i] = x[2 * i];
            }

            if (this.ndim % 2 == 0) {
                y.x[this.ndim / 2] = x[this.ndim - 1];
                y.y[this.ndim / 2] = 0.0D;
            }

        }
    }

    public void bt(double[] x) {
        if (x.length != this.ndim) {
            throw new IllegalArgumentException("The length of data can not match that of the wavetable");
        } else {
            this.rfftb(this.ndim, x, this.wavetable);
        }
    }

    public void bt(Complex1D x, double[] y) {
        if (this.ndim % 2 == 0) {
            if (x.x.length != this.ndim / 2 + 1) {
                throw new IllegalArgumentException("The length of data can not match that of the wavetable");
            }
        } else if (x.x.length != (this.ndim + 1) / 2) {
            throw new IllegalArgumentException("The length of data can not match that of the wavetable");
        }

        y[0] = x.x[0];

        for(int i = 1; i < (this.ndim + 1) / 2; ++i) {
            y[2 * i - 1] = x.x[i];
            y[2 * i] = x.y[i];
        }

        if (this.ndim % 2 == 0) {
            y[this.ndim - 1] = x.x[this.ndim / 2];
        }

        this.rfftb(this.ndim, y, this.wavetable);
    }
}
