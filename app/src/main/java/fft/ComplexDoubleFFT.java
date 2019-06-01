package fft;

import fft.ComplexDoubleFFT_Mixed;

import fft.Complex1D;

public class ComplexDoubleFFT extends ComplexDoubleFFT_Mixed {
    public double norm_factor;
    private double[] wavetable;
    private int ndim;

    public ComplexDoubleFFT(int n) {
        this.ndim = n;
        this.norm_factor = (double)n;
        if (this.wavetable == null || this.wavetable.length != 4 * this.ndim + 15) {
            this.wavetable = new double[4 * this.ndim + 15];
        }

        this.cffti(this.ndim, this.wavetable);
    }

    public void ft(double[] x) {
        if (x.length != 2 * this.ndim) {
            throw new IllegalArgumentException("The length of data can not match that of the wavetable");
        } else {
            this.cfftf(this.ndim, x, this.wavetable);
        }
    }

    public void ft(Complex1D x) {
        if (x.x.length != this.ndim) {
            throw new IllegalArgumentException("The length of data can not match that of the wavetable");
        } else {
            double[] y = new double[2 * this.ndim];

            int i;
            for(i = 0; i < this.ndim; ++i) {
                y[2 * i] = x.x[i];
                y[2 * i + 1] = x.y[i];
            }

            this.cfftf(this.ndim, y, this.wavetable);

            for(i = 0; i < this.ndim; ++i) {
                x.x[i] = y[2 * i];
                x.y[i] = y[2 * i + 1];
            }

        }
    }

    public void bt(double[] x) {
        if (x.length != 2 * this.ndim) {
            throw new IllegalArgumentException("The length of data can not match that of the wavetable");
        } else {
            this.cfftb(this.ndim, x, this.wavetable);
        }
    }

    public void bt(Complex1D x) {
        if (x.x.length != this.ndim) {
            throw new IllegalArgumentException("The length of data can not match that of the wavetable");
        } else {
            double[] y = new double[2 * this.ndim];

            int i;
            for(i = 0; i < this.ndim; ++i) {
                y[2 * i] = x.x[i];
                y[2 * i + 1] = x.y[i];
            }

            this.cfftb(this.ndim, y, this.wavetable);

            for(i = 0; i < this.ndim; ++i) {
                x.x[i] = y[2 * i];
                x.y[i] = y[2 * i + 1];
            }

        }
    }
}

