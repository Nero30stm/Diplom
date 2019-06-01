package fft;

public class RealDoubleFFT_Even_Odd extends RealDoubleFFT_Mixed {
    public double norm_factor;
    protected double[] wavetable;
    protected int ndim;

    public RealDoubleFFT_Even_Odd(int n) {
        this.ndim = n;
        this.norm_factor = (double)(4 * n);
        if (this.wavetable == null || this.wavetable.length != 3 * this.ndim + 15) {
            this.wavetable = new double[3 * this.ndim + 15];
        }

        this.cosqi(this.ndim, this.wavetable);
    }

    public void ft(double[] x) {
        this.cosqf(this.ndim, x, this.wavetable);
    }

    public void bt(double[] x) {
        this.cosqb(this.ndim, x, this.wavetable);
    }

    void cosqf1(int n, double[] x, double[] wtable) {
        int ns2 = (n + 1) / 2;
        int np2 = n + 2;

        int k;
        int kc;
        for(k = 1; k < ns2; ++k) {
            kc = n - k;
            wtable[k + n] = x[k] + x[kc];
            wtable[kc + n] = x[k] - x[kc];
        }

        int modn = n % 2;
        if (modn == 0) {
            wtable[ns2 + n] = x[ns2] + x[ns2];
        }

        for(k = 1; k < ns2; ++k) {
            kc = n - k;
            x[k] = wtable[k - 1] * wtable[kc + n] + wtable[kc - 1] * wtable[k + n];
            x[kc] = wtable[k - 1] * wtable[k + n] - wtable[kc - 1] * wtable[kc + n];
        }

        if (modn == 0) {
            x[ns2] = wtable[ns2 - 1] * wtable[ns2 + n];
        }

        this.rfftf1(n, x, wtable, n);

        for(int i = 2; i < n; i += 2) {
            double xim1 = x[i - 1] - x[i];
            x[i] += x[i - 1];
            x[i - 1] = xim1;
        }

    }

    void cosqb1(int n, double[] x, double[] wtable) {
        int ns2 = (n + 1) / 2;

        for(int i = 2; i < n; i += 2) {
            double xim1 = x[i - 1] + x[i];
            x[i] -= x[i - 1];
            x[i - 1] = xim1;
        }

        x[0] += x[0];
        int modn = n % 2;
        if (modn == 0) {
            x[n - 1] += x[n - 1];
        }

        this.rfftb1(n, x, wtable, n);

        int k;
        int kc;
        for(k = 1; k < ns2; ++k) {
            kc = n - k;
            wtable[k + n] = wtable[k - 1] * x[kc] + wtable[kc - 1] * x[k];
            wtable[kc + n] = wtable[k - 1] * x[k] - wtable[kc - 1] * x[kc];
        }

        if (modn == 0) {
            x[ns2] = wtable[ns2 - 1] * (x[ns2] + x[ns2]);
        }

        for(k = 1; k < ns2; ++k) {
            kc = n - k;
            x[k] = wtable[k + n] + wtable[kc + n];
            x[kc] = wtable[k + n] - wtable[kc + n];
        }

        x[0] += x[0];
    }

    void cosqf(int n, double[] x, double[] wtable) {
        double sqrt2 = 1.4142135623731D;
        if (n >= 2) {
            if (n == 2) {
                double tsqx = 1.4142135623731D * x[1];
                x[1] = x[0] - tsqx;
                x[0] += tsqx;
            } else {
                this.cosqf1(n, x, wtable);
            }

        }
    }

    void cosqb(int n, double[] x, double[] wtable) {
        double tsqrt2 = 2.82842712474619D;
        if (n < 2) {
            x[0] *= 4.0D;
        } else if (n == 2) {
            double x1 = 4.0D * (x[0] + x[1]);
            x[1] = 2.82842712474619D * (x[0] - x[1]);
            x[0] = x1;
        } else {
            this.cosqb1(n, x, wtable);
        }

    }

    void cosqi(int n, double[] wtable) {
        double pih = 1.5707963267948966D;
        double dt = 1.5707963267948966D / (double)n;

        for(int k = 0; k < n; ++k) {
            wtable[k] = Math.cos((double)(k + 1) * dt);
        }

        this.rffti1(n, wtable, n);
    }
}

