package fft;

public class RealDoubleFFT_Even extends RealDoubleFFT_Mixed {
    public double norm_factor;
    private double[] wavetable;
    private int ndim;

    public RealDoubleFFT_Even(int n) {
        this.ndim = n;
        this.norm_factor = (double)(2 * (n - 1));
        if (this.wavetable == null || this.wavetable.length != 3 * this.ndim + 15) {
            this.wavetable = new double[3 * this.ndim + 15];
        }

        this.costi(this.ndim, this.wavetable);
    }

    public void ft(double[] x) {
        this.cost(this.ndim, x, this.wavetable);
    }

    public void bt(double[] x) {
        this.cost(this.ndim, x, this.wavetable);
    }

    void cost(int n, double[] x, double[] wtable) {
        int nm1 = n - 1;
        int ns2 = n / 2;
        if (n - 2 >= 0) {
            if (n == 2) {
                double x1h = x[0] + x[1];
                x[1] = x[0] - x[1];
                x[0] = x1h;
            } else if (n == 3) {
                double x1p3 = x[0] + x[2];
                double tx2 = x[1] + x[1];
                x[1] = x[0] - x[2];
                x[0] = x1p3 + tx2;
                x[2] = x1p3 - tx2;
            } else {
                double c1 = x[0] - x[n - 1];
                x[0] += x[n - 1];

                for(int k = 1; k < ns2; ++k) {
                    int kc = nm1 - k;
                    double t1 = x[k] + x[kc];
                    double t2 = x[k] - x[kc];
                    c1 += wtable[kc] * t2;
                    t2 = wtable[k] * t2;
                    x[k] = t1 - t2;
                    x[kc] = t1 + t2;
                }

                int modn = n % 2;
                if (modn != 0) {
                    x[ns2] += x[ns2];
                }

                this.rfftf1(nm1, x, wtable, n);
                double xim2 = x[1];
                x[1] = c1;

                for(int i = 3; i < n; i += 2) {
                    double xi = x[i];
                    x[i] = x[i - 2] - x[i - 1];
                    x[i - 1] = xim2;
                    xim2 = xi;
                }

                if (modn != 0) {
                    x[n - 1] = xim2;
                }
            }

        }
    }

    void costi(int n, double[] wtable) {
        double pi = 3.141592653589793D;
        if (n > 3) {
            int ns2 = n / 2;
            double dt = 3.141592653589793D / (double)(n - 1);

            for(int k = 1; k < ns2; ++k) {
                int kc = n - k - 1;
                wtable[k] = 2.0D * Math.sin((double)k * dt);
                wtable[kc] = 2.0D * Math.cos((double)k * dt);
            }

            this.rffti1(n - 1, wtable, n);
        }
    }
}

