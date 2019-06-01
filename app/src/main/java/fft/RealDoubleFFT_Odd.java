package fft;

public class RealDoubleFFT_Odd extends RealDoubleFFT_Mixed {
    public double norm_factor;
    private double[] wavetable;
    private int ndim;

    public RealDoubleFFT_Odd(int n) {
        this.ndim = n;
        this.norm_factor = (double)(2 * (n + 1));
        int wtable_length = 2 * this.ndim + this.ndim / 2 + 3 + 15;
        if (this.wavetable == null || this.wavetable.length != wtable_length) {
            this.wavetable = new double[wtable_length];
        }

        this.sinti(this.ndim, this.wavetable);
    }

    public void ft(double[] x) {
        this.sint(this.ndim, x, this.wavetable);
    }

    public void bt(double[] x) {
        this.sint(this.ndim, x, this.wavetable);
    }

    void sint1(int n, double[] war, double[] wtable) {
        double sqrt3 = 1.73205080756888D;
        double[] wtable_p1 = new double[2 * (n + 1) + 15];
        int iw1 = n / 2;
        int iw2 = iw1 + n + 1;
        int iw3 = iw2 + n + 1;
        double[] x = new double[n + 1];

        int i;
        for(i = 0; i < n; ++i) {
            wtable[i + iw1] = war[i];
            war[i] = wtable[i + iw2];
        }

        if (n < 2) {
            wtable[0 + iw1] += wtable[0 + iw1];
        } else if (n == 2) {
            double xhold = 1.73205080756888D * (wtable[0 + iw1] + wtable[1 + iw1]);
            wtable[1 + iw1] = 1.73205080756888D * (wtable[0 + iw1] - wtable[1 + iw1]);
            wtable[0 + iw1] = xhold;
        } else {
            int np1 = n + 1;
            int ns2 = n / 2;
            wtable[0 + iw2] = 0.0D;

            for(int k = 0; k < ns2; ++k) {
                int kc = n - k - 1;
                double t1 = wtable[k + iw1] - wtable[kc + iw1];
                double t2 = wtable[k] * (wtable[k + iw1] + wtable[kc + iw1]);
                wtable[k + 1 + iw2] = t1 + t2;
                wtable[kc + 1 + iw2] = t2 - t1;
            }

            int modn = n % 2;
            if (modn != 0) {
                wtable[ns2 + 1 + iw2] = 4.0D * wtable[ns2 + iw1];
            }

            System.arraycopy(wtable, iw1, wtable_p1, 0, n + 1);
            System.arraycopy(war, 0, wtable_p1, n + 1, n);
            System.arraycopy(wtable, iw3, wtable_p1, 2 * (n + 1), 15);
            System.arraycopy(wtable, iw2, x, 0, n + 1);
            this.rfftf1(np1, x, wtable_p1, 0);
            System.arraycopy(x, 0, wtable, iw2, n + 1);
            wtable[0 + iw1] = 0.5D * wtable[0 + iw2];

            for(i = 2; i < n; i += 2) {
                wtable[i - 1 + iw1] = -wtable[i + iw2];
                wtable[i + iw1] = wtable[i - 2 + iw1] + wtable[i - 1 + iw2];
            }

            if (modn == 0) {
                wtable[n - 1 + iw1] = -wtable[n + iw2];
            }
        }

        for(i = 0; i < n; ++i) {
            wtable[i + iw2] = war[i];
            war[i] = wtable[i + iw1];
        }

    }

    void sint(int n, double[] x, double[] wtable) {
        this.sint1(n, x, wtable);
    }

    void sinti(int n, double[] wtable) {
        double pi = 3.141592653589793D;
        if (n > 1) {
            int ns2 = n / 2;
            double dt = 3.141592653589793D / (double)(n + 1);

            for(int k = 0; k < ns2; ++k) {
                wtable[k] = 2.0D * Math.sin((double)(k + 1) * dt);
            }

            this.rffti1(n + 1, wtable, ns2);
        }
    }
}
