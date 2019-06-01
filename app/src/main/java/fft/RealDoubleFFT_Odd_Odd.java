package fft;

public class RealDoubleFFT_Odd_Odd extends RealDoubleFFT_Even_Odd {
    public RealDoubleFFT_Odd_Odd(int n) {
        super(n);
    }

    public void ft(double[] x) {
        this.sinqf(this.ndim, x, this.wavetable);
    }

    public void bt(double[] x) {
        this.sinqb(this.ndim, x, this.wavetable);
    }

    void sinqf(int n, double[] x, double[] wtable) {
        if (n != 1) {
            int ns2 = n / 2;

            int k;
            for(k = 0; k < ns2; ++k) {
                int kc = n - k - 1;
                double xhold = x[k];
                x[k] = x[kc];
                x[kc] = xhold;
            }

            this.cosqf(n, x, wtable);

            for(k = 1; k < n; k += 2) {
                x[k] = -x[k];
            }

        }
    }

    void sinqb(int n, double[] x, double[] wtable) {
        if (n <= 1) {
            x[0] *= 4.0D;
        } else {
            int ns2 = n / 2;

            int k;
            for(k = 1; k < n; k += 2) {
                x[k] = -x[k];
            }

            this.cosqb(n, x, wtable);

            for(k = 0; k < ns2; ++k) {
                int kc = n - k - 1;
                double xhold = x[k];
                x[k] = x[kc];
                x[kc] = xhold;
            }

        }
    }
}

