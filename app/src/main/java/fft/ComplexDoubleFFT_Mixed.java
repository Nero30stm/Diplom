package fft;

class ComplexDoubleFFT_Mixed {
    ComplexDoubleFFT_Mixed() {
    }

    void passf2(int ido, int l1, double[] cc, double[] ch, double[] wtable, int offset, int isign) {
        int iw1 = offset;
        int k;
        int ah;
        int ac;
        if (ido <= 2) {
            for(k = 0; k < l1; ++k) {
                ah = k * ido;
                ac = 2 * k * ido;
                ch[ah] = cc[ac] + cc[ac + ido];
                ch[ah + ido * l1] = cc[ac] - cc[ac + ido];
                ch[ah + 1] = cc[ac + 1] + cc[ac + ido + 1];
                ch[ah + ido * l1 + 1] = cc[ac + 1] - cc[ac + ido + 1];
            }
        } else {
            for(k = 0; k < l1; ++k) {
                for(int i = 0; i < ido - 1; i += 2) {
                    ah = i + k * ido;
                    ac = i + 2 * k * ido;
                    ch[ah] = cc[ac] + cc[ac + ido];
                    double tr2 = cc[ac] - cc[ac + ido];
                    ch[ah + 1] = cc[ac + 1] + cc[ac + 1 + ido];
                    double ti2 = cc[ac + 1] - cc[ac + 1 + ido];
                    ch[ah + l1 * ido + 1] = wtable[i + iw1] * ti2 + (double)isign * wtable[i + 1 + iw1] * tr2;
                    ch[ah + l1 * ido] = wtable[i + iw1] * tr2 - (double)isign * wtable[i + 1 + iw1] * ti2;
                }
            }
        }

    }

    void passf3(int ido, int l1, double[] cc, double[] ch, double[] wtable, int offset, int isign) {
        double taur = -0.5D;
        double taui = 0.866025403784439D;
        int iw1 = offset;
        int iw2 = offset + ido;
        int k;
        int ac;
        int ah;
        double ci2;
        double ci3;
        double cr2;
        double cr3;
        double ti2;
        double tr2;
        if (ido == 2) {
            for(k = 1; k <= l1; ++k) {
                ac = (3 * k - 2) * ido;
                tr2 = cc[ac] + cc[ac + ido];
                cr2 = cc[ac - ido] + -0.5D * tr2;
                ah = (k - 1) * ido;
                ch[ah] = cc[ac - ido] + tr2;
                ti2 = cc[ac + 1] + cc[ac + ido + 1];
                ci2 = cc[ac - ido + 1] + -0.5D * ti2;
                ch[ah + 1] = cc[ac - ido + 1] + ti2;
                cr3 = (double)isign * 0.866025403784439D * (cc[ac] - cc[ac + ido]);
                ci3 = (double)isign * 0.866025403784439D * (cc[ac + 1] - cc[ac + ido + 1]);
                ch[ah + l1 * ido] = cr2 - ci3;
                ch[ah + 2 * l1 * ido] = cr2 + ci3;
                ch[ah + l1 * ido + 1] = ci2 + cr3;
                ch[ah + 2 * l1 * ido + 1] = ci2 - cr3;
            }
        } else {
            for(k = 1; k <= l1; ++k) {
                for(int i = 0; i < ido - 1; i += 2) {
                    ac = i + (3 * k - 2) * ido;
                    tr2 = cc[ac] + cc[ac + ido];
                    cr2 = cc[ac - ido] + -0.5D * tr2;
                    ah = i + (k - 1) * ido;
                    ch[ah] = cc[ac - ido] + tr2;
                    ti2 = cc[ac + 1] + cc[ac + ido + 1];
                    ci2 = cc[ac - ido + 1] + -0.5D * ti2;
                    ch[ah + 1] = cc[ac - ido + 1] + ti2;
                    cr3 = (double)isign * 0.866025403784439D * (cc[ac] - cc[ac + ido]);
                    ci3 = (double)isign * 0.866025403784439D * (cc[ac + 1] - cc[ac + ido + 1]);
                    double dr2 = cr2 - ci3;
                    double dr3 = cr2 + ci3;
                    double di2 = ci2 + cr3;
                    double di3 = ci2 - cr3;
                    ch[ah + l1 * ido + 1] = wtable[i + iw1] * di2 + (double)isign * wtable[i + 1 + iw1] * dr2;
                    ch[ah + l1 * ido] = wtable[i + iw1] * dr2 - (double)isign * wtable[i + 1 + iw1] * di2;
                    ch[ah + 2 * l1 * ido + 1] = wtable[i + iw2] * di3 + (double)isign * wtable[i + 1 + iw2] * dr3;
                    ch[ah + 2 * l1 * ido] = wtable[i + iw2] * dr3 - (double)isign * wtable[i + 1 + iw2] * di3;
                }
            }
        }

    }

    void passf4(int ido, int l1, double[] cc, double[] ch, double[] wtable, int offset, int isign) {
        int iw1 = offset;
        int iw2 = offset + ido;
        int iw3 = iw2 + ido;
        int k;
        int ac;
        int ah;
        double ti1;
        double ti2;
        double ti3;
        double ti4;
        double tr1;
        double tr2;
        double tr3;
        double tr4;
        if (ido == 2) {
            for(k = 0; k < l1; ++k) {
                ac = 4 * k * ido + 1;
                ti1 = cc[ac] - cc[ac + 2 * ido];
                ti2 = cc[ac] + cc[ac + 2 * ido];
                tr4 = cc[ac + 3 * ido] - cc[ac + ido];
                ti3 = cc[ac + ido] + cc[ac + 3 * ido];
                tr1 = cc[ac - 1] - cc[ac + 2 * ido - 1];
                tr2 = cc[ac - 1] + cc[ac + 2 * ido - 1];
                ti4 = cc[ac + ido - 1] - cc[ac + 3 * ido - 1];
                tr3 = cc[ac + ido - 1] + cc[ac + 3 * ido - 1];
                ah = k * ido;
                ch[ah] = tr2 + tr3;
                ch[ah + 2 * l1 * ido] = tr2 - tr3;
                ch[ah + 1] = ti2 + ti3;
                ch[ah + 2 * l1 * ido + 1] = ti2 - ti3;
                ch[ah + l1 * ido] = tr1 + (double)isign * tr4;
                ch[ah + 3 * l1 * ido] = tr1 - (double)isign * tr4;
                ch[ah + l1 * ido + 1] = ti1 + (double)isign * ti4;
                ch[ah + 3 * l1 * ido + 1] = ti1 - (double)isign * ti4;
            }
        } else {
            for(k = 0; k < l1; ++k) {
                for(int i = 0; i < ido - 1; i += 2) {
                    ac = i + 1 + 4 * k * ido;
                    ti1 = cc[ac] - cc[ac + 2 * ido];
                    ti2 = cc[ac] + cc[ac + 2 * ido];
                    ti3 = cc[ac + ido] + cc[ac + 3 * ido];
                    tr4 = cc[ac + 3 * ido] - cc[ac + ido];
                    tr1 = cc[ac - 1] - cc[ac + 2 * ido - 1];
                    tr2 = cc[ac - 1] + cc[ac + 2 * ido - 1];
                    ti4 = cc[ac + ido - 1] - cc[ac + 3 * ido - 1];
                    tr3 = cc[ac + ido - 1] + cc[ac + 3 * ido - 1];
                    ah = i + k * ido;
                    ch[ah] = tr2 + tr3;
                    double cr3 = tr2 - tr3;
                    ch[ah + 1] = ti2 + ti3;
                    double ci3 = ti2 - ti3;
                    double cr2 = tr1 + (double)isign * tr4;
                    double cr4 = tr1 - (double)isign * tr4;
                    double ci2 = ti1 + (double)isign * ti4;
                    double ci4 = ti1 - (double)isign * ti4;
                    ch[ah + l1 * ido] = wtable[i + iw1] * cr2 - (double)isign * wtable[i + 1 + iw1] * ci2;
                    ch[ah + l1 * ido + 1] = wtable[i + iw1] * ci2 + (double)isign * wtable[i + 1 + iw1] * cr2;
                    ch[ah + 2 * l1 * ido] = wtable[i + iw2] * cr3 - (double)isign * wtable[i + 1 + iw2] * ci3;
                    ch[ah + 2 * l1 * ido + 1] = wtable[i + iw2] * ci3 + (double)isign * wtable[i + 1 + iw2] * cr3;
                    ch[ah + 3 * l1 * ido] = wtable[i + iw3] * cr4 - (double)isign * wtable[i + 1 + iw3] * ci4;
                    ch[ah + 3 * l1 * ido + 1] = wtable[i + iw3] * ci4 + (double)isign * wtable[i + 1 + iw3] * cr4;
                }
            }
        }

    }

    void passf5(int ido, int l1, double[] cc, double[] ch, double[] wtable, int offset, int isign) {
        double tr11 = 0.309016994374947D;
        double ti11 = 0.951056516295154D;
        double tr12 = -0.809016994374947D;
        double ti12 = 0.587785252292473D;
        int iw1 = offset;
        int iw2 = offset + ido;
        int iw3 = iw2 + ido;
        int iw4 = iw3 + ido;
        double tr4;
        double tr5;
        int k;
        int ac;
        int ah;
        double ci2;
        double ci3;
        double ci4;
        double ci5;
        double cr2;
        double cr3;
        double cr5;
        double cr4;
        double ti2;
        double ti3;
        double ti4;
        double ti5;
        double tr2;
        double tr3;
        if (ido == 2) {
            for(k = 1; k <= l1; ++k) {
                ac = (5 * k - 4) * ido + 1;
                ti5 = cc[ac] - cc[ac + 3 * ido];
                ti2 = cc[ac] + cc[ac + 3 * ido];
                ti4 = cc[ac + ido] - cc[ac + 2 * ido];
                ti3 = cc[ac + ido] + cc[ac + 2 * ido];
                tr5 = cc[ac - 1] - cc[ac + 3 * ido - 1];
                tr2 = cc[ac - 1] + cc[ac + 3 * ido - 1];
                tr4 = cc[ac + ido - 1] - cc[ac + 2 * ido - 1];
                tr3 = cc[ac + ido - 1] + cc[ac + 2 * ido - 1];
                ah = (k - 1) * ido;
                ch[ah] = cc[ac - ido - 1] + tr2 + tr3;
                ch[ah + 1] = cc[ac - ido] + ti2 + ti3;
                cr2 = cc[ac - ido - 1] + 0.309016994374947D * tr2 + -0.809016994374947D * tr3;
                ci2 = cc[ac - ido] + 0.309016994374947D * ti2 + -0.809016994374947D * ti3;
                cr3 = cc[ac - ido - 1] + -0.809016994374947D * tr2 + 0.309016994374947D * tr3;
                ci3 = cc[ac - ido] + -0.809016994374947D * ti2 + 0.309016994374947D * ti3;
                cr5 = (double)isign * (0.951056516295154D * tr5 + 0.587785252292473D * tr4);
                ci5 = (double)isign * (0.951056516295154D * ti5 + 0.587785252292473D * ti4);
                cr4 = (double)isign * (0.587785252292473D * tr5 - 0.951056516295154D * tr4);
                ci4 = (double)isign * (0.587785252292473D * ti5 - 0.951056516295154D * ti4);
                ch[ah + l1 * ido] = cr2 - ci5;
                ch[ah + 4 * l1 * ido] = cr2 + ci5;
                ch[ah + l1 * ido + 1] = ci2 + cr5;
                ch[ah + 2 * l1 * ido + 1] = ci3 + cr4;
                ch[ah + 2 * l1 * ido] = cr3 - ci4;
                ch[ah + 3 * l1 * ido] = cr3 + ci4;
                ch[ah + 3 * l1 * ido + 1] = ci3 - cr4;
                ch[ah + 4 * l1 * ido + 1] = ci2 - cr5;
            }
        } else {
            for(k = 1; k <= l1; ++k) {
                for(int i = 0; i < ido - 1; i += 2) {
                    ac = i + 1 + (k * 5 - 4) * ido;
                    ti5 = cc[ac] - cc[ac + 3 * ido];
                    ti2 = cc[ac] + cc[ac + 3 * ido];
                    ti4 = cc[ac + ido] - cc[ac + 2 * ido];
                    ti3 = cc[ac + ido] + cc[ac + 2 * ido];
                    tr5 = cc[ac - 1] - cc[ac + 3 * ido - 1];
                    tr2 = cc[ac - 1] + cc[ac + 3 * ido - 1];
                    tr4 = cc[ac + ido - 1] - cc[ac + 2 * ido - 1];
                    tr3 = cc[ac + ido - 1] + cc[ac + 2 * ido - 1];
                    ah = i + (k - 1) * ido;
                    ch[ah] = cc[ac - ido - 1] + tr2 + tr3;
                    ch[ah + 1] = cc[ac - ido] + ti2 + ti3;
                    cr2 = cc[ac - ido - 1] + 0.309016994374947D * tr2 + -0.809016994374947D * tr3;
                    ci2 = cc[ac - ido] + 0.309016994374947D * ti2 + -0.809016994374947D * ti3;
                    cr3 = cc[ac - ido - 1] + -0.809016994374947D * tr2 + 0.309016994374947D * tr3;
                    ci3 = cc[ac - ido] + -0.809016994374947D * ti2 + 0.309016994374947D * ti3;
                    cr5 = (double)isign * (0.951056516295154D * tr5 + 0.587785252292473D * tr4);
                    ci5 = (double)isign * (0.951056516295154D * ti5 + 0.587785252292473D * ti4);
                    cr4 = (double)isign * (0.587785252292473D * tr5 - 0.951056516295154D * tr4);
                    ci4 = (double)isign * (0.587785252292473D * ti5 - 0.951056516295154D * ti4);
                    double dr3 = cr3 - ci4;
                    double dr4 = cr3 + ci4;
                    double di3 = ci3 + cr4;
                    double di4 = ci3 - cr4;
                    double dr5 = cr2 + ci5;
                    double dr2 = cr2 - ci5;
                    double di5 = ci2 - cr5;
                    double di2 = ci2 + cr5;
                    ch[ah + l1 * ido] = wtable[i + iw1] * dr2 - (double)isign * wtable[i + 1 + iw1] * di2;
                    ch[ah + l1 * ido + 1] = wtable[i + iw1] * di2 + (double)isign * wtable[i + 1 + iw1] * dr2;
                    ch[ah + 2 * l1 * ido] = wtable[i + iw2] * dr3 - (double)isign * wtable[i + 1 + iw2] * di3;
                    ch[ah + 2 * l1 * ido + 1] = wtable[i + iw2] * di3 + (double)isign * wtable[i + 1 + iw2] * dr3;
                    ch[ah + 3 * l1 * ido] = wtable[i + iw3] * dr4 - (double)isign * wtable[i + 1 + iw3] * di4;
                    ch[ah + 3 * l1 * ido + 1] = wtable[i + iw3] * di4 + (double)isign * wtable[i + 1 + iw3] * dr4;
                    ch[ah + 4 * l1 * ido] = wtable[i + iw4] * dr5 - (double)isign * wtable[i + 1 + iw4] * di5;
                    ch[ah + 4 * l1 * ido + 1] = wtable[i + iw4] * di5 + (double)isign * wtable[i + 1 + iw4] * dr5;
                }
            }
        }

    }

    void passfg(int[] nac, int ido, int ip, int l1, int idl1, double[] cc, double[] c1, double[] c2, double[] ch, double[] ch2, double[] wtable, int offset, int isign) {
        int iw1 = offset;
        int idot = ido / 2;
        int var10000 = ip * idl1;
        int ipph = (ip + 1) / 2;
        int idp = ip * ido;
        int i;
        int j;
        int k;
        int jc;
        if (ido >= l1) {
            for(j = 1; j < ipph; ++j) {
                jc = ip - j;

                for(k = 0; k < l1; ++k) {
                    for(i = 0; i < ido; ++i) {
                        ch[i + (k + j * l1) * ido] = cc[i + (j + k * ip) * ido] + cc[i + (jc + k * ip) * ido];
                        ch[i + (k + jc * l1) * ido] = cc[i + (j + k * ip) * ido] - cc[i + (jc + k * ip) * ido];
                    }
                }
            }

            for(k = 0; k < l1; ++k) {
                for(i = 0; i < ido; ++i) {
                    ch[i + k * ido] = cc[i + k * ip * ido];
                }
            }
        } else {
            for(j = 1; j < ipph; ++j) {
                jc = ip - j;

                for(i = 0; i < ido; ++i) {
                    for(k = 0; k < l1; ++k) {
                        ch[i + (k + j * l1) * ido] = cc[i + (j + k * ip) * ido] + cc[i + (jc + k * ip) * ido];
                        ch[i + (k + jc * l1) * ido] = cc[i + (j + k * ip) * ido] - cc[i + (jc + k * ip) * ido];
                    }
                }
            }

            for(i = 0; i < ido; ++i) {
                for(k = 0; k < l1; ++k) {
                    ch[i + k * ido] = cc[i + k * ip * ido];
                }
            }
        }

        int idl = 2 - ido;
        int inc = 0;

        int ik;
        for(int l = 1; l < ipph; ++l) {
            int lc = ip - l;
            idl += ido;

            for(ik = 0; ik < idl1; ++ik) {
                c2[ik + l * idl1] = ch2[ik] + wtable[idl - 2 + iw1] * ch2[ik + idl1];
                c2[ik + lc * idl1] = (double)isign * wtable[idl - 1 + iw1] * ch2[ik + (ip - 1) * idl1];
            }

            int idlj = idl;
            inc += ido;

            for(j = 2; j < ipph; ++j) {
                jc = ip - j;
                idlj += inc;
                if (idlj > idp) {
                    idlj -= idp;
                }

                double war = wtable[idlj - 2 + iw1];
                double wai = wtable[idlj - 1 + iw1];

                for(ik = 0; ik < idl1; ++ik) {
                    c2[ik + l * idl1] += war * ch2[ik + j * idl1];
                    c2[ik + lc * idl1] += (double)isign * wai * ch2[ik + jc * idl1];
                }
            }
        }

        for(j = 1; j < ipph; ++j) {
            for(ik = 0; ik < idl1; ++ik) {
                ch2[ik] += ch2[ik + j * idl1];
            }
        }

        for(j = 1; j < ipph; ++j) {
            jc = ip - j;

            for(ik = 1; ik < idl1; ik += 2) {
                ch2[ik - 1 + j * idl1] = c2[ik - 1 + j * idl1] - c2[ik + jc * idl1];
                ch2[ik - 1 + jc * idl1] = c2[ik - 1 + j * idl1] + c2[ik + jc * idl1];
                ch2[ik + j * idl1] = c2[ik + j * idl1] + c2[ik - 1 + jc * idl1];
                ch2[ik + jc * idl1] = c2[ik + j * idl1] - c2[ik - 1 + jc * idl1];
            }
        }

        nac[0] = 1;
        if (ido != 2) {
            nac[0] = 0;

            for(ik = 0; ik < idl1; ++ik) {
                c2[ik] = ch2[ik];
            }

            for(j = 1; j < ip; ++j) {
                for(k = 0; k < l1; ++k) {
                    c1[(k + j * l1) * ido + 0] = ch[(k + j * l1) * ido + 0];
                    c1[(k + j * l1) * ido + 1] = ch[(k + j * l1) * ido + 1];
                }
            }

            int idij;
            if (idot <= l1) {
                idij = 0;

                for(j = 1; j < ip; ++j) {
                    idij += 2;

                    for(i = 3; i < ido; i += 2) {
                        idij += 2;

                        for(k = 0; k < l1; ++k) {
                            c1[i - 1 + (k + j * l1) * ido] = wtable[idij - 2 + iw1] * ch[i - 1 + (k + j * l1) * ido] - (double)isign * wtable[idij - 1 + iw1] * ch[i + (k + j * l1) * ido];
                            c1[i + (k + j * l1) * ido] = wtable[idij - 2 + iw1] * ch[i + (k + j * l1) * ido] + (double)isign * wtable[idij - 1 + iw1] * ch[i - 1 + (k + j * l1) * ido];
                        }
                    }
                }
            } else {
                int idj = 2 - ido;

                for(j = 1; j < ip; ++j) {
                    idj += ido;

                    for(k = 0; k < l1; ++k) {
                        idij = idj;

                        for(i = 3; i < ido; i += 2) {
                            idij += 2;
                            c1[i - 1 + (k + j * l1) * ido] = wtable[idij - 2 + iw1] * ch[i - 1 + (k + j * l1) * ido] - (double)isign * wtable[idij - 1 + iw1] * ch[i + (k + j * l1) * ido];
                            c1[i + (k + j * l1) * ido] = wtable[idij - 2 + iw1] * ch[i + (k + j * l1) * ido] + (double)isign * wtable[idij - 1 + iw1] * ch[i - 1 + (k + j * l1) * ido];
                        }
                    }
                }
            }

        }
    }

    void cfftf1(int n, double[] c, double[] wtable, int isign) {
        int[] nac = new int[1];
        double[] ch = new double[2 * n];
        int iw1 = 2 * n;
        int iw2 = 4 * n;
        System.arraycopy(wtable, 0, ch, 0, 2 * n);
        nac[0] = 0;
        int nf = (int)wtable[1 + iw2];
        int na = 0;
        int l1 = 1;
        int iw = iw1;

        for(int k1 = 2; k1 <= nf + 1; ++k1) {
            int ip = (int)wtable[k1 + iw2];
            int l2 = ip * l1;
            int ido = n / l2;
            int idot = ido + ido;
            int idl1 = idot * l1;
            if (ip == 4) {
                if (na == 0) {
                    this.passf4(idot, l1, c, ch, wtable, iw, isign);
                } else {
                    this.passf4(idot, l1, ch, c, wtable, iw, isign);
                }

                na = 1 - na;
            } else if (ip == 2) {
                if (na == 0) {
                    this.passf2(idot, l1, c, ch, wtable, iw, isign);
                } else {
                    this.passf2(idot, l1, ch, c, wtable, iw, isign);
                }

                na = 1 - na;
            } else if (ip == 3) {
                if (na == 0) {
                    this.passf3(idot, l1, c, ch, wtable, iw, isign);
                } else {
                    this.passf3(idot, l1, ch, c, wtable, iw, isign);
                }

                na = 1 - na;
            } else if (ip == 5) {
                if (na == 0) {
                    this.passf5(idot, l1, c, ch, wtable, iw, isign);
                } else {
                    this.passf5(idot, l1, ch, c, wtable, iw, isign);
                }

                na = 1 - na;
            } else {
                if (na == 0) {
                    this.passfg(nac, idot, ip, l1, idl1, c, c, c, ch, ch, wtable, iw, isign);
                } else {
                    this.passfg(nac, idot, ip, l1, idl1, ch, ch, ch, c, c, wtable, iw, isign);
                }

                if (nac[0] != 0) {
                    na = 1 - na;
                }
            }

            l1 = l2;
            iw += (ip - 1) * idot;
        }

        if (na != 0) {
            for(int i = 0; i < 2 * n; ++i) {
                c[i] = ch[i];
            }

        }
    }

    void cfftf(int n, double[] c, double[] wtable) {
        this.cfftf1(n, c, wtable, -1);
    }

    void cfftb(int n, double[] c, double[] wtable) {
        this.cfftf1(n, c, wtable, 1);
    }

    void cffti1(int n, double[] wtable) {
        int[] ntryh = new int[]{3, 4, 2, 5};
        double twopi = 6.283185307179586D;
        int ntry = 0;
        int nl = n;
        int nf = 0;
        int j = 0;

        label66:
        while(true) {
            ++j;
            if (j <= 4) {
                ntry = ntryh[j - 1];
            } else {
                ntry += 2;
            }

            int i;
            int nq;
            do {
                nq = nl / ntry;
                int nr = nl - ntry * nq;
                if (nr != 0) {
                    continue label66;
                }

                ++nf;
                wtable[nf + 1 + 4 * n] = (double)ntry;
                nl = nq;
                if (ntry == 2 && nf != 1) {
                    for(i = 2; i <= nf; ++i) {
                        int ib = nf - i + 2;
                        wtable[ib + 1 + 4 * n] = wtable[ib + 4 * n];
                    }

                    wtable[2 + 4 * n] = 2.0D;
                }
            } while(nq != 1);

            wtable[0 + 4 * n] = (double)n;
            wtable[1 + 4 * n] = (double)nf;
            double argh = 6.283185307179586D / (double)n;
            i = 1;
            int l1 = 1;

            for(int k1 = 1; k1 <= nf; ++k1) {
                int ip = (int)wtable[k1 + 1 + 4 * n];
                int ld = 0;
                int l2 = l1 * ip;
                int ido = n / l2;
                int idot = ido + ido + 2;
                int ipm = ip - 1;

                for(j = 1; j <= ipm; ++j) {
                    int i1 = i;
                    wtable[i - 1 + 2 * n] = 1.0D;
                    wtable[i + 2 * n] = 0.0D;
                    ld += l1;
                    double fi = 0.0D;
                    double argld = (double)ld * argh;

                    for(int ii = 4; ii <= idot; ii += 2) {
                        i += 2;
                        ++fi;
                        double arg = fi * argld;
                        wtable[i - 1 + 2 * n] = Math.cos(arg);
                        wtable[i + 2 * n] = Math.sin(arg);
                    }

                    if (ip > 5) {
                        wtable[i - 1 + 2 * n] = wtable[i - 1 + 2 * n];
                        wtable[i1 + 2 * n] = wtable[i + 2 * n];
                    }
                }

                l1 = l2;
            }

            return;
        }
    }

    void cffti(int n, double[] wtable) {
        if (n != 1) {
            this.cffti1(n, wtable);
        }
    }
}