#include "fcwt.h"

Morlet::Morlet(float bandwidth) {
    four_wavelen = 0.9876f;
    fb = bandwidth;
    fb2 = 2.0f * fb * fb;
    ifb = 1.0f / fb;
    imag_frequency = false;
    doublesided = false;
    mother = nullptr;
}

void Morlet::generate(int size) {
    // Frequency domain, because we only need size. Default scale is always 2;
    width = size;

    float tmp1;
    float toradians = (2 * PI) / static_cast<float>(size);
    float norm = sqrt(2 * PI) * IPI4;

    mother = new float[width];

    // Calculate array
    for (int w = 0; w < width; w++) {
        tmp1 = (2.0f * static_cast<float>(w) * toradians * fb - 2.0f * PI * fb);
        tmp1 = -(tmp1 * tmp1) / 2;
        mother[w] = (norm * exp(tmp1));
    }
}

void Morlet::generate(float* real, float* imag, int size, float scale) {
    // Time domain because we know size from scale
    float tmp1, tmp2;
    width = getSupport(scale);
    float norm = static_cast<float>(size) * ifb * IPI4;

    for (int t = 0; t < width * 2 + 1; t++) {
        tmp1 = static_cast<float>(t - width) / scale;
        tmp2 = exp(-(tmp1 * tmp1) / fb2);

        real[t] = norm * tmp2 * cos(tmp1 * 2.0f * PI) / scale;
        imag[t] = norm * tmp2 * sin(tmp1 * 2.0f * PI) / scale;
    }
}

void Morlet::getWavelet(float scale, complex<float>* pwav, int pn) {
    int w = getSupport(scale);

    float* real = new float[max(w * 2 + 1, pn)];
    float* imag = new float[max(w * 2 + 1, pn)];
    for (int t = 0; t < max(w * 2 + 1, pn); t++) {
        real[t] = 0;
        imag[t] = 0;
    }

    generate(real, imag, pn, scale);

    for (int t = 0; t < pn; t++) {
        pwav[t].real(real[t]);
        pwav[t].imag(imag[t]);
    }

    delete[] real;
    delete[] imag;
}

//==============================================================//
//================== Scales =====================================//
//==============================================================//

Scales::Scales(Wavelet* wav, SCALETYPE st, int afs, float af0, float af1, int afn) {
    fs = afs;
    scales = new float[afn];
    fourwavl = wav->four_wavelen;
    nscales = afn;

    if (st == SCALETYPE::FCWT_LOGSCALES)
        calculate_logscale_array(2.0f, wav->four_wavelen, afs, af0, af1, afn);
    else if (st == SCALETYPE::FCWT_LINSCALES)
        calculate_linscale_array(wav->four_wavelen, afs, af0, af1, afn);
    else
        calculate_linfreq_array(wav->four_wavelen, afs, af0, af1, afn);
}

void Scales::getScales(float* pfreqs, int pnf) {
    for (int i = 0; i < pnf; i++) {
        pfreqs[i] = scales[i];
    }
}

void Scales::getFrequencies(float* pfreqs, int pnf) {
    for (int i = 0; i < pnf; i++) {
        pfreqs[i] = static_cast<float>(fs) / scales[i];
    }
}

void Scales::calculate_logscale_array(float base, float four_wavl, int fs, float f0, float f1, int fn) {
    // If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
    float nf0 = f0;
    float nf1 = f1;
    float s0 = (fs / nf1);
    float s1 = (fs / nf0);

    // Cannot pass the Nyquist frequency
    assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs / 2));

    float power0 = log(s0) / log(base);
    float power1 = log(s1) / log(base);
    float dpower = power1 - power0;

    for (int i = 0; i < fn; i++) {
        float power = power0 + (dpower / (fn - 1)) * i;
        scales[i] = pow(base, power);
    }
}

void Scales::calculate_linfreq_array(float four_wavl, int fs, float f0, float f1, int fn) {
    float nf0 = f0;
    float nf1 = f1;

    // Cannot pass the Nyquist frequency
    assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs / 2));

    float df = nf1 - nf0;
    for (int i = 0; i < fn; i++) {
        scales[fn - i - 1] = static_cast<float>(fs) / (nf0 + (df / fn) * i);
    }
}

void Scales::calculate_linscale_array(float four_wavl, int fs, float f0, float f1, int fn) {
    float nf0 = f0;
    float nf1 = f1;

    // Cannot pass the Nyquist frequency
    assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs / 2));

    float s0 = fs / nf1;
    float s1 = fs / nf0;
    float ds = s1 - s0;

    for (int i = 0; i < fn; i++) {
        scales[i] = s0 + (ds / fn) * i;
    }
}
