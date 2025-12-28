#ifndef ACQUISITION_H
#define ACQUISITION_H

#include <vector>
#include <complex>
#include <fftw3.h>

// These Settings come from the information of the recorded signal: https://gfix.dk/matlab-gnss-sdr-book/gnss-signal-records/
struct Settings {
    double samplingFreq = 16367600.0; //i'll need to convert this later if I run my own recorded GPS signal bc the rtl-sdr is dif. and also i thinnk it uses uint...
    double ifFreq = 4130400.0;
    double codeFreqBasis = 1023000.0;
    int codeLength = 1023;
    double acqSearchBand = 14.0; // This is >/< the freqeuncy to account for doppler shift from the satellite's motion
    double acqThreshold = 2.5;   // SNR threshold
};


struct AcqResult {
    int prn;
    double carrFreq;    // Hz
    double codePhase;      // alignment time in ms (time of arrival)
    double peakMetric;     // SNR value 
    bool isAcquired;
};

class Acquisition {
public:
    Acquisition(Settings settings);
    ~Acquisition();

    std::vector<AcqResult> process(const std::vector<int8_t>& rawData);
    //std::vector<int> generateCAcode(int prn);

private:
    Settings settings;
    int samplesPerCode;
    
    fftw_complex *fft_in;
    fftw_complex *fft_out;
    fftw_plan fft_plan;
    fftw_plan ifft_plan;

    std::vector<int> generateCAcode(int prn);
    std::vector<std::complex<double>> generateResampledCode(int prn);
};

#endif // ACQUISITION_H