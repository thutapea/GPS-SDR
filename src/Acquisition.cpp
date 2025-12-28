#include "Acquisition.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <complex>

// Helper macros for FFTW indexing
#define RE(x) ((x)[0])
#define IM(x) ((x)[1])

const double PI = 3.14159265358979323846;

Acquisition::Acquisition(Settings s) : settings(s) {
    samplesPerCode = std::round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength));
    fft_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samplesPerCode);
    fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * samplesPerCode);

    // Using ESTIMATE is often safer/faster for one-off planning than MEASURE
    fft_plan = fftw_plan_dft_1d(samplesPerCode, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    ifft_plan = fftw_plan_dft_1d(samplesPerCode, fft_in, fft_out, FFTW_BACKWARD, FFTW_ESTIMATE);
}

Acquisition::~Acquisition() {
    fftw_destroy_plan(fft_plan);
    fftw_destroy_plan(ifft_plan);
    fftw_free(fft_in);
    fftw_free(fft_out);
}

std::vector<int> Acquisition::generateCAcode(int prn) {
    static const int g2_taps[32][2] = {
        {2,6}, {3,7}, {4,8}, {5,9}, {1,9}, {2,10}, {1,8}, {2,9}, {3,10}, {2,3},
        {3,4}, {5,6}, {6,7}, {7,8}, {8,9}, {9,10}, {1,4}, {2,5}, {3,6}, {4,7},
        {5,8}, {6,9}, {1,3}, {4,6}, {5,7}, {6,8}, {7,9}, {8,10}, {1,6}, {2,7},
        {3,8}, {4,9}
    };

    std::vector<int> code(1023);
    int G1[10] = {1,1,1,1,1,1,1,1,1,1};
    int G2[10] = {1,1,1,1,1,1,1,1,1,1};

    for (int i = 0; i < 1023; i++) {
        int tap1 = g2_taps[prn - 1][0] - 1;
        int tap2 = g2_taps[prn-1][1] - 1;
        code[i] = G1[9] ^ G2[tap1] ^ G2[tap2];

        int feedG1 = G1[9] ^ G1[2];
        for (int j = 9; j>0; j--) G1[j] = G1[j-1];
        G1[0] = feedG1;

        int feedG2 = G2[9] ^ G2[8] ^ G2[7] ^ G2[5] ^ G2[2] ^ G2[1];
        for (int j = 9; j > 0; j--) G2[j] = G2[j-1];
        G2[0] = feedG2;
    }
    return code;
}

std::vector<std::complex<double>> Acquisition::generateResampledCode(int prn) {
    std::vector<int> caCode = generateCAcode(prn);
    std::vector<std::complex<double>> resampled(samplesPerCode);
    double ts = 1.0 / settings.samplingFreq;
    double tc = 1.0 / settings.codeFreqBasis;
    
    for (int i = 0; i < samplesPerCode; i++) {
        int codeIdx = (int)floor((i * ts) / tc) % 1023;
        double val = (caCode[codeIdx] == 0) ? -1.0 : 1.0;
        resampled[i] = std::complex<double>(val, 0.0);
    }
    return resampled;
}

std::vector<AcqResult> Acquisition::process(const std::vector<int8_t>& rawData) {
    std::vector<AcqResult> results;
    std::cout << "Starting Acquisition (Samples: " << samplesPerCode << ")..." << std::endl;

    // --- Pre-calculation ---
    std::vector<double> phasePoints(samplesPerCode);
    double ts = 1.0 / settings.samplingFreq;
    for (int i = 0; i < samplesPerCode; i++){
        phasePoints[i] = i * 2.0 * PI * ts;
    } 

    std::vector<double> signal(samplesPerCode);
    for (int i = 0; i < samplesPerCode; i++){
        signal[i] = static_cast<double>(rawData[i]);
    }

    int numFrqBins = (int)(settings.acqSearchBand * 2) + 1;

    // --- Outer Loop: Satellites (PRN) ---
    for (int prn = 1; prn <= 32; prn++) {
        
        // 1. Convert Code to Frequency Domain
        std::vector<std::complex<double>> localCode = generateResampledCode(prn);
        for (int k = 0; k < samplesPerCode; k++){
            RE(fft_in[k]) = localCode[k].real();
            IM(fft_in[k]) = localCode[k].imag();
        }
        
        fftw_execute(fft_plan); // Code -> Freq Domain (result in fft_out)

        // 2. Store Conjugate of Code
        std::vector<std::complex<double>> caCodeFreqDomConj(samplesPerCode);
        for (int k = 0; k < samplesPerCode; k++){
            // FIX: Read from fft_out, not fft_in
            caCodeFreqDomConj[k] = std::complex<double>(RE(fft_out[k]), -IM(fft_out[k]));
        }

        double maxPeak = 0.0;
        double secondPeak = 0.0;
        double bestDoppler = 0.0;
        double bestCodePhase = 0.0;

        // --- Inner Loop: Frequency Bins ---
        for (int i = 0; i < numFrqBins; i++) {
            // FIX: Correct variable name freqStep
            double freqStep = 500.0;
            double carrierFreq = settings.ifFreq - (settings.acqSearchBand * 1000.0/2.0) + (i * freqStep);

            // 3. Carrier Wipeoff
            for (int k = 0; k < samplesPerCode; k++){
                double angle = carrierFreq * phasePoints[k];
                RE(fft_in[k]) = signal[k] * sin(angle);
                IM(fft_in[k]) = signal[k] * cos(angle);
            }

            // 4. Signal to Frequency Domain
            fftw_execute(fft_plan); // Result in fft_out

            // 5. Multiply (Correlation in Freq Domain)
            for (int k = 0; k < samplesPerCode; k++) {
                std::complex<double> sig(RE(fft_out[k]), IM(fft_out[k]));
                std::complex<double> prod = sig * caCodeFreqDomConj[k];
                
                // FIX: Write result to fft_in for the IFFT!
                RE(fft_in[k]) = prod.real();
                IM(fft_in[k]) = prod.imag();
            }

            // 6. Return to Time Domain
            fftw_execute(ifft_plan); // Result in fft_out

            // 7. Find Peak in this frequency bin
            std::vector<double> magSq(samplesPerCode);
            double currentMax = 0.0;
            int currentMaxIdx = 0; // FIX: Added semicolon
            
            for (int k = 0; k < samplesPerCode; k++) {
                magSq[k] = (RE(fft_out[k]) * RE(fft_out[k])) + (IM(fft_out[k]) * IM(fft_out[k]));
                if (magSq[k] > currentMax) {
                    currentMax = magSq[k];
                    currentMaxIdx = k;
                }
            }

            // 8. Update global stats for this PRN
            if (currentMax > maxPeak) {
                maxPeak = currentMax;
                bestDoppler = carrierFreq;
                bestCodePhase = currentMaxIdx;
                
                // Exclusion zone logic
                int samplesPerChip = (int)std::round(settings.samplingFreq / settings.codeFreqBasis);
                int excludeStart = currentMaxIdx - samplesPerChip;
                int excludeEnd = currentMaxIdx + samplesPerChip;
                
                double currentSecondPeak = 0.0;
                for (int k = 0; k < samplesPerCode; k++) {
                    bool inExclusion = false;
                    if (excludeStart < 0) {
                        if (k >= (samplesPerCode + excludeStart) || k <= excludeEnd) inExclusion = true;
                    } else if (excludeEnd >= samplesPerCode) {
                        if (k >= excludeStart || k <= (excludeEnd - samplesPerCode)) inExclusion = true;
                    } else {
                        if (k >= excludeStart && k <= excludeEnd) inExclusion = true;
                    }

                    if (!inExclusion && magSq[k] > currentSecondPeak) {
                        currentSecondPeak = magSq[k];
                    }
                }
                secondPeak = currentSecondPeak;
            }
        } // End Frequency Loop

        // FIX: Ratio check happens HERE (Outside Freq loop, Inside PRN loop)
        if (secondPeak > 0) {
            double ratio = maxPeak / secondPeak;
            
            if (ratio > settings.acqThreshold) {
                std::cout << "Detected PRN " << prn 
                          << " | Doppler: " << bestDoppler 
                          << " | Phase: " << bestCodePhase 
                          << " | Ratio: " << ratio << std::endl;

                AcqResult result;
                result.prn = prn;
                result.carrFreq = bestDoppler;
                result.codePhase = bestCodePhase;
                result.peakMetric = ratio;
                results.push_back(result);
            }
        }

    } // End PRN Loop

    std::cout << "Acquisition complete." << std::endl;
    return results;
}