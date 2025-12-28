#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Acquisition.h"

// Helper function to load binary data
std::vector<int8_t> loadData(const std::string& filename, size_t numSamples) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return {};
    }

    std::vector<int8_t> data(numSamples);
    // Cast the pointer to char* because ifstream reads bytes
    file.read(reinterpret_cast<char*>(data.data()), numSamples);
    
    size_t bytesRead = file.gcount();
    if (bytesRead < numSamples) {
        std::cerr << "Warning: Requested " << numSamples << " samples, but only " << bytesRead << "." << std::endl;
        data.resize(bytesRead);
    }
    
    return data;
}

int main() {
    //  Configured for "gioveAandB_short.bin"
    Settings settings;
    settings.samplingFreq = 16367600.0;    // 16.3676 MHz
    settings.ifFreq = 4130400.0;           // 4.1304 MHz
    settings.codeFreqBasis = 1023000.0;    //  GPS L1 C/A 
    settings.codeLength = 1023;            
    settings.acqSearchBand = 14.0;         // Search +/- 7 kHz (14 kHz total width)
    settings.acqThreshold = 2.5;           // Detection threshold

    // Instantiate
    Acquisition acq(settings);

    //Load Data
    size_t samplesToLoad = (size_t)(settings.samplingFreq * 0.040); 
    std::string filePath = "../data/gioveAandB_short.bin";

    std::cout << "Loading " ;
    std::vector<int8_t> rawSignal = loadData(filePath, samplesToLoad);



    // Acquisition
    std::cout << "Processing" << std::endl;
    std::vector<AcqResult> results = acq.process(rawSignal);

    //  Print Results
    std::cout << "ACQUISITION RESULTS" << std::endl;
    
    if (results.empty()) {
        std::cout << "No satellites detected." << std::endl;
    } else {
        std::cout << "Satellites Found: " << results.size() << std::endl;
        std::cout << "PRN\tDoppler(Hz)\tCodePhase\tMetric" << std::endl;
        for (const auto& res : results) {
            std::cout << res.prn 
                      << "\t" << res.carrFreq 
                      << "\t\t" << res.codePhase 
                      << "\t\t" << res.peakMetric << std::endl;
        }
    }

    return 0;
}