#include <iostream>
#include <fstream>
#include <sstream>
#include <TCanvas.h>
#include <TGraph.h>

void plotWaveformFromTextFile(const char* filename, int lineNumber) {
    // Open the text file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    // Find the specified line number
    std::string line;
    for (int i = 0; i < lineNumber; ++i) {
        if (!std::getline(inputFile, line)) {
            std::cerr << "Error: Specified line number is out of range." << std::endl;
            inputFile.close();
            return;
        }
    }

    // Close the file
    inputFile.close();

    // Read ADC values from the line
    std::istringstream iss(line);
    std::vector<float> adcValues;
    float adcValue;
    while (iss >> adcValue) {
        adcValues.push_back(adcValue);
    }

    // Check if no ADC values were read
    if (adcValues.empty()) {
        std::cerr << "Error: No ADC values found in line " << lineNumber << std::endl;
        return;
    }

    // Create a TGraph from the ADC values
    int numPoints = adcValues.size();
    std::cout << numPoints << std::endl;
    double* x = new double[numPoints];
    double* y = new double[numPoints];
    for (int i = 0; i < numPoints; ++i) {
        x[i] = i;  // Assuming x-axis is just the sample index
        y[i] = adcValues[i];
        cout << adcValues[i] << endl;
    }

    TCanvas* canvas = new TCanvas("canvas", "Waveform Plot", 800, 600);
    TGraph* graph = new TGraph(numPoints, x, y);
    graph->SetTitle("Waveform Plot");
    graph->GetXaxis()->SetTitle("Sample Index");
    graph->GetYaxis()->SetTitle("ADC Value");
    graph->Draw("AL"); // Draw the graph with lines connecting points
    canvas->Draw();
}

int plot_root() {
    int number; 
    for (int i = 0; i < 100; ++i){
        plotWaveformFromTextFile("waveforms.txt", i); // Plot the first waveform
        std::cout << "Enter to continue" << std::endl;
        std::cin >> number; 
   } 

    return 0;
}
