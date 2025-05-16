//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TSystem.h>
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <stdio.h>
#include <string.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <bits/stdc++.h>
#include "Strlen.h"
#include "TAttMarker.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TVirtualX.h"
#include "TVirtualPadEditor.h"
#include "TColor.h"
#include "Rtypes.h"

using std::cout; 
using std::endl;
using namespace std;

/** @brief Properties recorded for each detected pulse in a waveform */
struct pulse_temp {
    double start;  /* Start time of pulse (10% peak) in waveform (ns) */
    double end;    /* End time of pulse (reach baseline) in waveform (ns) */
    double peak;   /* Max amplitude of pulse (photo-electrons) */
    double energy; /* Energy (integral) of pulse (photo-electrons) */
};

/** @brief Properties recorded for each detected pulse in a waveform */
struct pulse {
  double start;   /* Universal start time of pulse (10% peak) in waveform (ns) */
  double end;     /* Universal end time of pulse (reach baseline) in waveform (ns) */
  double peak;    /* Max amplitude of pulse (photo-electrons) */
  double energy;  /* Energy (integral) of pulse (photo-electrons) */
  double number;  /* Number of channels in which we see pulse (photo-electrons) */
  bool single;    /* Is pulse (photo-electrons) timing consistent across all channels */
  bool beam;      /* Tracks whether beam is on or off */
  double trigger; /* Tracks whether trigger is external (2) or internal (16) */
  double length;  /* Length of waveform in number of bins */
  double side_vp_energy; /* Energy (integral) of pulse (photo-electrons) in SIDE veto panels */
  double top_vp_energy; /* Energy (integral) of pulse (photo-electrons) in TOP veto panel */
};

/** @brief Constants for pulse and pulse-edge detection */
const int PULSE_THRESHOLD = 30;     /* Pulse detected if read above this value */
const int BS_UNCERTAINTY = 5;       /* Baseline uncertainty */
const int PULSE_PE_THRESHOLD = 40;  /* Large pulse detected if read above this value (given in units of photoelectrons) */
const int EV61_THRESHOLD = 1200;    /* SNS Beam assumed to be on if Event 61 read above this value (given in units of ADC) */

/** @brief Maximum number of waveforms to process from input root file */
const int MAX_NUM_ENTRIES = 1900000;        // run4144: 2700000 ; run4176: 2500000 ; run4193: 1900000

/** @brief These values should change for each PMT */
const char* inputFilePath = "";
const char* inputFileName = "processed";
const char* outputFilePath = "";                            // What does "" mean?
const char* outputFileName = "PMTWaveformAnalysis";         // Why is output file going into ~ directory?
const char* outputStatsName = "PMTAnalysisStats";

std::vector<double> amplitudeToPE{ 26,28,30,29,27,12,30,28,30,30,24,30 };                // 12/3: run4176

template<typename T>
double getAverage(std::vector<T> const& v) {

    if (v.empty()) {

        return 0;

    }

    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

template<typename T>
int mostFrequent(std::vector<T> const& v) {

    if (v.empty()) {

        return 0;

    }

    int maxcount = 0;
    int element_having_max_freq;
    for (int i = 0; i < v.size(); i++) {
        int count = 0;
        for (int j = 0; j < v.size(); j++) {
            if (v[i] == v[j]) {
                count++;
            }
        }

        if (count > maxcount) {
            maxcount = count;
            element_having_max_freq = v[i];
        }
    }

    if (maxcount > 1) {

        return element_having_max_freq;

    }

    else {

        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();

    }
    
}

template<typename T>
int mostFrequentIndex(std::vector<T> const& v) {

    if (v.empty()) {

        return -1;

    }

    int maxcount = 0;
    int index_having_max_freq;
    for (int i = 0; i < v.size(); i++) {
        int count = 0;
        for (int j = 0; j < v.size(); j++) {
            if (v[i] == v[j]) {
                count++;
            }
        }

        if (count > maxcount) {
            maxcount = count;
            index_having_max_freq = i;
        }
    }

    if (maxcount > 1) {

        return index_having_max_freq;

    }

    else {

        return -2;

    }
    
}

template<typename T>
T variance(const std::vector<T>& vec) {

    const size_t sz = vec.size();

    if (vec.empty()) {
        return 0.0;
    } else if (sz == 1) {
        return 0.0;
    }

    // Calculate the mean
    const T mean = std::accumulate(vec.begin(), vec.end(), 0.0) / sz;

    // Now calculate the variance
    auto variance_func = [&mean, &sz](T accumulator, const T& val) {
        return accumulator + ((val - mean) * (val - mean) / (sz - 1));      // No - 1?
    };

    return std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

template<typename T>
double rmsValue(std::vector<T>& vec)
{
    double square = 0.0;
    double mean = 0.0, root = 0.0;
 
    // Calculate square.
    for (int i = 0; i < vec.size(); i++) {
        square += pow(vec[i], 2);
    }
 
    // Calculate Mean.
    mean = (square / vec.size());
 
    // Calculate Root.
    root = sqrt(mean);
 
    return root;
}

int main(int argc, char *argv[])
{
    int run {0};
    int last_run {0};
    if (argc == 3) {
        sscanf(argv[1], "%i", &run);
        sscanf(argv[2], "%i", &last_run);
    } else if(argc == 2){ 
        sscanf(argv[1], "%i", &run);
        last_run = run;
    } else {
        cout <<"Usase: "<< argv[0] <<" [run]" << " [last run] (optional)" << endl;
        return -1; 
    }   
    cout<<"run: "<< run <<" last run: " << last_run <<endl;

    int ADCSIZE = 45;
    TH1D *h_wf = new TH1D("h_wf", "Waveform", ADCSIZE, 0, ADCSIZE);

    /* Create vectors for plotting histograms */

    int ev_61_dt_max = pow(10, 7);

    int muon_dt_max = pow(10, 5);

    bool before_ev61_peak = true;

    bool bt_HL_n_LL = false;

    int num_events_low_dt = 0;

    int num_events_low_dt_cent_spike = 0;

    int num_events_low_dt_no_spike_neg = 0;

    int num_events_low_dt_no_spike_pos = 0;

    int num_muons_found = 0;

    /* Initialize histograms */                 // NUMBER OF BINS SHOULD BE A MULTIPLE OF THE SMALLEST AXIS UNIT!!! (OR LENGTH OF AXIS DIVIDED BY SMALLEST AXIS UNIT?)

    TH1D* h_dt_61_16 = new TH1D("h_dt_61_16", "Distribution of Time Separations Between Event 61 and Veto Panel 1 Peaks", 3000, - 5 * pow(10, 7), 5 * pow(10, 7));

    TH1D* h_dt_61_17 = new TH1D("h_dt_61_17", "Distribution of Time Separations Between Event 61 and Veto Panel 2 Peaks", 3000, - 5 * pow(10, 7), 5 * pow(10, 7));

    TH1D* h_dt_61_18 = new TH1D("h_dt_61_18", "Distribution of Time Separations Between Event 61 and Veto Panel 3 Peaks", 3000, - 5 * pow(10, 7), 5 * pow(10, 7));

    TH1D* h_dt_61_19 = new TH1D("h_dt_61_19", "Distribution of Time Separations Between Event 61 and Veto Panel 4 Peaks", 3000, - 5 * pow(10, 7), 5 * pow(10, 7));

    TH1D* h_dt_61_20 = new TH1D("h_dt_61_20", "Distribution of Time Separations Between Event 61 and Veto Panel 5 Peaks", 3000, - 5 * pow(10, 7), 5 * pow(10, 7));

    TH1D* h_dt_61_21 = new TH1D("h_dt_61_21", "Distribution of Time Separations Between Event 61 and Veto Panel 6 Peaks", 3000, - 5 * pow(10, 7), 5 * pow(10, 7));

    TH1D* h_dt_61_22 = new TH1D("h_dt_61_22", "Distribution of Time Separations Between Event 61 and Veto Panel 7 Peaks", 3000, - 5 * pow(10, 7), 5 * pow(10, 7));

    TH1D* h_dt_61_23 = new TH1D("h_dt_61_23", "Distribution of Time Separations Between Event 61 and Veto Panel 8 Peaks", 3000, - 5 * pow(10, 7), 5 * pow(10, 7));

    TH1D* h_dt_61_top = new TH1D("h_dt_61_top", "Distribution of Time Separations Between Event 61 and Top Veto Panel Peaks", 3000, - 5 * pow(10, 7), 5 * pow(10, 7));

    TH2D* h_dt_61_v_int = new TH2D("h_dt_61_v_int", "Veto Panel Delta-T vs Veto Panel Integral Values", 3000, - 5 * pow(10, 7), 5 * pow(10, 7), 51000, -1000, 50000);

    TH1D* h_dt_61_16_simple = new TH1D("h_dt_61_16_simple", "Distribution of Time Separations Between Event 61 and Veto Panel 1 Peaks", 170000, 0, 17 * pow(10, 6));

    TH1D* h_dt_61_17_simple = new TH1D("h_dt_61_17_simple", "Distribution of Time Separations Between Event 61 and Veto Panel 2 Peaks", 170000, 0, 17 * pow(10, 6));

    TH1D* h_dt_61_18_simple = new TH1D("h_dt_61_18_simple", "Distribution of Time Separations Between Event 61 and Veto Panel 3 Peaks", 170000, 0, 17 * pow(10, 6));

    TH1D* h_dt_61_19_simple = new TH1D("h_dt_61_19_simple", "Distribution of Time Separations Between Event 61 and Veto Panel 4 Peaks", 170000, 0, 17 * pow(10, 6));

    TH1D* h_dt_61_20_simple = new TH1D("h_dt_61_20_simple", "Distribution of Time Separations Between Event 61 and Veto Panel 5 Peaks", 170000, 0, 17 * pow(10, 6));

    TH1D* h_dt_61_21_simple = new TH1D("h_dt_61_21_simple", "Distribution of Time Separations Between Event 61 and Veto Panel 6 Peaks", 170000, 0, 17 * pow(10, 6));

    TH1D* h_dt_61_22_simple = new TH1D("h_dt_61_22_simple", "Distribution of Time Separations Between Event 61 and Veto Panel 7 Peaks", 170000, 0, 17 * pow(10, 6));

    TH1D* h_dt_61_23_simple = new TH1D("h_dt_61_23_simple", "Distribution of Time Separations Between Event 61 and Veto Panel 8 Peaks", 170000, 0, 17 * pow(10, 6));

    TH1D* h_dt_61_top_simple = new TH1D("h_dt_61_top_simple", "Distribution of Time Separations Between Event 61 and Top Veto Panel Peaks", 170000, 0, 17 * pow(10, 6));

    TH2D* h_dt_61_v_int_simple = new TH2D("h_dt_61_v_int_simple", "Veto Panel Delta-T vs Veto Panel Integral Values", 17000, 0,  17 * pow(10, 6), 51000, -1000, 50000);

    // TH1D* h_int = new TH1D("h_int", "Distribution of All Event Integral Values Between High Light and Low Light LED Events", 5500, -500, 5000);

    // TH1D* h_int_cent_spike = new TH1D("h_int_cent_spike", "Distribution of All Central Spike Detector Integral Values", 5500, -500, 5000);

    // TH1D* h_int_pos_low_dt = new TH1D("h_int_pos_low_dt", "Distribution of Low dt Detector Integral Values", 5500, -500, 5000);

    // TH1D* h_int_neg_low_dt = new TH1D("h_int_neg_low_dt", "Distribution of Low dt Detector Integral Values", 5500, -500, 5000);

    // TH1D* h_svpint_cent_spike = new TH1D("h_svpint_cent_spike", "Distribution of All Central Spike Side Veto Panel Integral Values", 240, -400, 2000);

    // TH1D* h_svpint_pos_low_dt = new TH1D("h_svpint_pos_low_dt", "Distribution of Low dt Side Veto Panel Integral Values", 240, -400, 2000);

    // TH1D* h_svpint_neg_low_dt = new TH1D("h_svpint_neg_low_dt", "Distribution of Low dt Side Veto Panel Integral Values", 240, -400, 2000);

    // TH1D* h_tvpint_cent_spike = new TH1D("h_tvpint_cent_spike", "Distribution of All Central Spike Top Veto Panel Integral Values", 60, -200, 400);

    // TH1D* h_tvpint_pos_low_dt = new TH1D("h_tvpint_pos_low_dt", "Distribution of Low dt Top Veto Panel Integral Values", 60, -200, 400);

    // TH1D* h_tvpint_neg_low_dt = new TH1D("h_tvpint_neg_low_dt", "Distribution of Low dt Top Veto Panel Integral Values", 60, -200, 400);

    // TH1D* h_int_HL = new TH1D("h_int_HL", "Distribution of All High Light LED Integral Values", 5500, -500, 5000);

    // TH2D* h_dt_v_tmuon = new TH2D("h_dt_v_tmuon", "Delta-T vs Most Recent Muon Time Difference", 1000, - 2 * pow(10, 7), 2 * pow(10, 7), 1000, 0, pow(10, 8));

    for(int iRun = run; iRun <= last_run ; iRun++){

        TFile *f;
        //your root file location here
///        if(gSystem->AccessPathName(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v4.root", iRun))){
        if (gSystem->AccessPathName(Form("/data9/coherent/data/d2o/processedData/run%i_processed_v4.root", iRun))) {
///        if(gSystem->AccessPathName(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v4.root", iRun))){
            cout << "Could not open file" << endl;
            return -1; 
        } else{
///            f = new TFile(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/run%i_processed_v4.root", iRun));
            f = new TFile(Form("/data9/coherent/data/d2o/processedData/run%i_processed_v4.root", iRun));
///            f = new TFile(Form("/mnt/c/Users/eliwa/OneDrive/Documents/Detector_Data_Analysis/run%i_processed_v4.root", iRun));
        }

        // Create a text string, which is used to output the text file

        std::vector<double> integralToPE;

        integralToPE.clear();

        string int_val;

        // Read from the text file

        ifstream ReadSPIFile(Form("Single_Phe_Integral_Values_Run%i.txt", iRun));

        // Use a while loop together with the getline() function to read the file line by line

        while (getline(ReadSPIFile, int_val)) {

            // Output the text from the file

            integralToPE.push_back(stod(int_val));
        }

        // Close the file

        ReadSPIFile.close();

///        TTree *t = (TTree *) f->Get("wf");
        TTree* t = (TTree*)f->Get("tree");

        // Declaration of leaf types
        Int_t           eventID;
        Int_t           nSamples[32];
///        UInt_t          adcTime;
///        Int_t           adcSize[32];
        Short_t         adcVal[32][45];
///        Int_t           trigPattern;
///        ULong64_t       timeStamp_extTrig;
        Double_t        baselineMean[32];
        Double_t        baselineRMS[32];
///        Int_t           peakBin[32];
        Double_t        pulseH[32];
        Int_t           peakPosition[32];
        Double_t        area[32];
        Long64_t        nsTime;
        Int_t           triggerBits;

        // List of branches
        TBranch        *b_eventID;
        TBranch        *b_nSamples;
///        TBranch        *b_adcTime;
///        TBranch        *b_adcSize;
        TBranch        *b_adcVal;
///        TBranch        *b_trigPattern;
///        TBranch        *b_timeStamp_extTrig;
        TBranch        *b_baselineMean;
        TBranch        *b_baselineRMS;
///        TBranch        *b_peakBin;
        TBranch        *b_pulseH;
        TBranch        *b_peakPosition;
        TBranch        *b_area;
        TBranch        *b_nsTime;
        TBranch        *b_triggerBits;

        t->SetBranchAddress("eventID", &eventID, &b_eventID);
        t->SetBranchAddress("nSamples", &nSamples, &b_nSamples);
///        t->SetBranchAddress("adcTime", &adcTime, &b_adcTime);
///        t->SetBranchAddress("adcSize", adcSize, &b_adcSize);
        t->SetBranchAddress("adcVal", adcVal, &b_adcVal);
///        t->SetBranchAddress("trigPattern", &trigPattern, &b_trigPattern);
///        t->SetBranchAddress("timeStamp_extTrig", &timeStamp_extTrig, &b_timeStamp_extTrig);
        t->SetBranchAddress("baselineMean", baselineMean, &b_baselineMean);
        t->SetBranchAddress("baselineRMS", baselineRMS, &b_baselineRMS);
///        t->SetBranchAddress("peakBin", peakBin, &b_peakBin);
        t->SetBranchAddress("pulseH", pulseH, &b_pulseH);
        t->SetBranchAddress("peakPosition", &peakPosition, &b_peakPosition);
        t->SetBranchAddress("area", area, &b_area);
        t->SetBranchAddress("nsTime", &nsTime, &b_nsTime);
        t->SetBranchAddress("triggerBits", &triggerBits, &b_triggerBits);

	    // Create output root file of processed waveform data
	    TString outputFile;
	    if (strcmp(outputFilePath, "") == 0) {
	      outputFile.Form("%s.root", outputFileName);
	    } else {
	      outputFile.Form("%s/%s.root", outputFilePath, outputFileName);
	    }
	    TFile *fileOut = new TFile(outputFile, "RECREATE");
	    TTree *michelTree = new TTree("michelTree", "michelTree");
	    michelTree->SetDirectory(fileOut);

	    // Data to keep track of for each michelTree entry
	    int *br_entry;  /* Entry # from inputFile TTree T */
	    double *br_e1;  /* Energy of first pulse (photo-electrons) */
	    double *br_p1;  /* Peak (max amplitude) of first pulse (photo-electrons) */
	    double *br_t1;  /* Universal start time (10% peak) of first pulse (nanoseconds) */
	    double *br_d1;  /* Duration of first pulse (nanoseconds) */
	    double *br_e2;  /* Energy of second pulse (photo-electrons) */
	    double *br_p2;  /* Peak (max amplitude) of second pulse (photo-electrons) */
	    double *br_t2;  /* Universal start time (10% peak) of second pulse (nanoseconds) */
	    double *br_d2;  /* Duration of second pulse (nanoseconds) */
	    double *br_dt;  /* Time separation between pulse onsets (nanoseconds) */
	    bool *br_issue; /* Flag to keep track of unusual michelTree entries */
        double* br_n1;  /* Number of channels in which we see first pulse (photo-electrons) */
        bool* br_s1;    /* Is first pulse (photo-electrons) timing consistent across all channels */
        bool* br_b1;    /* Tracks whether beam is on or off for first pulse */
        double* br_tr1; /* Tracks whether trigger is external (2) or internal (16) for first pulse */
        double* br_l1;  /* Length of first waveform in number of bins */
        double* br_n2;  /* Number of channels in which we see second pulse (photo-electrons) */
        bool* br_s2;    /* Is second pulse (photo-electrons) timing consistent across all channels */
        bool* br_b2;    /* Tracks whether beam is on or off for second pulse */
        double* br_tr2; /* Tracks whether trigger is external (2) or internal (16) for second pulse */
        double* br_l2;  /* Length of second waveform in number of bins */

	    michelTree->Branch("entry", &br_entry, "entry/I");
	    michelTree->Branch("e1", &br_e1, "e1/d");
	    michelTree->Branch("p1", &br_p1, "p1/d");
	    michelTree->Branch("t1", &br_t1, "t1/d");
	    michelTree->Branch("d1", &br_d1, "d1/d");
	    michelTree->Branch("e2", &br_e2, "e2/d");
	    michelTree->Branch("p2", &br_p2, "p2/d");
	    michelTree->Branch("t2", &br_t2, "t2/d");
	    michelTree->Branch("d2", &br_d2, "d2/d");
	    michelTree->Branch("dt", &br_dt, "dt/d");
	    michelTree->Branch("issue", &br_issue, "issue/O");
        michelTree->Branch("n1", &br_n1, "n1/d");
        michelTree->Branch("s1", &br_s1, "s1/O");
        michelTree->Branch("b1", &br_b1, "b1/O");
        michelTree->Branch("tr1", &br_tr1, "tr1/d");
        michelTree->Branch("l1", &br_l1, "l1/d");
        michelTree->Branch("n2", &br_n2, "n2/d");
        michelTree->Branch("s2", &br_s2, "s2/O");
        michelTree->Branch("b2", &br_b2, "b2/O");
        michelTree->Branch("tr2", &br_tr2, "tr2/d");
        michelTree->Branch("l2", &br_l2, "l2/d");

        // Get statistics for up to 1 million entries from fileIn TTree T
        int numEntries = std::min((int)t->GetEntries(), MAX_NUM_ENTRIES);

        std::vector<double> chan_lengths;

        std::vector<double> peak_pos_RMS;

        double avg_peak_pos_RMS = 0.0;

        double var_peak_pos_RMS = 0.0;

        int ch15_on = 0;

        double last_muon_time = 0.0;

        double last_vp_time = 0.0;

        double last_ev61_time = 0.0;

        double last_muon_time_simple = 0.0;

        double last_vp_time_simple = 0.0;

        double last_ev61_time_simple = 0.0;

        double last_HL_time = 0.0;

        double last_LL_time = 0.0;

        double last_vp_time_16 = 0.0;

        double last_vp_time_17 = 0.0;

        double last_vp_time_18 = 0.0;

        double last_vp_time_19 = 0.0;

        double last_vp_time_20 = 0.0;

        double last_vp_time_21 = 0.0;

        double last_vp_time_22 = 0.0;

        double last_vp_time_23 = 0.0;

        double last_vp_time_top = 0.0;

        double last_vp_time_16_simple = 0.0;

        double last_vp_time_17_simple = 0.0;

        double last_vp_time_18_simple = 0.0;

        double last_vp_time_19_simple = 0.0;

        double last_vp_time_20_simple = 0.0;

        double last_vp_time_21_simple = 0.0;

        double last_vp_time_22_simple = 0.0;

        double last_vp_time_23_simple = 0.0;

        double last_vp_time_top_simple = 0.0;

        std::vector<double> all_chan_start;
        std::vector<double> all_chan_end;
        std::vector<double> all_chan_peak;
        std::vector<double> all_chan_energy;
        std::vector<double> peak_energy;
        std::vector<double> ten_chan_peak;
        std::vector<double> all_chan_peakbin;
        std::vector<double> ten_chan_peakbin;
        std::vector<double> side_veto_panel_energy;
        std::vector<double> top_veto_panel_energy;
        std::vector<double> all_chan_start_vp;
        std::vector<double> all_chan_start_vp_iChan;
    
        // Can use either "t->GetEntries()" or "numEntries" or number
        for (int iEnt = 0; iEnt < t->GetEntries(); iEnt++) {

            std::cout << "\n" << "Processing event " << iEnt + 1 << " of " << t->GetEntries() << " in Run " << iRun << "\n";

            Long64_t tentry = t->LoadTree(iEnt);

            b_eventID->GetEntry(tentry);
            b_nSamples->GetEntry(tentry);
            ///            b_adcTime->GetEntry(tentry);
            ///            b_adcSize->GetEntry(tentry);
            b_adcVal->GetEntry(tentry);
            ///            b_trigPattern->GetEntry(tentry);
            ///            b_timeStamp_extTrig->GetEntry(tentry);
            b_baselineMean->GetEntry(tentry);
            b_baselineRMS->GetEntry(tentry);
            ///            b_peakBin->GetEntry(tentry);
            b_pulseH->GetEntry(tentry);
            b_peakPosition->GetEntry(tentry);
            b_area->GetEntry(tentry);
            b_nsTime->GetEntry(tentry);
            b_triggerBits->GetEntry(tentry);

            // Create vector to hold info of all pulses detected
            // std::vector<struct pulse> pulses;
            std::vector<struct pulse_temp> pulses_temp;

            // Digitizer cannot read above this so flag entries w/ p1,p2 > maxPeak
            double maxPeak = 15776 / amplitudeToPE[0];

            // Create variables to hold info of each current pulse
            bool onPulse = false, onLastPulseTail = false;
            int thresholdBin = 0, peakBin = 0;
            double tailWindow = 250;
            double peak = 0.;
            double pulseEnergy = 0.;
            double AllPulseEnergy = 0.;
            double Ev61Energy = 0.;
            double chan_15_start = 0.;

            bool all_chan_beam = false;

            bool top_vp_event = false;

            int num_chan = 0;

            int num_chan_over_1phe = 0;

            int num_chan_simple = 0;

            double num_chan_vp = 0;

            double p_int_16 = 0;
            double p_int_17 = 0;
            double p_int_18 = 0;
            double p_int_19 = 0;
            double p_int_20 = 0;
            double p_int_21 = 0;
            double p_int_22 = 0;
            double p_int_23 = 0;
            double p_int_24 = 0;
            double p_int_25 = 0;

            bool peak_in_chan_16 = false;
            bool peak_in_chan_17 = false;
            bool peak_in_chan_18 = false;
            bool peak_in_chan_19 = false;
            bool peak_in_chan_20 = false;
            bool peak_in_chan_21 = false;
            bool peak_in_chan_22 = false;
            bool peak_in_chan_23 = false;
            bool peak_in_chan_24 = false;
            bool peak_in_chan_25 = false;

            chan_lengths.clear();

            for (int iChan = 0; iChan < 26; iChan++) {

                if (iChan == 12 || iChan == 13 || iChan == 14) {

                    continue;

                }

                for (int i = 0; i < ADCSIZE; i++) {

                    h_wf->SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]);

                }

                // h_wf->GetXaxis()->SetTitle("Time (1 bin = 16 ns)");
                // h_wf->GetYaxis()->SetTitle("Counts");
                // h_wf->Draw();
                // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", iRun, iEnt + 1, iChan));
                // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/wf_image/wf_run_%i_event_%i_ch_%i.png", iRun, iEnt + 1, iChan));

                if (iChan == 15) {

                    for (int iBin = 1; iBin <= h_wf->GetNbinsX(); iBin++) {

                        double iBinContent = h_wf->GetBinContent(iBin);

                        Ev61Energy += iBinContent;

                    }

                    if (Ev61Energy > EV61_THRESHOLD) {

                        all_chan_beam = true;

                        ch15_on += 1;

                    }

                    Ev61Energy = 0.;

                }

                pulses_temp.clear();

                onLastPulseTail = false;

                bool peak_in_chan = false;

                bool peak_in_chan_simple = false;

                for (int iBin = 1; iBin <= h_wf->GetNbinsX(); iBin++) {

                    double iBinContent = h_wf->GetBinContent(iBin);
                    // std::cout << iBinContent << "\t";

                    if (iChan <= 11 && !peak_in_chan_simple && iBinContent >= PULSE_THRESHOLD) {

                        peak_in_chan_simple = true;

                        num_chan_simple += 1;

                    }
                    
                    if (iChan == 16 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_16 = true;}

                    if (iChan == 17 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_17 = true;}

                    if (iChan == 18 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_18 = true;}

                    if (iChan == 19 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_19 = true;}

                    if (iChan == 20 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_20 = true;}

                    if (iChan == 21 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_21 = true;}

                    if (iChan == 22 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_22 = true;}

                    if (iChan == 23 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_23 = true;}

                    if (iChan == 24 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_24 = true;}

                    if (iChan == 25 && iBin <= ADCSIZE - 3 && iBinContent >= PULSE_THRESHOLD + 1) {peak_in_chan_25 = true;}
                    
                    if (iBin > 15) {

                        AllPulseEnergy += iBinContent;

                    }

                    if (pulses_temp.size() > 0) {    // Could get rid of this check b/c event window is so small

                        onLastPulseTail = iBin * 16 - pulses_temp.back().start < tailWindow;

                    }

                    // Find pulse
                    if (!onPulse && !onLastPulseTail && iBinContent >= PULSE_THRESHOLD) {

                        onPulse = true;
                        thresholdBin = iBin;
                        peakBin = iBin;
                        peak = iBinContent;
                        pulseEnergy += iBinContent;

                    }

                    // Pulse found. Find pulse duration & energy
                    else if (onPulse) {

                        // Accumlate energy of pulse after threshold bin
                        pulseEnergy += iBinContent;

                        // Update pulse's peak
                        if (peak < iBinContent) {
                            peak = iBinContent;
                            peakBin = iBin;
                        }

                        // Search for end of pulse (falls below noiselevel)
                        // Assumes no pulse pileup
                        if (iBinContent < BS_UNCERTAINTY) {

                            // std::cout << "\n" << "Found a pulse in Event " << iEnt + 1 << " Channel " << iChan << "\n";

                            // Create pulse info
                            struct pulse_temp p;        // Clear these struct variables?
                            p.start = thresholdBin * 16.;
                            p.peak = peak / amplitudeToPE[iChan];
                            // Record end of pulse
                            p.end = iBin * 16.;
                            // Accumlate energy of pulse before threshold bin
                            for (int j = peakBin - 1; BS_UNCERTAINTY < h_wf->GetBinContent(j); j--)
                            {
                                if (j < thresholdBin) {     // Technically, this if statement is a bit time inefficient
                                    pulseEnergy += h_wf->GetBinContent(j);
                                }
                                // Record start of pulse (10% of peak)
                                if (h_wf->GetBinContent(j) > peak * 0.1) {
                                    p.start = j * 16.;
                                }
                            }

                            if (iChan <= 11) {

                                // Record energy of pulse
                                p.energy = pulseEnergy / integralToPE[iChan];

                                all_chan_start.push_back(p.start);      // This is triggering for multiple bins above threshold? Even within single peak?
                                all_chan_end.push_back(p.end);
                                peak_energy.push_back(p.energy);
                                // all_chan_peak.push_back(p.peak);

                                if (!peak_in_chan) {

                                    num_chan += 1;

                                    all_chan_peak.push_back(p.peak);

                                    ten_chan_peak.push_back(p.peak);

                                    all_chan_peakbin.push_back(peakBin);

                                    ten_chan_peakbin.push_back(peakBin);

                                    // all_chan_energy.push_back(p.energy);

                                    peak_in_chan = true;

                                }

                                // Add pulse info to vector of pulses
                                pulses_temp.push_back(p);

                            }

                            if (iChan == 15) {

                                chan_15_start = p.start;

                            }

                            if (iChan >= 16) {

                                all_chan_start_vp.push_back(p.start);

                                all_chan_start_vp_iChan.push_back(iChan);

                                // std::cout << "\n" << "p.start = " << p.start << "\n";

                            }

                            // Clear current pulse variables to look for new pulse
                            peak = 0.;
                            peakBin = 0;
                            pulseEnergy = 0.;
                            thresholdBin = 0;
                            onPulse = false;
                        }

                        if (iChan <= 11 && iBin == ADCSIZE) {
                            
                            // Clear current pulse variables to look for new pulse
                            peak = 0.;
                            peakBin = 0;
                            pulseEnergy = 0.;
                            thresholdBin = 0;
                            onPulse = false;

                        }

                    }

                }

                if (iChan <= 11) {

                    all_chan_energy.push_back(AllPulseEnergy / integralToPE[iChan]);

                    if (AllPulseEnergy / integralToPE[iChan] > 1) {

                        num_chan_over_1phe += 1;

                    }

                    if (!peak_in_chan) {

                        ten_chan_peak.push_back(0);

                        ten_chan_peakbin.push_back(0);

                    }

                }

                if (iChan >= 16 && iChan <= 23) {

                    if (AllPulseEnergy > 200) {

                        num_chan_vp += 1;

                    }

                    if (iChan == 16) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_16 = AllPulseEnergy;

                    }

                    else if (iChan == 17) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_17 = AllPulseEnergy;

                    }

                    else if (iChan == 18) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_18 = AllPulseEnergy;

                    }

                    else if (iChan == 19) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_19 = AllPulseEnergy;

                    }

                    else if (iChan == 20) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_20 = AllPulseEnergy;

                    }

                    else if (iChan == 21) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_21 = AllPulseEnergy;

                    }

                    else if (iChan == 22) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_22 = AllPulseEnergy;

                    }

                    else if (iChan == 23) {

                        side_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_23 = AllPulseEnergy;

                    }

                }

                if (iChan >= 24) {

                    if (!top_vp_event && AllPulseEnergy > 200) {

                        num_chan_vp += 1;

                        top_vp_event = true;
                        
                    }

                    if (iChan == 24) {

                        top_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_24 = 1.07809 * AllPulseEnergy;

                    }

                    else if (iChan == 25) {

                        top_veto_panel_energy.push_back(AllPulseEnergy);

                        p_int_25 = AllPulseEnergy;

                    }

                }

                AllPulseEnergy = 0.;

                h_wf->Reset();

                chan_lengths.push_back(nSamples[iChan]);

            }   // Channel loop

            double pulse_start_time = mostFrequent(all_chan_start);

            double vp_start_time = nsTime + mostFrequent(all_chan_start_vp);

            if (num_chan_over_1phe >= 10) {

                int pulse_start_index = mostFrequentIndex(all_chan_start);

            }

            //if (pulse_start_time == 0) {

                //continue;

            //}

            struct pulse avg_pulse;

            avg_pulse.start = nsTime + pulse_start_time;                                                       // Need a way to weed out secondary pulses - use weighted average?
            avg_pulse.end = nsTime + mostFrequent(all_chan_end);
            avg_pulse.energy = std::accumulate(all_chan_energy.begin(), all_chan_energy.end(), 0.0);           // Integrals still sum over secondary peaks discarded in .start and .end values 
            avg_pulse.peak = std::accumulate(all_chan_peak.begin(), all_chan_peak.end(), 0.0);
            avg_pulse.number = num_chan_over_1phe;
            avg_pulse.beam = all_chan_beam;
            avg_pulse.trigger = triggerBits;
            avg_pulse.length = getAverage(chan_lengths);
            avg_pulse.side_vp_energy = std::accumulate(side_veto_panel_energy.begin(), side_veto_panel_energy.end(), 0.0);
            avg_pulse.top_vp_energy = std::accumulate(top_veto_panel_energy.begin(), top_veto_panel_energy.end(), 0.0);

            // std::cout << "\n" << "Run Number = " << iRun << ", Event Number = " << iEnt + 1 << ", triggerBits = " << avg_pulse.trigger << "\n";
            /*
            if (avg_pulse.beam) {

                std::cout << "\n" << "Event 61 event" << "\n";

            }

            if (avg_pulse.trigger == 4) {

                std::cout << "\n" << "Min bias LED event (tB == 4)" << "\n";

            }

            else if (avg_pulse.trigger == 8) {

                std::cout << "\n" << "High light LED event (tB == 8)" << "\n";

            }

            else if (avg_pulse.trigger == 16) {

                std::cout << "\n" << "Low light LED event (tB == 16)" << "\n";

            }
            */
            std::vector<double> chan_start_no_outliers;

            for (int iPeak = 0; iPeak < all_chan_start.size(); iPeak++) {

                if (all_chan_start[iPeak] < (pulse_start_time + 10 * 16) && all_chan_start[iPeak] > (pulse_start_time - 10 * 16)) {

                    chan_start_no_outliers.push_back(all_chan_start[iPeak]);

                }

                else {

                    // Check if pulse outside this range occurs more than once

                    for (int jPeak = 0; jPeak < all_chan_start.size(); jPeak++) {

                        if (iPeak == jPeak) {

                            continue;

                        }

                        else if (all_chan_start[iPeak] < (all_chan_start[jPeak] + 1 * 16) && all_chan_start[iPeak] > (all_chan_start[jPeak] - 1 * 16)) {

                            chan_start_no_outliers.push_back(all_chan_start[iPeak]);

                        }

                    }

                }
            }

            double var_val = variance(chan_start_no_outliers);

            // h_var->Fill(var_val);

            if (var_val < 5 * 16) {            // Pulses largely not consistent, need to weed out outlier and secondary pulses

                avg_pulse.single = true;

            }

            else {

                avg_pulse.single = false;

            }

            num_chan = 0;
            num_chan_over_1phe = 0;
            num_chan_simple = 0;

            all_chan_beam = false;

            top_vp_event = false;

            all_chan_start.clear();
            all_chan_end.clear();
            all_chan_peak.clear();
            all_chan_energy.clear();
            peak_energy.clear();
            ten_chan_peak.clear();
            all_chan_peakbin.clear();
            ten_chan_peakbin.clear();

            chan_lengths.clear();

            // Event 61 - Veto Panel dt Plot: Simple Version

            if (avg_pulse.beam) {last_ev61_time_simple = nsTime + chan_15_start;}

            if (last_ev61_time_simple != 0 && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                if (peak_in_chan_16 || peak_in_chan_17 || peak_in_chan_18 || peak_in_chan_19 || peak_in_chan_20 || peak_in_chan_21 || peak_in_chan_22 || peak_in_chan_23 || peak_in_chan_24 || peak_in_chan_25) {

                    last_vp_time_simple = vp_start_time;

                    std::cout << "\n" << "all_chan_start_vp_iChan.size() = " << all_chan_start_vp_iChan.size() << "\n";

                    for (int iVec = 0; iVec < all_chan_start_vp_iChan.size(); iVec++) {

                        if (all_chan_start_vp_iChan[iVec] == 16 && peak_in_chan_16) {last_vp_time_16_simple = nsTime + all_chan_start_vp[iVec]; h_dt_61_16_simple->Fill(last_vp_time_16_simple - last_ev61_time_simple);}

                        else if (all_chan_start_vp_iChan[iVec] == 17 && peak_in_chan_17) {last_vp_time_17_simple = nsTime + all_chan_start_vp[iVec]; h_dt_61_17_simple->Fill(last_vp_time_17_simple - last_ev61_time_simple);}

                        else if (all_chan_start_vp_iChan[iVec] == 18 && peak_in_chan_18) {last_vp_time_18_simple = nsTime + all_chan_start_vp[iVec]; h_dt_61_18_simple->Fill(last_vp_time_18_simple - last_ev61_time_simple);}

                        else if (all_chan_start_vp_iChan[iVec] == 19 && peak_in_chan_19) {last_vp_time_19_simple = nsTime + all_chan_start_vp[iVec]; h_dt_61_19_simple->Fill(last_vp_time_19_simple - last_ev61_time_simple);}

                        else if (all_chan_start_vp_iChan[iVec] == 20 && peak_in_chan_20) {last_vp_time_20_simple = nsTime + all_chan_start_vp[iVec]; h_dt_61_20_simple->Fill(last_vp_time_20_simple - last_ev61_time_simple);}

                        else if (all_chan_start_vp_iChan[iVec] == 21 && peak_in_chan_21) {last_vp_time_21_simple = nsTime + all_chan_start_vp[iVec]; h_dt_61_21_simple->Fill(last_vp_time_21_simple - last_ev61_time_simple);}

                        else if (all_chan_start_vp_iChan[iVec] == 22 && peak_in_chan_22) {last_vp_time_22_simple = nsTime + all_chan_start_vp[iVec]; h_dt_61_22_simple->Fill(last_vp_time_22_simple - last_ev61_time_simple);}

                        else if (all_chan_start_vp_iChan[iVec] == 23 && peak_in_chan_23) {last_vp_time_23_simple = nsTime + all_chan_start_vp[iVec]; h_dt_61_23_simple->Fill(last_vp_time_23_simple - last_ev61_time_simple);}

                        else if ((all_chan_start_vp_iChan[iVec] == 24 && all_chan_start_vp_iChan[iVec] == 25) && (peak_in_chan_24 || peak_in_chan_25)) {last_vp_time_top_simple = nsTime + (all_chan_start_vp[iVec] + all_chan_start_vp[iVec]) / 2; h_dt_61_top_simple->Fill(last_vp_time_top_simple - last_ev61_time_simple);}

                    }

                    h_dt_61_v_int_simple->Fill((last_vp_time_simple - last_ev61_time_simple), (avg_pulse.side_vp_energy + avg_pulse.top_vp_energy));

                }

            }

            /*
            
            New Ev61-Det dt plot code loop
            
            */

            // Only looking for events between high light LED pulse and low light LED pulse

            // if (avg_pulse.trigger == 8) {bt_HL_n_LL = true; last_HL_time = vp_start_time; h_int_HL->Fill(avg_pulse.energy);}

            // else if (avg_pulse.trigger == 16) {bt_HL_n_LL = false; before_ev61_peak = true; last_muon_time = 0.0; last_vp_time = 0.0; last_ev61_time = 0.0; last_LL_time = vp_start_time; last_vp_time_16 = 0.0; last_vp_time_17 = 0.0; last_vp_time_18 = 0.0; last_vp_time_19 = 0.0; last_vp_time_20 = 0.0; last_vp_time_21 = 0.0; last_vp_time_22 = 0.0; last_vp_time_23 = 0.0; last_vp_time_top = 0.0;}

            // Stop looking (only) for detector pulses after Event 61 pulses

            // if (vp_start_time - last_ev61_time > ev_61_dt_max) {before_ev61_peak = true;}

            // Is this event an internally triggered detector pulse?

            // if (bt_HL_n_LL) {         // Should not limit when/where I am searching for internally triggered detector peaks? No "bt_HL_n_LL" check here?
                
            if (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                // if (peak_in_chan_16 || peak_in_chan_17 || peak_in_chan_18 || peak_in_chan_19 || peak_in_chan_20 || peak_in_chan_21 || peak_in_chan_22 || peak_in_chan_23 || peak_in_chan_24 || peak_in_chan_25) {

                // std::cout << "\n" << "peak_in_chan_16-25 = " << peak_in_chan_16 << ", " << peak_in_chan_17 << ", " << peak_in_chan_18 << ", " << peak_in_chan_19 << ", " << peak_in_chan_20 << ", " << peak_in_chan_21 << ", " << peak_in_chan_22 << ", " << peak_in_chan_23 << ", " << peak_in_chan_24 << ", " << peak_in_chan_25 << "\n";

                if (p_int_16 >= 200 || p_int_17 >= 200 || p_int_18 >= 200 || p_int_19 >= 200 || p_int_20 >= 200 || p_int_21 >= 200 || p_int_22 >= 200 || p_int_23 >= 200 || p_int_24 >= 200 || p_int_25 >= 200) {

                    // h_int->Fill(avg_pulse.energy);

                    last_vp_time = vp_start_time;

                    for (int iVec = 0; iVec < all_chan_start_vp_iChan.size(); iVec++) {

                        if (all_chan_start_vp_iChan[iVec] == 16 && p_int_16 >= 200) {last_vp_time_16 = nsTime + all_chan_start_vp[iVec];}

                        else if (all_chan_start_vp_iChan[iVec] == 17 && p_int_17 >= 200) {last_vp_time_17 = nsTime + all_chan_start_vp[iVec];}

                        else if (all_chan_start_vp_iChan[iVec] == 18 && p_int_18 >= 200) {last_vp_time_18 = nsTime + all_chan_start_vp[iVec];}

                        else if (all_chan_start_vp_iChan[iVec] == 19 && p_int_19 >= 200) {last_vp_time_19 = nsTime + all_chan_start_vp[iVec];}

                        else if (all_chan_start_vp_iChan[iVec] == 20 && p_int_20 >= 200) {last_vp_time_20 = nsTime + all_chan_start_vp[iVec];}

                        else if (all_chan_start_vp_iChan[iVec] == 21 && p_int_21 >= 200) {last_vp_time_21 = nsTime + all_chan_start_vp[iVec];}

                        else if (all_chan_start_vp_iChan[iVec] == 22 && p_int_22 >= 200) {last_vp_time_22 = nsTime + all_chan_start_vp[iVec];}

                        else if (all_chan_start_vp_iChan[iVec] == 23 && p_int_23 >= 200) {last_vp_time_23 = nsTime + all_chan_start_vp[iVec];}

                        else if ((all_chan_start_vp_iChan[iVec] == 24 && all_chan_start_vp_iChan[iVec] == 25) && (p_int_24 >= 200 || p_int_25 >= 200)) {last_vp_time_top = nsTime + (all_chan_start_vp[iVec] + all_chan_start_vp[iVec]) / 2;}

                    }
                    
                    // Is this event ALSO an Event 61 event?

                    if (avg_pulse.beam) {

                        last_ev61_time = nsTime + chan_15_start;

                        before_ev61_peak = false;
                        
                    }

                    if (before_ev61_peak) {

                        continue;

                    }

                    else if (!before_ev61_peak) {           // Need another cut for first detector check to see if there was recent Event 61 peak?

                        if (last_ev61_time != 0.0 && last_vp_time_16 != 0.0 && (last_vp_time_16 - last_ev61_time) >= -720 && (last_vp_time_16 - last_ev61_time) <= ev_61_dt_max) {h_dt_61_16->Fill(last_vp_time_16 - last_ev61_time);}

                        if (last_ev61_time != 0.0 && last_vp_time_17 != 0.0 && (last_vp_time_17 - last_ev61_time) >= -720 && (last_vp_time_17 - last_ev61_time) <= ev_61_dt_max) {h_dt_61_17->Fill(last_vp_time_17 - last_ev61_time);}

                        if (last_ev61_time != 0.0 && last_vp_time_18 != 0.0 && (last_vp_time_18 - last_ev61_time) >= -720 && (last_vp_time_18 - last_ev61_time) <= ev_61_dt_max) {h_dt_61_18->Fill(last_vp_time_18 - last_ev61_time);}

                        if (last_ev61_time != 0.0 && last_vp_time_19 != 0.0 && (last_vp_time_19 - last_ev61_time) >= -720 && (last_vp_time_19 - last_ev61_time) <= ev_61_dt_max) {h_dt_61_19->Fill(last_vp_time_19 - last_ev61_time);}

                        if (last_ev61_time != 0.0 && last_vp_time_20 != 0.0 && (last_vp_time_20 - last_ev61_time) >= -720 && (last_vp_time_20 - last_ev61_time) <= ev_61_dt_max) {h_dt_61_20->Fill(last_vp_time_20 - last_ev61_time);}

                        if (last_ev61_time != 0.0 && last_vp_time_21 != 0.0 && (last_vp_time_21 - last_ev61_time) >= -720 && (last_vp_time_21 - last_ev61_time) <= ev_61_dt_max) {h_dt_61_21->Fill(last_vp_time_21 - last_ev61_time);}

                        if (last_ev61_time != 0.0 && last_vp_time_22 != 0.0 && (last_vp_time_22 - last_ev61_time) >= -720 && (last_vp_time_22 - last_ev61_time) <= ev_61_dt_max) {h_dt_61_22->Fill(last_vp_time_22 - last_ev61_time);}

                        if (last_ev61_time != 0.0 && last_vp_time_23 != 0.0 && (last_vp_time_23 - last_ev61_time) >= -720 && (last_vp_time_23 - last_ev61_time) <= ev_61_dt_max) {h_dt_61_23->Fill(last_vp_time_23 - last_ev61_time);}

                        if (last_ev61_time != 0.0 && last_vp_time_top != 0.0 && (last_vp_time_top - last_ev61_time) >= -720 && (last_vp_time_top - last_ev61_time) <= ev_61_dt_max) {h_dt_61_top->Fill(last_vp_time_top - last_ev61_time);}

                        if (last_ev61_time != 0.0 && last_vp_time != 0.0 && (last_vp_time - last_ev61_time) >= -720 && (last_vp_time - last_ev61_time) <= ev_61_dt_max) {h_dt_61_v_int->Fill((last_vp_time - last_ev61_time), (avg_pulse.side_vp_energy + avg_pulse.top_vp_energy));}

                        before_ev61_peak = true;           // Delete this?

                        continue;

                    }

                }

            }

            // }

            // Is this event an Event 61 pulse?

            if (before_ev61_peak && avg_pulse.beam) {           // bt_HL_n_LL &&             // Does before_ev61_peak check make me miss ev61 events?

                // h_int->Fill(avg_pulse.energy);

                last_ev61_time = nsTime + chan_15_start;

                if (last_vp_time_16 != 0.0 && (last_ev61_time - last_vp_time_16) >= -720 && (last_ev61_time - last_vp_time_16) <= ev_61_dt_max) {h_dt_61_16->Fill(last_vp_time_16 - last_ev61_time);}

                if (last_vp_time_17 != 0.0 && (last_ev61_time - last_vp_time_17) >= -720 && (last_ev61_time - last_vp_time_17) <= ev_61_dt_max) {h_dt_61_17->Fill(last_vp_time_17 - last_ev61_time);}

                if (last_vp_time_18 != 0.0 && (last_ev61_time - last_vp_time_18) >= -720 && (last_ev61_time - last_vp_time_18) <= ev_61_dt_max) {h_dt_61_18->Fill(last_vp_time_18 - last_ev61_time);}

                if (last_vp_time_19 != 0.0 && (last_ev61_time - last_vp_time_19) >= -720 && (last_ev61_time - last_vp_time_19) <= ev_61_dt_max) {h_dt_61_19->Fill(last_vp_time_19 - last_ev61_time);}

                if (last_vp_time_20 != 0.0 && (last_ev61_time - last_vp_time_20) >= -720 && (last_ev61_time - last_vp_time_20) <= ev_61_dt_max) {h_dt_61_20->Fill(last_vp_time_20 - last_ev61_time);}

                if (last_vp_time_21 != 0.0 && (last_ev61_time - last_vp_time_21) >= -720 && (last_ev61_time - last_vp_time_21) <= ev_61_dt_max) {h_dt_61_21->Fill(last_vp_time_21 - last_ev61_time);}

                if (last_vp_time_22 != 0.0 && (last_ev61_time - last_vp_time_22) >= -720 && (last_ev61_time - last_vp_time_22) <= ev_61_dt_max) {h_dt_61_22->Fill(last_vp_time_22 - last_ev61_time);}

                if (last_vp_time_23 != 0.0 && (last_ev61_time - last_vp_time_23) >= -720 && (last_ev61_time - last_vp_time_23) <= ev_61_dt_max) {h_dt_61_23->Fill(last_vp_time_23 - last_ev61_time);}

                if (last_vp_time_top != 0.0 && (last_ev61_time - last_vp_time_top) >= -720 && (last_ev61_time - last_vp_time_top) <= ev_61_dt_max) {h_dt_61_top->Fill(last_vp_time_top - last_ev61_time);}

                if (last_vp_time != 0.0 && (last_ev61_time - last_vp_time) >= -720 && (last_ev61_time - last_vp_time) <= ev_61_dt_max) {h_dt_61_v_int->Fill((last_vp_time - last_ev61_time), (avg_pulse.side_vp_energy + avg_pulse.top_vp_energy));}

                before_ev61_peak = false;           // Delete this?
                
                continue;
                
            }

            all_chan_start_vp.clear();
            all_chan_start_vp_iChan.clear();

            num_chan_vp = 0;

            p_int_16 = 0;
            p_int_17 = 0;
            p_int_18 = 0;
            p_int_19 = 0;
            p_int_20 = 0;
            p_int_21 = 0;
            p_int_22 = 0;
            p_int_23 = 0;
            p_int_24 = 0;
            p_int_25 = 0;

            // Ev61Energy = 0.;

        } // Event loop

        fileOut->Close();

        ch15_on = 0;

        integralToPE.clear();

        avg_peak_pos_RMS = getAverage(peak_pos_RMS);

        var_peak_pos_RMS = variance(peak_pos_RMS);

        peak_pos_RMS.clear();

        // std::cout << "\n" << "Completed run " << iRun << "\n";

    } // Run loop

    std::cout << "\n" << "Number of low dt events = " << num_events_low_dt << "\n";

    std::cout << "\n" << "Number of central spike events = " << num_events_low_dt_cent_spike << "\n";

    std::cout << "\n" << "Number of negative low dt events, excluding central spike = " << num_events_low_dt_no_spike_neg << "\n";

    std::cout << "\n" << "Number of positive low dt events, excluding central spike = " << num_events_low_dt_no_spike_pos << "\n";

    std::cout << "\n" << "Number of muons found = " << num_muons_found << "\n";

    /* Fill & plot histograms */

    TCanvas* c_dt_61_16 = new TCanvas("c_dt_61_16", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61_16->SetLogy();
    c_dt_61_16->cd();
    h_dt_61_16->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_16->GetYaxis()->SetTitle("Counts");
    h_dt_61_16->Draw();
    c_dt_61_16->SaveAs("h_dt_61_16.png");
    c_dt_61_16->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_16.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_16.png"));
    h_dt_61_16->Reset();

    TCanvas* c_dt_61_16_simple = new TCanvas("c_dt_61_16_simple", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_16_simple->SetLogy();
    c_dt_61_16_simple->cd();
    h_dt_61_16_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_16_simple->GetYaxis()->SetTitle("Counts");
    h_dt_61_16_simple->Draw();
    h_dt_61_16_simple->SetMaximum(900);
    h_dt_61_16_simple->SetMinimum(600);
    c_dt_61_16_simple->SaveAs("h_dt_61_16_simple.png");
    c_dt_61_16_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_16_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_16_simple.png"));
    h_dt_61_16_simple->Reset();

    TCanvas* c_dt_61_17 = new TCanvas("c_dt_61_17", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61_17->SetLogy();
    c_dt_61_17->cd();
    h_dt_61_17->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_17->GetYaxis()->SetTitle("Counts");
    h_dt_61_17->Draw();
    c_dt_61_17->SaveAs("h_dt_61_17.png");
    c_dt_61_17->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_17.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_17.png"));
    h_dt_61_17->Reset();

    TCanvas* c_dt_61_17_simple = new TCanvas("c_dt_61_17_simple", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_17_simple->SetLogy();
    c_dt_61_17_simple->cd();
    h_dt_61_17_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_17_simple->GetYaxis()->SetTitle("Counts");
    h_dt_61_17_simple->Draw();
    h_dt_61_17_simple->SetMaximum(1200);
    h_dt_61_17_simple->SetMinimum(900);
    c_dt_61_17_simple->SaveAs("h_dt_61_17_simple.png");
    c_dt_61_17_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_17_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_17_simple.png"));
    h_dt_61_17_simple->Reset();

    TCanvas* c_dt_61_18 = new TCanvas("c_dt_61_18", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61_18->SetLogy();
    c_dt_61_18->cd();
    h_dt_61_18->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_18->GetYaxis()->SetTitle("Counts");
    h_dt_61_18->Draw();
    c_dt_61_18->SaveAs("h_dt_61_18.png");
    c_dt_61_18->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_18.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_18.png"));
    h_dt_61_18->Reset();

    TCanvas* c_dt_61_18_simple = new TCanvas("c_dt_61_18_simple", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_18_simple->SetLogy();
    c_dt_61_18_simple->cd();
    h_dt_61_18_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_18_simple->GetYaxis()->SetTitle("Counts");
    h_dt_61_18_simple->Draw();
    h_dt_61_18_simple->SetMaximum(1900);
    h_dt_61_18_simple->SetMinimum(1500);
    c_dt_61_18_simple->SaveAs("h_dt_61_18_simple.png");
    c_dt_61_18_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_18_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_18_simple.png"));
    h_dt_61_18_simple->Reset();

    TCanvas* c_dt_61_19 = new TCanvas("c_dt_61_19", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61_19->SetLogy();
    c_dt_61_19->cd();
    h_dt_61_19->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_19->GetYaxis()->SetTitle("Counts");
    h_dt_61_19->Draw();
    c_dt_61_19->SaveAs("h_dt_61_19.png");
    c_dt_61_19->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_19.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_19.png"));
    h_dt_61_19->Reset();

    TCanvas* c_dt_61_19_simple = new TCanvas("c_dt_61_19_simple", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_19_simple->SetLogy();
    c_dt_61_19_simple->cd();
    h_dt_61_19_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_19_simple->GetYaxis()->SetTitle("Counts");
    h_dt_61_19_simple->Draw();
    h_dt_61_19_simple->SetMaximum(1200);
    h_dt_61_19_simple->SetMinimum(900);
    c_dt_61_19_simple->SaveAs("h_dt_61_19_simple.png");
    c_dt_61_19_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_19_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_19_simple.png"));
    h_dt_61_19_simple->Reset();

    TCanvas* c_dt_61_20 = new TCanvas("c_dt_61_20", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61_20->SetLogy();
    c_dt_61_20->cd();
    h_dt_61_20->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_20->GetYaxis()->SetTitle("Counts");
    h_dt_61_20->Draw();
    c_dt_61_20->SaveAs("h_dt_61_20.png");
    c_dt_61_20->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_20.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_20.png"));
    h_dt_61_20->Reset();

    TCanvas* c_dt_61_20_simple = new TCanvas("c_dt_61_20_simple", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_20_simple->SetLogy();
    c_dt_61_20_simple->cd();
    h_dt_61_20_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_20_simple->GetYaxis()->SetTitle("Counts");
    h_dt_61_20_simple->Draw();
    h_dt_61_20_simple->SetMaximum(1400);
    h_dt_61_20_simple->SetMinimum(1100);
    c_dt_61_20_simple->SaveAs("h_dt_61_20_simple.png");
    c_dt_61_20_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_20_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_20_simple.png"));
    h_dt_61_20_simple->Reset();

    TCanvas* c_dt_61_21 = new TCanvas("c_dt_61_21", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61_21->SetLogy();
    c_dt_61_21->cd();
    h_dt_61_21->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_21->GetYaxis()->SetTitle("Counts");
    h_dt_61_21->Draw();
    c_dt_61_21->SaveAs("h_dt_61_21.png");
    c_dt_61_21->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_21.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_21.png"));
    h_dt_61_21->Reset();

    TCanvas* c_dt_61_21_simple = new TCanvas("c_dt_61_21_simple", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_21_simple->SetLogy();
    c_dt_61_21_simple->cd();
    h_dt_61_21_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_21_simple->GetYaxis()->SetTitle("Counts");
    h_dt_61_21_simple->Draw();
    h_dt_61_21_simple->SetMaximum(1800);
    h_dt_61_21_simple->SetMinimum(1400);
    c_dt_61_21_simple->SaveAs("h_dt_61_21_simple.png");
    c_dt_61_21_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_21_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_21_simple.png"));
    h_dt_61_21_simple->Reset();

    TCanvas* c_dt_61_22 = new TCanvas("c_dt_61_22", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61_22->SetLogy();
    c_dt_61_22->cd();
    h_dt_61_22->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_22->GetYaxis()->SetTitle("Counts");
    h_dt_61_22->Draw();
    c_dt_61_22->SaveAs("h_dt_61_22.png");
    c_dt_61_22->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_22.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_22.png"));
    h_dt_61_22->Reset();

    TCanvas* c_dt_61_22_simple = new TCanvas("c_dt_61_22_simple", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_22_simple->SetLogy();
    c_dt_61_22_simple->cd();
    h_dt_61_22_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_22_simple->GetYaxis()->SetTitle("Counts");
    h_dt_61_22_simple->Draw();
    h_dt_61_22_simple->SetMaximum(2100);
    h_dt_61_22_simple->SetMinimum(1600);
    c_dt_61_22_simple->SaveAs("h_dt_61_22_simple.png");
    c_dt_61_22_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_22_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_22_simple.png"));
    h_dt_61_22_simple->Reset();

    TCanvas* c_dt_61_23 = new TCanvas("c_dt_61_23", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61_23->SetLogy();
    c_dt_61_23->cd();
    h_dt_61_23->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_23->GetYaxis()->SetTitle("Counts");
    h_dt_61_23->Draw();
    c_dt_61_23->SaveAs("h_dt_61_23.png");
    c_dt_61_23->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_23.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_23.png"));
    h_dt_61_23->Reset();

    TCanvas* c_dt_61_23_simple = new TCanvas("c_dt_61_23_simple", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_23_simple->SetLogy();
    c_dt_61_23_simple->cd();
    h_dt_61_23_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_23_simple->GetYaxis()->SetTitle("Counts");
    h_dt_61_23_simple->Draw();
    h_dt_61_23_simple->SetMaximum(1600);
    h_dt_61_23_simple->SetMinimum(1200);
    c_dt_61_23_simple->SaveAs("h_dt_61_23_simple.png");
    c_dt_61_23_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_23_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_23_simple.png"));
    h_dt_61_23_simple->Reset();

    TCanvas* c_dt_61_top = new TCanvas("c_dt_61_top", "Delta-T Between Event 61 and Detector Events", 900, 700);
    c_dt_61_top->SetLogy();
    c_dt_61_top->cd();
    h_dt_61_top->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_top->GetYaxis()->SetTitle("Counts");
    h_dt_61_top->Draw();
    c_dt_61_top->SaveAs("h_dt_61_top.png");
    c_dt_61_top->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_top.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_top.png"));
    h_dt_61_top->Reset();

    TCanvas* c_dt_61_top_simple = new TCanvas("c_dt_61_top_simple", "Delta-T Between Event 61 and Detector Events", 900, 700);
    // c_dt_61_top_simple->SetLogy();
    c_dt_61_top_simple->cd();
    h_dt_61_top_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_top_simple->GetYaxis()->SetTitle("Counts");
    h_dt_61_top_simple->Draw();
    h_dt_61_top_simple->SetMaximum(1900);
    h_dt_61_top_simple->SetMinimum(1400);
    c_dt_61_top_simple->SaveAs("h_dt_61_top_simple.png");
    c_dt_61_top_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_top_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_top_simple.png"));
    h_dt_61_top_simple->Reset();

    TCanvas* c_dt_61_v_int = new TCanvas("c_dt_61_v_int", "Veto Panel Event 61 dt vs Veto Panel Integral Value", 1200, 700);
    // c_dt_61_v_int->SetLogy();
    c_dt_61_v_int->cd();
    h_dt_61_v_int->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_v_int->GetYaxis()->SetTitle("Integral Value (ADC)");
    h_dt_61_v_int->SetMarkerStyle(7);
    h_dt_61_v_int->Draw();
    c_dt_61_v_int->SaveAs("h_dt_61_v_int.png");
    c_dt_61_v_int->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_v_int.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_v_int.png"));
    h_dt_61_v_int->Reset();

    TCanvas* c_dt_61_v_int_simple = new TCanvas("c_dt_61_v_int_simple", "Veto Panel Event 61 dt vs Veto Panel Integral Value", 1200, 700);
    // c_dt_61_v_int_simple->SetLogy();
    c_dt_61_v_int_simple->cd();
    h_dt_61_v_int_simple->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_61_v_int_simple->GetYaxis()->SetTitle("Integral Value (ADC)");
    h_dt_61_v_int_simple->SetMarkerStyle(7);
    h_dt_61_v_int_simple->Draw();
    c_dt_61_v_int_simple->SaveAs("h_dt_61_v_int_simple.png");
    c_dt_61_v_int_simple->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_v_int_simple.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_v_int_simple.png"));
    h_dt_61_v_int_simple->Reset();
    /*
    TCanvas* c_int = new TCanvas("c_int", "Distribution of All Event Integral Values", 1200, 700);
    //c_int->SetLogy();
    c_int->cd();
    h_int->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int->GetYaxis()->SetTitle("Counts");
    h_int->Draw();
    c_int->SaveAs("h_int.png");
    c_int->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int.png"));
    h_int->Reset();

    TCanvas* c_int_HL = new TCanvas("c_int_HL", "Distribution of All High Light Integral Values", 1200, 700);
    c_int_HL->SetLogy();
    c_int_HL->cd();
    h_int_HL->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int_HL->GetYaxis()->SetTitle("Counts");
    h_int_HL->Draw();
    c_int_HL->SaveAs("h_int_HL.png");
    c_int_HL->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_HL.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_HL.png"));
    h_int_HL->Reset();

    TCanvas* c_int_cent_spike = new TCanvas("c_int_cent_spike", "Distribution of All Central Spike Integral Values", 1200, 700);
    //c_int_cent_spike->SetLogy();
    c_int_cent_spike->cd();
    h_int_cent_spike->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int_cent_spike->GetYaxis()->SetTitle("Counts");
    h_int_cent_spike->Draw();
    c_int_cent_spike->SaveAs("h_int_cent_spike.png");
    c_int_cent_spike->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_cent_spike.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_cent_spike.png"));
    h_int_cent_spike->Reset();

    TCanvas* c_int_neg_low_dt = new TCanvas("c_int_neg_low_dt", "Distribution of All Low dt Detector Integral Values", 1200, 700);
    c_int_neg_low_dt->SetLogy();
    c_int_neg_low_dt->cd();
    h_int_neg_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
    h_int_neg_low_dt->GetYaxis()->SetTitle("Counts");
    h_int_neg_low_dt->Draw("same");
    h_int_neg_low_dt->SetLineColor(kRed);
    h_int_pos_low_dt->Draw("same");
    h_int_pos_low_dt->SetLineColor(kGreen);
    TLegend *leg_int_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
    leg_int_neg_low_dt->AddEntry(h_int_pos_low_dt, "Positive Low dt Events", "l");
    leg_int_neg_low_dt->AddEntry(h_int_neg_low_dt, "Negative Low dt Events", "l");
    leg_int_neg_low_dt->Draw();
    c_int_neg_low_dt->SaveAs("h_int_low_dt.png");
    c_int_neg_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_neg_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_neg_low_dt.png"));
    h_int_neg_low_dt->Reset();
    h_int_neg_low_dt->Reset();

    TCanvas* c_svpint_neg_low_dt = new TCanvas("c_svpint_neg_low_dt", "Distribution of All Low dt Side Veto Panel Integral Values", 1200, 700);
    c_svpint_neg_low_dt->SetLogy();
    c_svpint_neg_low_dt->cd();
    h_svpint_neg_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
    h_svpint_neg_low_dt->GetYaxis()->SetTitle("Counts");
    h_svpint_neg_low_dt->Draw("same");
    h_svpint_neg_low_dt->SetLineColor(kRed);
    h_svpint_pos_low_dt->Draw("same");
    h_svpint_pos_low_dt->SetLineColor(kGreen);
    h_svpint_cent_spike->Draw("same");
    h_svpint_cent_spike->SetLineColor(kBlue);
    TLegend *leg_svpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
    leg_svpint_neg_low_dt->AddEntry(h_svpint_cent_spike, "Central Spike Events", "l");
    leg_svpint_neg_low_dt->AddEntry(h_svpint_pos_low_dt, "Positive Low dt Events", "l");
    leg_svpint_neg_low_dt->AddEntry(h_svpint_neg_low_dt, "Negative Low dt Events", "l");
    leg_svpint_neg_low_dt->Draw();
    c_svpint_neg_low_dt->SaveAs("h_svpint_low_dt.png");
    c_svpint_neg_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_svpint_neg_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_svpint_neg_low_dt.png"));
    h_svpint_neg_low_dt->Reset();
    h_svpint_neg_low_dt->Reset();

    TCanvas* c_tvpint_neg_low_dt = new TCanvas("c_tvpint_neg_low_dt", "Distribution of All Low dt Top Veto Panel Integral Values", 1200, 700);
    c_tvpint_neg_low_dt->SetLogy();
    c_tvpint_neg_low_dt->cd();
    h_tvpint_neg_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
    h_tvpint_neg_low_dt->GetYaxis()->SetTitle("Counts");
    h_tvpint_neg_low_dt->Draw("same");
    h_tvpint_neg_low_dt->SetLineColor(kRed);
    h_tvpint_pos_low_dt->Draw("same");
    h_tvpint_pos_low_dt->SetLineColor(kGreen);
    h_tvpint_cent_spike->Draw("same");
    h_tvpint_cent_spike->SetLineColor(kBlue);
    TLegend *leg_tvpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
    leg_tvpint_neg_low_dt->AddEntry(h_tvpint_cent_spike, "Central Spike Events", "l");
    leg_tvpint_neg_low_dt->AddEntry(h_tvpint_pos_low_dt, "Positive Low dt Events", "l");
    leg_tvpint_neg_low_dt->AddEntry(h_tvpint_neg_low_dt, "Negative Low dt Events", "l");
    leg_tvpint_neg_low_dt->Draw();
    c_tvpint_neg_low_dt->SaveAs("h_tvpint_low_dt.png");
    c_tvpint_neg_low_dt->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tvpint_neg_low_dt.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tvpint_neg_low_dt.png"));
    h_tvpint_neg_low_dt->Reset();
    h_tvpint_neg_low_dt->Reset();

    TCanvas* c_dt_v_tmuon = new TCanvas("c_dt_v_tmuon", "Event 61 Detector dt vs Muon dt", 1200, 700);
    // c_dt_v_tmuon->SetLogy();
    c_dt_v_tmuon->cd();
    h_dt_v_tmuon->GetXaxis()->SetTitle("Delta-T (ns)");
    h_dt_v_tmuon->GetYaxis()->SetTitle("Time Since Last Muon (ns)");
    h_dt_v_tmuon->SetMarkerStyle(7);
    h_dt_v_tmuon->Draw();
    c_dt_v_tmuon->SaveAs("h_dt_v_tmuon.png");
    c_dt_v_tmuon->Close();
    // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_tmuon.png"));
    // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_tmuon.png"));
    h_dt_v_tmuon->Reset();
    */
    std::cout << "\n" << "End of code." << "\n";

    // _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG); 
    // _CrtDumpMemoryLeaks();

    return 0;

}