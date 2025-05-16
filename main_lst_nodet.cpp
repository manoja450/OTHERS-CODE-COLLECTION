#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <TCanvas.h>
#include <TSystem.h>
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
// #include "TSpectrum.h"
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
  double top_vp_energy;  /* Energy (integral) of pulse (photo-electrons) in TOP veto panel */
  double all_vp_energy;  /* Energy (integral) of pulse (photo-electrons) in ALL veto panels */
  double last_muon_time; /* Time value of most recent detected muon event */
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
    cout << "run: " << run << " last run: " << last_run << endl;

    int ADCSIZE = 45;
    TH1D *h_wf = new TH1D("h_wf", "Waveform", ADCSIZE, 0, ADCSIZE);

    /* Create vectors for plotting histograms */

    int ev_61_dt_max = pow(10, 5);

    int muon_dt_max = 15 * pow(10, 3);

    int sample_rate = pow(10, 5);

    int take_sample_count = 0;

    int before_cent_count = 0;

    int after_cent_count = 0;

    int det_event_count = 0;

    int cent_event_muon_cut_count = 0;

    int cent_event_else_count = 0;

    int det_event_muon_cut_count = 0;

    int det_event_else_count = 0;

    bool cent_event_found = false;

    bool take_sample = false;

    bool bt_LL_n_HL = false;

    int num_events_low_dt = 0;

    int num_events_low_dt_cent_spike = 0;

    int num_events_low_dt_no_spike_neg = 0;

    int num_events_low_dt_no_spike_pos = 0;

    int num_muons_found = 0;

    double LL_time = 0.0;

    /* Initialize histograms */                 // NUMBER OF BINS SHOULD BE A MULTIPLE OF THE SMALLEST AXIS UNIT!!! (OR LENGTH OF AXIS DIVIDED BY SMALLEST AXIS UNIT?)

    TH1D* h_dt_61 = new TH1D("h_dt_61", "Distribution of Time Separations Between Event 61 and Detector Peaks", 1250, - 2 * pow(10, 5), 2 * pow(10, 5));           // 1 bin = 320 ns

    TH1D* h_dt_61_1 = new TH1D("h_dt_61_1", "Distribution of Time Separations Between Event 61 and Detector Peaks", 625, - 1 * pow(10, 4), 1 * pow(10, 4));         // 1 bin = 32 ns

    TH1D* h_dt_61_2 = new TH1D("h_dt_61_2", "Distribution of Time Separations Between Event 61 and Detector Peaks", 1250, - 2 * pow(10, 4), 2 * pow(10, 4));        // 1 bin = 32 ns

    TH1D* h_dt_61_3 = new TH1D("h_dt_61_3", "Distribution of Time Separations Between Event 61 and Detector Peaks", 3125, - 5 * pow(10, 4), 5 * pow(10, 4));        // 1 bin = 32 ns

    TH1D* h_dt_61_4 = new TH1D("h_dt_61_4", "Distribution of Time Separations Between Event 61 and Detector Peaks", 1250, - 10 * pow(10, 4), 10 * pow(10, 4));      // 1 bin = 160 ns

    TH1D* h_detint = new TH1D("h_detint", "Distribution of All Ev61Det Detector Integral Values", 200, 0, 2000);

    TH1D* h_detint_low_dt = new TH1D("h_detint_low_dt", "Distribution of Low dt Ev61Det Detector Integral Values", 80, 0, 800);

    TH2D* h_dt_v_detint = new TH2D("h_dt_v_detint", "Delta-T vs Detector Integral Value", 1250, - ev_61_dt_max, ev_61_dt_max, 200, 0, 2000);                        // 1 bin = 160 ns

    TH1D* h_int = new TH1D("h_int", "Distribution of All Event Integral Values Between High Light and Low Light LED Events", 200, 0, 2000);

    TH1D* h_int_cent_spike = new TH1D("h_int_cent_spike", "Distribution of All Central Spike Detector Integral Values", 200, 0, 2000);

    TH1D* h_int_pos_low_dt = new TH1D("h_int_pos_low_dt", "Distribution of Low dt Detector Integral Values", 200, 0, 2000);

    TH1D* h_int_neg_low_dt = new TH1D("h_int_neg_low_dt", "Distribution of Low dt Detector Integral Values", 200, 0, 2000);

    TH1D* h_svpint_cent_spike = new TH1D("h_svpint_cent_spike", "Distribution of All Central Spike Side Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_svpint_pos_low_dt = new TH1D("h_svpint_pos_low_dt", "Distribution of Low dt Side Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_svpint_neg_low_dt = new TH1D("h_svpint_neg_low_dt", "Distribution of Low dt Side Veto Panel Integral Values", 240, -400, 2000);

    TH1D* h_tvpint_cent_spike = new TH1D("h_tvpint_cent_spike", "Distribution of All Central Spike Top Veto Panel Integral Values", 60, -200, 400);

    TH1D* h_tvpint_pos_low_dt = new TH1D("h_tvpint_pos_low_dt", "Distribution of Low dt Top Veto Panel Integral Values", 60, -200, 400);

    TH1D* h_tvpint_neg_low_dt = new TH1D("h_tvpint_neg_low_dt", "Distribution of Low dt Top Veto Panel Integral Values", 60, -200, 400);

    TH1D* h_avpint_cent_spike = new TH1D("h_avpint_cent_spike", "Distribution of All Central Spike Veto Panel Integral Values", 60, -400, 2000);

    TH1D* h_avpint_pos_low_dt = new TH1D("h_avpint_pos_low_dt", "Distribution of Low dt Veto Panel Integral Values", 60, -400, 2000);

    TH1D* h_avpint_neg_low_dt = new TH1D("h_avpint_neg_low_dt", "Distribution of Low dt Veto Panel Integral Values", 60, -400, 2000);

    TH1D* h_avpint_low_dt = new TH1D("h_avpint_low_dt", "Distribution of Low dt Veto Panel Integral Values", 60, -400, 2000);

    TH1D* h_int_HL = new TH1D("h_int_HL", "Distribution of All High Light LED Integral Values", 5500, -500, 5000);

    TH2D* h_dt_v_tmuon = new TH2D("h_dt_v_tmuon", "Delta-T vs Most Recent Muon Time Difference", 1000, - 2 * ev_61_dt_max, 2 * ev_61_dt_max, 1000, 0, 2 * pow(10, 7));

    TH1D* h_allplsvec_size = new TH1D("h_allplsvec_size", "Size of Vector of All Events Between Low Light and High Light Events", 10, 0, 10);

    TH1D* h_inttrigvec_size = new TH1D("h_inttrigvec_size", "Size of Vector of Internally Triggered Events Between Low Light and High Light Events", 10, 0, 10);

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

        std::vector<struct pulse> pulses_vec;

        std::vector<struct pulse> int_trig_vec;

        std::vector<double> chan_lengths;

        std::vector<double> peak_pos_RMS;

        std::vector<double> chan_start_no_outliers;

        bool record_pulses = false;

        double avg_peak_pos_RMS = 0.0;

        double var_peak_pos_RMS = 0.0;

        int ch15_on = 0;

        double last_muon_time = 0.0;

        double last_take_sample_time = 0.0;
    
        // Can use either "t->GetEntries()" or "numEntries" or number
        for (int iEnt = 0; iEnt < t->GetEntries(); iEnt++) {

            std::cout << "\n" << "Processing event " << iEnt + 1 << " of " << t->GetEntries() << "\n";

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

            bool all_chan_beam = false;

            bool top_vp_event = false;

            bool pulse_at_end = false;

            int pulse_at_end_count = 0;

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

                if (iChan == 15) {          // Move this to line 641?

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
                    /*
                    if (iChan >= 16 && !peak_in_chan_veto && iBinContent >= PULSE_THRESHOLD) {

                        peak_in_chan_veto = true;

                        num_chan_veto += 1;

                    }
                    */
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

                        // Check if pulse is at end of waveform
                        if (iChan <= 11 && iBin == ADCSIZE && iBinContent > 100) {pulse_at_end_count += 1; if (pulse_at_end_count >= 10) {pulse_at_end = true;}}            // std::cout << "\n" << "pulse_at_end_count = " << pulse_at_end_count << ", Run, Event, triggerBits = " << iRun << " " << iEnt + 1 << " " << triggerBits << "\n";}}

                        // if (iChan <= 11 && iRun == 7894 && iEnt + 1 == 389830) {std::cout << "\n" << "Run 7894, Event 389830: pulse_at_end_count = " << pulse_at_end_count << "; iChan = " << iChan << "; iBin = " << iBin << "; iBinContent = " << iBinContent << "\n";}

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

                                all_chan_start.push_back(p.start);
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

                            }

                            // Clear current pulse variables to look for new pulse
                            peak = 0.;
                            peakBin = 0;
                            pulseEnergy = 0.;
                            thresholdBin = 0;
                            onPulse = false;
                        }

                        if (iBin == ADCSIZE) {
                            
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

            // std::cout << "\n" << "Run Number = " << iRun << ", Event Number = " << iEnt + 1 << ", triggerBits = " << triggerBits << "\n";

            // if (all_chan_beam) {std::cout << "\n" << "Event 61 event" << "\n";}

            // if (triggerBits == 4) {std::cout << "\n" << "Min bias LED event (tB == 4)" << "\n";}

            // else if (triggerBits == 8) {std::cout << "\n" << "High light LED event (tB == 8)" << "\n";}

            // else if (triggerBits == 16) {std::cout << "\n" << "Low light LED event (tB == 16)" << "\n";}

            double pulse_start_time = mostFrequent(all_chan_start);

            double vp_start_time = nsTime + mostFrequent(all_chan_start_vp);

            if (num_chan_over_1phe >= 10) {int pulse_start_index = mostFrequentIndex(all_chan_start);}

            // if (pulse_start_time == 0) {continue;}

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
            avg_pulse.all_vp_energy = avg_pulse.side_vp_energy + avg_pulse.top_vp_energy;

            for (int iPeak = 0; iPeak < all_chan_start.size(); iPeak++) {

                if (all_chan_start[iPeak] < (pulse_start_time + 10 * 16) && all_chan_start[iPeak] > (pulse_start_time - 10 * 16)) {chan_start_no_outliers.push_back(all_chan_start[iPeak]);}

                else {

                    // Check if pulse outside this range occurs more than once

                    for (int jPeak = 0; jPeak < all_chan_start.size(); jPeak++) {

                        if (iPeak == jPeak) {continue;}

                        else if (all_chan_start[iPeak] < (all_chan_start[jPeak] + 1 * 16) && all_chan_start[iPeak] > (all_chan_start[jPeak] - 1 * 16)) {chan_start_no_outliers.push_back(all_chan_start[iPeak]);}

                    }

                }
                
            }

            double var_val = variance(chan_start_no_outliers);

            // h_var->Fill(var_val);

            if (var_val < 5 * 16) {avg_pulse.single = true;}            // Pulses largely not consistent, need to weed out outlier and secondary pulses

            else {avg_pulse.single = false;}

            // Is this event a muon? (Looking at events with peak at the very end of the waveform)

            if (pulse_at_end && avg_pulse.energy > 20 / 2 && (p_int_16 > 750 / 2 || p_int_17 > 950 / 2 || p_int_18 > 1200 / 2 || p_int_19 > 1375 / 2 || p_int_20 > 525 / 2 || p_int_21 > 700 / 2 || p_int_22 > 700 / 2 || p_int_23 > 500 / 2 || avg_pulse.top_vp_energy > 450 / 2)) {

                avg_pulse.last_muon_time = vp_start_time;

                num_muons_found += 1;

            }

            // Is this event a muon? (Looking at "normal" events)

            else if (avg_pulse.energy > 20 && (p_int_16 > 750 || p_int_17 > 950 || p_int_18 > 1200 || p_int_19 > 1375 || p_int_20 > 525 || p_int_21 > 700 || p_int_22 > 700 || p_int_23 > 500 || avg_pulse.top_vp_energy > 450)) {

                avg_pulse.last_muon_time = vp_start_time;

                num_muons_found += 1;

            }

            /*
            
            New Rand-Det dt plot code loop
            
            */

            // After every low light LED event, enable recording of vector pulses

            if (avg_pulse.trigger == 16) {record_pulses = true; LL_time = avg_pulse.start; pulses_vec.clear(); int_trig_vec.clear();}
            
            // After every high light LED event, disable recording of vector pulses
            
            else if (avg_pulse.trigger == 8) {

                // Scan through time in 100 us intervals, from time of last low light event to time of current high light event

                int iTime = ceil(LL_time);

                while (iTime < floor(avg_pulse.start)) {
                    
                    iTime = iTime + sample_rate;

                    // Iterate over all events in compiled vector

                    for (int iVec = 0; iVec < int_trig_vec.size(); iVec++) {

                        take_sample_count += 1;         // Debugging counter

                    }
                    
                }
                
                record_pulses = false; h_allplsvec_size->Fill(pulses_vec.size()); h_inttrigvec_size->Fill(int_trig_vec.size()); pulses_vec.clear(); int_trig_vec.clear();
                
            }

            // If recording of vector pulses enabled, put average pulse variables into a vector 

            if (record_pulses && avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {pulses_vec.push_back(avg_pulse);}

            // Is this event an internally triggered detector pulse? All if statements are just cuts on which events I record

            if (record_pulses && avg_pulse.number >= 10 && avg_pulse.energy >= 100) {
                
                if (avg_pulse.trigger != 0 && avg_pulse.trigger != 4 && avg_pulse.trigger != 8 && avg_pulse.trigger != 16) {

                    if (p_int_16 <= 750 && p_int_17 <= 950 && p_int_18 <= 1200 && p_int_19 <= 1375 && p_int_20 <= 525 && p_int_21 <= 700 && p_int_22 <= 700 && p_int_23 <= 500 && avg_pulse.top_vp_energy <= 450) {

                        if (avg_pulse.start - avg_pulse.last_muon_time > muon_dt_max) {

                            det_event_count += 1;

                            int_trig_vec.push_back(avg_pulse);

                        }

                    }

                }

            }

            h_int->Fill(avg_pulse.energy);

            // Reset variables at end of loop

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
            num_chan = 0;
            num_chan_vp = 0;
            num_chan_simple = 0;
            num_chan_over_1phe = 0;

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
            chan_start_no_outliers.clear();
            side_veto_panel_energy.clear();
            top_veto_panel_energy.clear();
            chan_lengths.clear();
            all_chan_start_vp.clear();

        } // Event loop

        fileOut->Close();

        ch15_on = 0;

        integralToPE.clear();

        avg_peak_pos_RMS = getAverage(peak_pos_RMS);

        var_peak_pos_RMS = variance(peak_pos_RMS);

        peak_pos_RMS.clear();

        // std::cout << "\n" << "Completed run " << iRun << "\n";

    } // Run loop

    if (run == 7894) {

        std::cout << "\n" << "Number of low dt events, Runs 7894-7918 = " << num_events_low_dt << "\n";

        std::cout << "\n" << "Number of central spike events, Runs 7894-7918 = " << num_events_low_dt_cent_spike << "\n";

        std::cout << "\n" << "Number of negative low dt events, excluding central spike, Runs 7894-7918 = " << num_events_low_dt_no_spike_neg << "\n";

        std::cout << "\n" << "Number of positive low dt events, excluding central spike, Runs 7894-7918 = " << num_events_low_dt_no_spike_pos << "\n";

        std::cout << "\n" << "Number of muons found, Runs 7894-7918 = " << num_muons_found << "\n";

        std::cout << "\n" << "Number of central samples taken, Runs 7894-7918 = " << take_sample_count << "\n";

        std::cout << "\n" << "Number of detector events found, Runs 7894-7918 = " << det_event_count << "\n";

        std::cout << "\n" << "Number of histogram events BEFORE central sample, Runs 7894-7918 = " << before_cent_count << "\n";

        std::cout << "\n" << "Number of histogram events AFTER central sample, Runs 7894-7918 = " << after_cent_count << "\n" << "\n";

    }

    if (run == 11608) {

        std::cout << "\n" << "Number of low dt events, Runs 11608-11632 = " << num_events_low_dt << "\n";

        std::cout << "\n" << "Number of central spike events, Runs 11608-11632 = " << num_events_low_dt_cent_spike << "\n";

        std::cout << "\n" << "Number of negative low dt events, excluding central spike, Runs 11608-11632 = " << num_events_low_dt_no_spike_neg << "\n";

        std::cout << "\n" << "Number of positive low dt events, excluding central spike, Runs 11608-11632 = " << num_events_low_dt_no_spike_pos << "\n";

        std::cout << "\n" << "Number of muons found, Runs 11608-11632 = " << num_muons_found << "\n";

        std::cout << "\n" << "Number of central samples taken, Runs 11608-11632 = " << take_sample_count << "\n";

        std::cout << "\n" << "Number of detector events found, Runs 11608-11632 = " << det_event_count << "\n";

        std::cout << "\n" << "Number of histogram events BEFORE central sample, Runs 11608-11632 = " << before_cent_count << "\n";

        std::cout << "\n" << "Number of histogram events AFTER central sample, Runs 11608-11632 = " << after_cent_count << "\n" << "\n";

    }

    if (run == 12636) {

        std::cout << "\n" << "Number of low dt events, Runs 12636-12660 = " << num_events_low_dt << "\n";

        std::cout << "\n" << "Number of central spike events, Runs 12636-12660 = " << num_events_low_dt_cent_spike << "\n";

        std::cout << "\n" << "Number of negative low dt events, excluding central spike, Runs 12636-12660 = " << num_events_low_dt_no_spike_neg << "\n";

        std::cout << "\n" << "Number of positive low dt events, excluding central spike, Runs 12636-12660 = " << num_events_low_dt_no_spike_pos << "\n";

        std::cout << "\n" << "Number of muons found, Runs 12636-12660 = " << num_muons_found << "\n";

        std::cout << "\n" << "Number of central samples taken, Runs 12636-12660 = " << take_sample_count << "\n";

        std::cout << "\n" << "Number of detector events found, Runs 12636-12660 = " << det_event_count << "\n";

        std::cout << "\n" << "Number of histogram events BEFORE central sample, Runs 12636-12660 = " << before_cent_count << "\n";

        std::cout << "\n" << "Number of histogram events AFTER central sample, Runs 12636-12660 = " << after_cent_count << "\n" << "\n";

    }

    /* Fill & plot histograms */

    if (run == 7894) {

        TCanvas* c_dt_61 = new TCanvas("c_dt_61", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61->SetLogy();
        c_dt_61->cd();
        h_dt_61->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61->GetYaxis()->SetTitle("Counts");
        h_dt_61->Draw();
        c_dt_61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61.png");
        c_dt_61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61.png"));
        h_dt_61->Reset();

        TCanvas* c_dt_61_1 = new TCanvas("c_dt_61_1", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_1->SetLogy();
        c_dt_61_1->cd();
        h_dt_61_1->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_1->GetYaxis()->SetTitle("Counts");
        h_dt_61_1->Draw();
        c_dt_61_1->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_1.png");
        c_dt_61_1->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_1.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_1.png"));
        h_dt_61_1->Reset();

        TCanvas* c_dt_61_2 = new TCanvas("c_dt_61_2", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_2->SetLogy();
        c_dt_61_2->cd();
        h_dt_61_2->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_2->GetYaxis()->SetTitle("Counts");
        h_dt_61_2->Draw();
        c_dt_61_2->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_2.png");
        c_dt_61_2->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_2.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_2.png"));
        h_dt_61_2->Reset();

        TCanvas* c_dt_61_3 = new TCanvas("c_dt_61_3", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_3->SetLogy();
        c_dt_61_3->cd();
        h_dt_61_3->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_3->GetYaxis()->SetTitle("Counts");
        h_dt_61_3->Draw();
        c_dt_61_3->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_3.png");
        c_dt_61_3->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_3.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_3.png"));
        h_dt_61_3->Reset();

        TCanvas* c_dt_61_4 = new TCanvas("c_dt_61_4", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_4->SetLogy();
        c_dt_61_4->cd();
        h_dt_61_4->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_4->GetYaxis()->SetTitle("Counts");
        h_dt_61_4->Draw();
        c_dt_61_4->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_61_4.png");
        c_dt_61_4->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_4.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_4.png"));
        h_dt_61_4->Reset();

        TCanvas* c_detint = new TCanvas("c_detint", "Distribution of All Ev61Det Integral Values", 1200, 700);
        //c_detint->SetLogy();
        c_detint->cd();
        h_detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint->GetYaxis()->SetTitle("Counts");
        h_detint->Draw();
        c_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint.png");
        c_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint.png"));
        h_detint->Reset();

        TCanvas* c_detint_low_dt = new TCanvas("c_detint_low_dt", "Distribution of Low dt Ev61Det Integral Values", 1200, 700);
        //c_detint_low_dt->SetLogy();
        c_detint_low_dt->cd();
        h_detint_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint_low_dt->GetYaxis()->SetTitle("Counts");
        h_detint_low_dt->Draw();
        c_detint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_detint_low_dt.png");
        c_detint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_low_dt.png"));
        h_detint_low_dt->Reset();

        TCanvas* c_dt_v_detint = new TCanvas("c_dt_v_detint", "Event 61 Detector dt vs Integral Value", 1200, 700);
        // c_dt_v_detint->SetLogy();
        c_dt_v_detint->cd();
        h_dt_v_detint->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_detint->GetYaxis()->SetTitle("Integral (Ph.e.)");
        h_dt_v_detint->SetMarkerStyle(7);
        h_dt_v_detint->Draw();
        c_dt_v_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_v_detint.png");
        c_dt_v_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_detint.png"));
        h_dt_v_detint->Reset();

        TCanvas* c_int = new TCanvas("c_int", "Distribution of All Event Integral Values", 1200, 700);
        c_int->SetLogy();
        c_int->cd();
        h_int->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int->GetYaxis()->SetTitle("Counts");
        h_int->Draw();
        c_int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_int.png");
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
        c_int_HL->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_int_HL.png");
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
        c_int_cent_spike->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_int_cent_spike.png");
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
        c_int_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_int_low_dt.png");
        c_int_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_neg_low_dt.png"));
        h_int_pos_low_dt->Reset();
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
        c_svpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_svpint_low_dt.png");
        c_svpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_svpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_svpint_neg_low_dt.png"));
        h_svpint_pos_low_dt->Reset();
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
        c_tvpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_tvpint_low_dt.png");
        c_tvpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tvpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tvpint_neg_low_dt.png"));
        h_tvpint_pos_low_dt->Reset();
        h_tvpint_neg_low_dt->Reset();

        TCanvas* c_avpint_neg_low_dt = new TCanvas("c_avpint_neg_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        // c_avpint_neg_low_dt->SetLogy();
        c_avpint_neg_low_dt->cd();
        TF1* func1 = new TF1("func1", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -200, 200);
        func1->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func1->SetParameters(1, 10, 0, 50);
        h_avpint_pos_low_dt->Fit("func1", "R");
        h_avpint_neg_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_neg_low_dt->GetYaxis()->SetTitle("Counts");
        // h_avpint_neg_low_dt->GetYaxis()->SetRange(0, 10);
        h_avpint_neg_low_dt->Draw();
        h_avpint_neg_low_dt->SetLineColor(kRed);
        h_avpint_pos_low_dt->Draw("same");
        h_avpint_pos_low_dt->SetLineColor(kGreen);
        h_avpint_cent_spike->Draw("same");
        h_avpint_cent_spike->SetLineColor(kBlue);
        func1->Draw("same");
        func1->SetLineColor(1);
        TLegend *leg_avpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_avpint_neg_low_dt->AddEntry(h_avpint_cent_spike, "Central Spike Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_avpint_neg_low_dt->Draw();
        c_avpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_avpint_low_dt.png");
        c_avpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_neg_low_dt.png"));
        h_avpint_pos_low_dt->Reset();
        h_avpint_neg_low_dt->Reset();

        TCanvas* c_avpint_low_dt = new TCanvas("c_avpint_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        // c_avpint_low_dt->SetLogy();
        c_avpint_low_dt->cd();
        TF1* func2 = new TF1("func2", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -300, 300);
        func2->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func2->SetParameters(1, 10, 0, 50);
        h_avpint_low_dt->Fit("func2", "R");
        h_avpint_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_low_dt->GetYaxis()->SetTitle("Counts");
        // h_avpint_low_dt->GetYaxis()->SetRange(0, 50);
        h_avpint_low_dt->Draw();
        func2->Draw("same");
        func2->SetLineColor(2);
        c_avpint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_avpint_low_dt_1c.png");
        c_avpint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_low_dt.png"));
        h_avpint_low_dt->Reset();

        TCanvas* c_dt_v_tmuon = new TCanvas("c_dt_v_tmuon", "Event 61 Detector dt vs Muon dt", 1200, 700);
        // c_dt_v_tmuon->SetLogy();
        c_dt_v_tmuon->cd();
        h_dt_v_tmuon->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_tmuon->GetYaxis()->SetTitle("Time Since Last Muon (ns)");
        h_dt_v_tmuon->SetMarkerStyle(7);
        h_dt_v_tmuon->Draw();
        c_dt_v_tmuon->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_dt_v_tmuon.png");
        c_dt_v_tmuon->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_tmuon.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_tmuon.png"));
        h_dt_v_tmuon->Reset();

        TCanvas* c_allplsvec_size = new TCanvas("c_allplsvec_size", "Number of Events b/t LL and HL", 1200, 700);
        c_allplsvec_size->SetLogy();
        c_allplsvec_size->cd();
        h_allplsvec_size->GetXaxis()->SetTitle("Number of Events");
        h_allplsvec_size->GetYaxis()->SetTitle("Counts");
        h_allplsvec_size->Draw();
        c_allplsvec_size->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_allplsvec_size.png");
        c_allplsvec_size->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_allplsvec_size.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_allplsvec_size.png"));
        h_allplsvec_size->Reset();

        TCanvas* c_inttrigvec_size = new TCanvas("c_inttrigvec_size", "Number of Int. Trig. Events b/t LL and HL", 1200, 700);
        c_inttrigvec_size->SetLogy();
        c_inttrigvec_size->cd();
        h_inttrigvec_size->GetXaxis()->SetTitle("Number of Events");
        h_inttrigvec_size->GetYaxis()->SetTitle("Counts");
        h_inttrigvec_size->Draw();
        c_inttrigvec_size->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_inttrigvec_size.png");
        c_inttrigvec_size->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_inttrigvec_size.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_inttrigvec_size.png"));
        h_inttrigvec_size->Reset();

    }

    if (run == 11608) {

        TCanvas* c_dt_61 = new TCanvas("c_dt_61", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61->SetLogy();
        c_dt_61->cd();
        h_dt_61->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61->GetYaxis()->SetTitle("Counts");
        h_dt_61->Draw();
        c_dt_61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61.png");
        c_dt_61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61.png"));
        h_dt_61->Reset();

        TCanvas* c_dt_61_1 = new TCanvas("c_dt_61_1", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_1->SetLogy();
        c_dt_61_1->cd();
        h_dt_61_1->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_1->GetYaxis()->SetTitle("Counts");
        h_dt_61_1->Draw();
        c_dt_61_1->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_1.png");
        c_dt_61_1->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_1.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_1.png"));
        h_dt_61_1->Reset();

        TCanvas* c_dt_61_2 = new TCanvas("c_dt_61_2", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_2->SetLogy();
        c_dt_61_2->cd();
        h_dt_61_2->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_2->GetYaxis()->SetTitle("Counts");
        h_dt_61_2->Draw();
        c_dt_61_2->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_2.png");
        c_dt_61_2->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_2.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_2.png"));
        h_dt_61_2->Reset();

        TCanvas* c_dt_61_3 = new TCanvas("c_dt_61_3", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_3->SetLogy();
        c_dt_61_3->cd();
        h_dt_61_3->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_3->GetYaxis()->SetTitle("Counts");
        h_dt_61_3->Draw();
        c_dt_61_3->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_3.png");
        c_dt_61_3->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_3.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_3.png"));
        h_dt_61_3->Reset();

        TCanvas* c_dt_61_4 = new TCanvas("c_dt_61_4", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_4->SetLogy();
        c_dt_61_4->cd();
        h_dt_61_4->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_4->GetYaxis()->SetTitle("Counts");
        h_dt_61_4->Draw();
        c_dt_61_4->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_61_4.png");
        c_dt_61_4->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_4.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_4.png"));
        h_dt_61_4->Reset();

        TCanvas* c_detint = new TCanvas("c_detint", "Distribution of All Ev61Det Integral Values", 1200, 700);
        //c_detint->SetLogy();
        c_detint->cd();
        h_detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint->GetYaxis()->SetTitle("Counts");
        h_detint->Draw();
        c_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_detint.png");
        c_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint.png"));
        h_detint->Reset();

        TCanvas* c_detint_low_dt = new TCanvas("c_detint_low_dt", "Distribution of Low dt Ev61Det Integral Values", 1200, 700);
        //c_detint_low_dt->SetLogy();
        c_detint_low_dt->cd();
        h_detint_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint_low_dt->GetYaxis()->SetTitle("Counts");
        h_detint_low_dt->Draw();
        c_detint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_detint_low_dt.png");
        c_detint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_low_dt.png"));
        h_detint_low_dt->Reset();

        TCanvas* c_dt_v_detint = new TCanvas("c_dt_v_detint", "Event 61 Detector dt vs Integral Value", 1200, 700);
        // c_dt_v_detint->SetLogy();
        c_dt_v_detint->cd();
        h_dt_v_detint->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_detint->GetYaxis()->SetTitle("Integral (Ph.e.)");
        h_dt_v_detint->SetMarkerStyle(7);
        h_dt_v_detint->Draw();
        c_dt_v_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_v_detint.png");
        c_dt_v_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_detint.png"));
        h_dt_v_detint->Reset();

        TCanvas* c_int = new TCanvas("c_int", "Distribution of All Event Integral Values", 1200, 700);
        c_int->SetLogy();
        c_int->cd();
        h_int->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int->GetYaxis()->SetTitle("Counts");
        h_int->Draw();
        c_int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_int.png");
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
        c_int_HL->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_int_HL.png");
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
        c_int_cent_spike->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_int_cent_spike.png");
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
        c_int_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_int_low_dt.png");
        c_int_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_neg_low_dt.png"));
        h_int_pos_low_dt->Reset();
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
        c_svpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_svpint_low_dt.png");
        c_svpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_svpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_svpint_neg_low_dt.png"));
        h_svpint_pos_low_dt->Reset();
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
        c_tvpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_tvpint_low_dt.png");
        c_tvpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tvpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tvpint_neg_low_dt.png"));
        h_tvpint_pos_low_dt->Reset();
        h_tvpint_neg_low_dt->Reset();

        TCanvas* c_avpint_neg_low_dt = new TCanvas("c_avpint_neg_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        // c_avpint_neg_low_dt->SetLogy();
        c_avpint_neg_low_dt->cd();
        TF1* func1 = new TF1("func1", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -200, 200);
        func1->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func1->SetParameters(1, 10, 0, 50);
        h_avpint_pos_low_dt->Fit("func1", "R");
        h_avpint_neg_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_neg_low_dt->GetYaxis()->SetTitle("Counts");
        // h_avpint_neg_low_dt->GetYaxis()->SetRange(0, 10);
        h_avpint_neg_low_dt->Draw();
        h_avpint_neg_low_dt->SetLineColor(kRed);
        h_avpint_pos_low_dt->Draw("same");
        h_avpint_pos_low_dt->SetLineColor(kGreen);
        h_avpint_cent_spike->Draw("same");
        h_avpint_cent_spike->SetLineColor(kBlue);
        func1->Draw("same");
        func1->SetLineColor(1);
        TLegend *leg_avpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_avpint_neg_low_dt->AddEntry(h_avpint_cent_spike, "Central Spike Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_avpint_neg_low_dt->Draw();
        c_avpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_avpint_low_dt.png");
        c_avpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_neg_low_dt.png"));
        h_avpint_pos_low_dt->Reset();
        h_avpint_neg_low_dt->Reset();

        TCanvas* c_avpint_low_dt = new TCanvas("c_avpint_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        // c_avpint_low_dt->SetLogy();
        c_avpint_low_dt->cd();
        TF1* func2 = new TF1("func2", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -300, 300);
        func2->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func2->SetParameters(1, 10, 0, 50);
        h_avpint_low_dt->Fit("func2", "R");
        h_avpint_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_low_dt->GetYaxis()->SetTitle("Counts");
        // h_avpint_low_dt->GetYaxis()->SetRange(0, 50);
        h_avpint_low_dt->Draw();
        func2->Draw("same");
        func2->SetLineColor(2);
        c_avpint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_avpint_low_dt_1c.png");
        c_avpint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_low_dt.png"));
        h_avpint_low_dt->Reset();

        TCanvas* c_dt_v_tmuon = new TCanvas("c_dt_v_tmuon", "Event 61 Detector dt vs Muon dt", 1200, 700);
        // c_dt_v_tmuon->SetLogy();
        c_dt_v_tmuon->cd();
        h_dt_v_tmuon->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_tmuon->GetYaxis()->SetTitle("Time Since Last Muon (ns)");
        h_dt_v_tmuon->SetMarkerStyle(7);
        h_dt_v_tmuon->Draw();
        c_dt_v_tmuon->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots11608/h_dt_v_tmuon.png");
        c_dt_v_tmuon->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_tmuon.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_tmuon.png"));
        h_dt_v_tmuon->Reset();

        TCanvas* c_allplsvec_size = new TCanvas("c_allplsvec_size", "Number of Events b/t LL and HL", 1200, 700);
        c_allplsvec_size->SetLogy();
        c_allplsvec_size->cd();
        h_allplsvec_size->GetXaxis()->SetTitle("Number of Events");
        h_allplsvec_size->GetYaxis()->SetTitle("Counts");
        h_allplsvec_size->Draw();
        c_allplsvec_size->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_allplsvec_size.png");
        c_allplsvec_size->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_allplsvec_size.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_allplsvec_size.png"));
        h_allplsvec_size->Reset();

        TCanvas* c_inttrigvec_size = new TCanvas("c_inttrigvec_size", "Number of Int. Trig. Events b/t LL and HL", 1200, 700);
        c_inttrigvec_size->SetLogy();
        c_inttrigvec_size->cd();
        h_inttrigvec_size->GetXaxis()->SetTitle("Number of Events");
        h_inttrigvec_size->GetYaxis()->SetTitle("Counts");
        h_inttrigvec_size->Draw();
        c_inttrigvec_size->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_inttrigvec_size.png");
        c_inttrigvec_size->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_inttrigvec_size.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_inttrigvec_size.png"));
        h_inttrigvec_size->Reset();

    }

    if (run == 12636) {

        TCanvas* c_dt_61 = new TCanvas("c_dt_61", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61->SetLogy();
        c_dt_61->cd();
        h_dt_61->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61->GetYaxis()->SetTitle("Counts");
        h_dt_61->Draw();
        c_dt_61->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61.png");
        c_dt_61->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61.png"));
        h_dt_61->Reset();

        TCanvas* c_dt_61_1 = new TCanvas("c_dt_61_1", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_1->SetLogy();
        c_dt_61_1->cd();
        h_dt_61_1->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_1->GetYaxis()->SetTitle("Counts");
        h_dt_61_1->Draw();
        c_dt_61_1->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_1.png");
        c_dt_61_1->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_1.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_1.png"));
        h_dt_61_1->Reset();

        TCanvas* c_dt_61_2 = new TCanvas("c_dt_61_2", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_2->SetLogy();
        c_dt_61_2->cd();
        h_dt_61_2->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_2->GetYaxis()->SetTitle("Counts");
        h_dt_61_2->Draw();
        c_dt_61_2->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_2.png");
        c_dt_61_2->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_2.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_2.png"));
        h_dt_61_2->Reset();

        TCanvas* c_dt_61_3 = new TCanvas("c_dt_61_3", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_3->SetLogy();
        c_dt_61_3->cd();
        h_dt_61_3->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_3->GetYaxis()->SetTitle("Counts");
        h_dt_61_3->Draw();
        c_dt_61_3->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_3.png");
        c_dt_61_3->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_3.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_3.png"));
        h_dt_61_3->Reset();

        TCanvas* c_dt_61_4 = new TCanvas("c_dt_61_4", "Delta-T Between Event 61 and Detector Events", 900, 700);
        // c_dt_61_4->SetLogy();
        c_dt_61_4->cd();
        h_dt_61_4->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_61_4->GetYaxis()->SetTitle("Counts");
        h_dt_61_4->Draw();
        c_dt_61_4->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_61_4.png");
        c_dt_61_4->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_61_4.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_61_4.png"));
        h_dt_61_4->Reset();

        TCanvas* c_detint = new TCanvas("c_detint", "Distribution of All Ev61Det Integral Values", 1200, 700);
        //c_detint->SetLogy();
        c_detint->cd();
        h_detint->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint->GetYaxis()->SetTitle("Counts");
        h_detint->Draw();
        c_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint.png");
        c_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint.png"));
        h_detint->Reset();

        TCanvas* c_detint_low_dt = new TCanvas("c_detint_low_dt", "Distribution of Low dt Ev61Det Integral Values", 1200, 700);
        //c_detint_low_dt->SetLogy();
        c_detint_low_dt->cd();
        h_detint_low_dt->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_detint_low_dt->GetYaxis()->SetTitle("Counts");
        h_detint_low_dt->Draw();
        c_detint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_detint_low_dt.png");
        c_detint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_detint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_detint_low_dt.png"));
        h_detint_low_dt->Reset();

        TCanvas* c_dt_v_detint = new TCanvas("c_dt_v_detint", "Event 61 Detector dt vs Integral Value", 1200, 700);
        // c_dt_v_detint->SetLogy();
        c_dt_v_detint->cd();
        h_dt_v_detint->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_detint->GetYaxis()->SetTitle("Integral (Ph.e.)");
        h_dt_v_detint->SetMarkerStyle(7);
        h_dt_v_detint->Draw();
        c_dt_v_detint->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_v_detint.png");
        c_dt_v_detint->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_detint.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_detint.png"));
        h_dt_v_detint->Reset();

        TCanvas* c_int = new TCanvas("c_int", "Distribution of All Event Integral Values", 1200, 700);
        c_int->SetLogy();
        c_int->cd();
        h_int->GetXaxis()->SetTitle("Integral (Ph.e.)");
        h_int->GetYaxis()->SetTitle("Counts");
        h_int->Draw();
        c_int->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_int.png");
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
        c_int_HL->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_int_HL.png");
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
        c_int_cent_spike->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_int_cent_spike.png");
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
        c_int_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_int_low_dt.png");
        c_int_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_int_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_int_neg_low_dt.png"));
        h_int_pos_low_dt->Reset();
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
        c_svpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_svpint_low_dt.png");
        c_svpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_svpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_svpint_neg_low_dt.png"));
        h_svpint_pos_low_dt->Reset();
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
        c_tvpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_tvpint_low_dt.png");
        c_tvpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_tvpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_tvpint_neg_low_dt.png"));
        h_tvpint_pos_low_dt->Reset();
        h_tvpint_neg_low_dt->Reset();

        TCanvas* c_avpint_neg_low_dt = new TCanvas("c_avpint_neg_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        // c_avpint_neg_low_dt->SetLogy();
        c_avpint_neg_low_dt->cd();
        TF1* func1 = new TF1("func1", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -200, 200);
        func1->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func1->SetParameters(1, 10, 0, 50);
        h_avpint_pos_low_dt->Fit("func1", "R");
        h_avpint_neg_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_neg_low_dt->GetYaxis()->SetTitle("Counts");
        // h_avpint_neg_low_dt->GetYaxis()->SetRange(0, 10);
        h_avpint_neg_low_dt->Draw();
        h_avpint_neg_low_dt->SetLineColor(kRed);
        h_avpint_pos_low_dt->Draw("same");
        h_avpint_pos_low_dt->SetLineColor(kGreen);
        h_avpint_cent_spike->Draw("same");
        h_avpint_cent_spike->SetLineColor(kBlue);
        func1->Draw("same");
        func1->SetLineColor(1);
        TLegend *leg_avpint_neg_low_dt = new TLegend(0.9, 0.65, 0.65, 0.45);
        leg_avpint_neg_low_dt->AddEntry(h_avpint_cent_spike, "Central Spike Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_pos_low_dt, "Positive Low dt Events", "l");
        leg_avpint_neg_low_dt->AddEntry(h_avpint_neg_low_dt, "Negative Low dt Events", "l");
        leg_avpint_neg_low_dt->Draw();
        c_avpint_neg_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_avpint_low_dt.png");
        c_avpint_neg_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_neg_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_neg_low_dt.png"));
        h_avpint_pos_low_dt->Reset();
        h_avpint_neg_low_dt->Reset();

        TCanvas* c_avpint_low_dt = new TCanvas("c_avpint_low_dt", "Distribution of All Low dt Veto Panel Integral Values", 1200, 700);
        // c_avpint_low_dt->SetLogy();
        c_avpint_low_dt->cd();
        TF1* func2 = new TF1("func2", "[0]+[1]*exp(-1*pow((x-[2]),2)/(2*pow([3],2)))", -300, 300);
        func2->SetParNames("Baseline", "Amplitude", "Peak Center", "Standard Deviation");
        func2->SetParameters(1, 10, 0, 50);
        h_avpint_low_dt->Fit("func2", "R");
        h_avpint_low_dt->GetXaxis()->SetTitle("Integral (ADC)");
        h_avpint_low_dt->GetYaxis()->SetTitle("Counts");
        // h_avpint_low_dt->GetYaxis()->SetRange(0, 50);
        h_avpint_low_dt->Draw();
        func2->Draw("same");
        func2->SetLineColor(2);
        c_avpint_low_dt->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_avpint_low_dt_1c.png");
        c_avpint_low_dt->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_avpint_low_dt.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_avpint_low_dt.png"));
        h_avpint_low_dt->Reset();

        TCanvas* c_dt_v_tmuon = new TCanvas("c_dt_v_tmuon", "Event 61 Detector dt vs Muon dt", 1200, 700);
        // c_dt_v_tmuon->SetLogy();
        c_dt_v_tmuon->cd();
        h_dt_v_tmuon->GetXaxis()->SetTitle("Delta-T (ns)");
        h_dt_v_tmuon->GetYaxis()->SetTitle("Time Since Last Muon (ns)");
        h_dt_v_tmuon->SetMarkerStyle(7);
        h_dt_v_tmuon->Draw();
        c_dt_v_tmuon->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots12636/h_dt_v_tmuon.png");
        c_dt_v_tmuon->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_dt_v_tmuon.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_dt_v_tmuon.png"));
        h_dt_v_tmuon->Reset();

        TCanvas* c_allplsvec_size = new TCanvas("c_allplsvec_size", "Number of Events b/t LL and HL", 1200, 700);
        c_allplsvec_size->SetLogy();
        c_allplsvec_size->cd();
        h_allplsvec_size->GetXaxis()->SetTitle("Number of Events");
        h_allplsvec_size->GetYaxis()->SetTitle("Counts");
        h_allplsvec_size->Draw();
        c_allplsvec_size->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_allplsvec_size.png");
        c_allplsvec_size->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_allplsvec_size.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_allplsvec_size.png"));
        h_allplsvec_size->Reset();

        TCanvas* c_inttrigvec_size = new TCanvas("c_inttrigvec_size", "Number of Int. Trig. Events b/t LL and HL", 1200, 700);
        c_inttrigvec_size->SetLogy();
        c_inttrigvec_size->cd();
        h_inttrigvec_size->GetXaxis()->SetTitle("Number of Events");
        h_inttrigvec_size->GetYaxis()->SetTitle("Counts");
        h_inttrigvec_size->Draw();
        c_inttrigvec_size->SaveAs("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/plots7894/h_inttrigvec_size.png");
        c_inttrigvec_size->Close();
        // gPad->Print(Form("/data7/coherent/data/d2o/emward/Detector_Data_Analysis/h_inttrigvec_size.png"));
        // gPad->Print(Form("/mnt/c/Users/eliwa/Documents/Detector_Data_Analysis/plots/h_inttrigvec_size.png"));
        h_inttrigvec_size->Reset();

    }

    std::cout << "\n" << "End of code." << "\n";

    return 0;

}