/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sim;

import cern.jet.random.Poisson;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;
import com.jmatio.io.MatFileReader;
import com.jmatio.types.MLDouble;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import java.util.Arrays;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.Math.signum;
import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import java.util.ArrayList;

/**
 *
 * @author Logan
 */
public class Network {
    /*===Simulation time parameters===*/

    public double t = 0; //ms
    public double dt = 0.1; //ms
    
    /*===Conductance dt (for stretching or shrinking curves)===*/
    public double dtExc = dt;
    public double dtInh = dt;

    /*===Basic network layout parameters===*/
    public int numHCRows, numHC, numExcitatoryPerHC, numInhibitoryPerHC; //HC = "hypercolumn"; model consists of a square array of hypercolumns
    public int numExcitatory, numInhibitory;
    public int numLGNBlockRows, numLGNBlocks; //LGN is organized into a square array of blocks
    public double numLGNPerRowPerBlock, numLGNPerColPerBlock;
    public int numLGNPerBlock, numLGN;
    public int HCWidth, patchWidth;
    public MyPoint[] locExc, locInh, locLGNON, locLGNOFF;

    /*===Connectivity parameters and matrices===*/
    public int[] numEPostExc, numIPostExc, numEPostInh, numIPostInh, numEPreExc, numIPreExc, numEPreInh, numIPreInh; //numEPostExc[i] = number of E cells postsynaptic to the ith E cell; numEPreExc[i] = number of E cells presynaptic to ith E cell
    public int[][] cEExc, cIExc, cEInh, cIInh, cPreEExc, cPreIExc, cPreEInh, cPreIInh; //cIExc[i][j] = index of jth postsynaptic I cell to the E cell at index i; cPreIExc[i][j] = index of jth presynaptic I cell to the E cell at index i
    public double[] randomPicksExc;
    public double[] randomPicksInh;
    public double randomPicksEEConnectivity;

    /*===Connectivity build variables===*/
    public double[] egal;
    public double peakCEI;
    public double peakCIE;
    public double peakCII;
    public double numPreDistCutoff, numPreDistCutoffFB;
    
    public double[] sEgal;

    /*===Connection strengths===*/
    public double sLGN, sLGNInh, sEE, sIE, sEI, sII, sIEOverSEE, sEIOverSEE, sELGNOverSEE, sILGNOverSEE;
    public double sEELowMultiplier, sEEHighMultiplier;
    public double sEEL6Low, sEEL6High, sIEL6;
    public double sIEL6OverSEEL6FracL4;
    public double aChSELGNFac = 1.0;
    public double vertTemplateMod = 1.0;
    public double horizTemplateMod = 1.0;
    public double[] targetMod;
    public double sELGNOverSEERnd;
    /*===Strength variability===*/
    public double sModSTD = 0.1;
    public double[] sModSeedExc, sModSeedInh, sModSeedLGNON, sModSeedLGNOFF;

    /*===Layer 6 recording parameters===*/
    /*===Layer 6 recording parameters===*/
    public boolean recording = false;
    public boolean firstRecordingRound = false;
    public double firstRecordingStart;
    public double firstRecordingEnd;
    public double recordingPeriod = 250.0;
    public double[] fbSpikeTimes;
    public int[] fbSpikes;
    public int numFBSpikes, numFB;
    public int[][] excPostFB, inhPostFB;
    public int[][] fBPreExc, fBPreInh;
    public int[] numFBPreExc, numFBPreInh;
    public int[] numSimplePreExc;

    public double[] noiseFBSpikeTimes;
    public int[] noiseFBSpikes;
    public int numNoiseFBSpikes;
    public int currentNoiseFBSpike;
    public double noiseFBClock;
    public double noiseFBRandEnd;
    public double noiseFBRandStart;
    public double[] noiseFBPComplex;
    
    public double fullContrastLowCpxRate = 0.0;
    public double fullContrastHighCpxRate = 0.0;

    public double sigmaFB;
    public double sigmaNoiseFB;

    public double[] timeLastSpike;
    public double[] timeLastNoiseSpike;

    public int[][] l4L6Map, l6L4Map;
    public int[] numL4L6Map, numL6L4Map;
    public double l4L6Scaling = 0.9;
    public int[][] pastAvgCoupling, presentAvgCoupling;
    public int[] numPastAvgCoupling, numPresentAvgCoupling;
    public double[] pastAvgAngles;
    public ArrayList<ArrayList<Integer>> regionList;
    public double minDistFromCenter;
    public double minDistFromBoundary;
    public int[][] l4L6InstFBMap;
    public int[] numL4L6InstFBMap;


    public double l4L6MapWidth      = 75.0;
    public double instBlockWidth       = 75.0;
    public double numL4InL4L6InnerMap      = 30;
    public double numL4InL4L6OuterMap = 15;
    
    public double numGaussianL4L6Coupling = 100;
    public double fracOriginalL4L6Coupling = 0;
    public double numL4L6GaussianSDs = 1.5;
    public double l4L6GaussianSD = 150.0;
    
    public double totalL4L6Couplings;

    public boolean l6HearsFromSimpleOnly = false;

    public int[] numExcPostFB, numInhPostFB;
    public double[][] locFB;
    public MyPoint[] locFBAsPoints;
    public double fbClock;
    public double fbStart;
    public double fbEnd;
    public double fbRandEnd;
    public double fbRandStart;
    public String fbInputFile;
    public int currentFBSpike = 0;
    public double[] fbPComplex;
     public double[] fbPComplexPast;
    public double sigmaFBPast = 1000.0;
    public double xL4Start = 6.0;
    public double xL4End = 15.0;
    public double minL6 = 0.0;
    public int minL6RefreshPeriod = 100;
    public double minL6MaxL6Ratio = .175;
    public double peakLGNMinL6 = 0.0;
    public double largestMinL6 = 0.0;
    
    public double fbStdSqr = 4.0;    
    public double fbIPeak = .6;
    public double[] pValues;
    
    public double[] l6L4xs = {0.0, 8.0, 13.5, 18.0, 30.0, 1000.0};  
    public double[] l6L4ys = {5.5, 5.5, 20.0, 48.0, 60.0, 60.0 };
    public double[] l6L4slopes = {(5.5-5.5)/(8.0-0.0), (20.0-5.5)/(13.5-8.0), (48.0-20.0)/(18.0-13.5), (60.0-48.0)/(30.0-18.0), 0.0};
    public double[] fbPoissonxs = {0.0, 5.0, 100.0};
    public double[] fbPoissonys = {0.0, 1.33*5.0, 1.33*5.0};
    public double[] fbPoissonSlopes = {1.33, 0.0, 0.0};
    public double fbBGRate = 0.0; // kept for back compat.
    
    public double[] fbWideEgal;
    public double[] fbNarrowEgal;
    public double fbWideNumSDsE;
    public double fbWideSDE;   
    public double fbNarrowNumSDsE;
    public double fbNarrowSDE;
    public double fbWideNumSDsI;
    public double fbWideSDI;   
    public double fbNarrowNumSDsI;
    public double fbNarrowSDI;
    
    public double fbIExpNumPreWide;
    public double fbIExpNumPreNarrow;
    
    public int[][] longRangeCoupling;

    public double[] facxs;
    public double[] facys;
    public double[] facslopes;
    
    public double constFacilitationFac = 0.0;
    
    public double facIOverFacE = 0.0;
    
    /*===FB spike spreading mechanism===*/
    public double[] nextFBSpikeTime;
    public int[] numWaitingFBSpikes;
    public double fbSpreadLow, fbSpreadHigh;
    public int[] fbSpikeSpreadProbChanges; //example: = { 5, 8, 1000 };
    public double[] fbSpikeSpreadProbs; //example: = {1.0, p1, p2};
    
    /*===FB anti-spike mechanism===*/
    public double[] antiSpikeFB;
    public double bgRecRate = 3.0;
    public double fullContrastAntiSpikeFrac = 0.0;
    public double noContrastAntiSpikeFrac = 0.0;
    
    /*===L6 double buffer===*/
    public int[] dBFBSpikes;
    public double[] dBFBSpikeTimes;
    public int numDBFBSpikes;
    public double recordingSmearSD;
    public double l6L4FRRatio;
    public double l6L4FRRatioStart = 1.8;
    public double constL6L4FRRatio;
    //public double l6PoissonSpPerSec;
    
    /*===Default Network parameters===*/
    public double gLExc = 50.0 / 1000.0;
    public double gLInh = gLExc / 0.75;
    public double eRev = 4.67;
    public double iRev = -0.67;
    public double refractoryE = 2.0; //ms
    public double refractoryI = 1.0;
    public double tauAMPAToEa = 1;
    public double tauAMPAToIa = 1;
    public double tauAMPAToEb = 1;
    public double tauAMPAToIb = 1;
    public double tauGABAa = 1.67;
    public double tauGABAb = 1.67;
    public double tauGABAc = 1.67;
    public double tauGABAd = 1.67;
    public double tauAmbToEa, tauAmbToEb, tauAmbToIa, tauAmbToIb;
    public double weightToSecondGABA = 0.0;
    public double tauNMDA1 = 80.0;
    public double tauNMDA2 = 2.0;
    public double eLenScale = 200.0 / sqrt(2);
    public double iLenScale = 125.0 / sqrt(2);

    /*===Neuron state variables===*/
    public double[] vExc, vInh, sleepExc, sleepInh; //E and I membrane potentials, and refractory period clocks
    public double[][] gLGNExc, gLGNInh,
            gAmbExc, gAmbInh,
            gAMPAFBExc, gAMPAFBInh,
            gAMPAExc, gAMPAInh,
            gGABAExc, gGABAInh,
            gNMDAFBExc, gNMDAFBInh,
            gNMDAExc, gNMDAInh; //conductances
    
    /*===Frac NMDA===*/
    public double fracNMDAEE = 10.0; //hopefully to make things obviously wrong if not set
    public double fracNMDAIE = 10.0;
    
    /*===LGN -- Cortex connectivity===*/
    public int[][] cELGNON; //cELGNON[i][j] = index of jth E cell postsynaptic to ith LGN ON cell
    public int[][] cILGNON;
    public int[][] cELGNOFF;
    public int[][] cILGNOFF;
    public int[] nELGNON; //nELGNON[i] = number of postsynaptic E connections from ith LGN ON cell
    public int[] nILGNON;
    public int[] nELGNOFF;
    public int[] nILGNOFF;
    public int[] numInputsExc; //numInputsExc[i] = number of LGN ON and OFF inputs to excitatory cell i
    public int[] numInputsInh;
    public String LGNInputFile; //LGN connectivity data is placed in a MATLAB .mat file, the name of which is stored here

    /*===egalitarian feedback===*/


    /*===feedback simple vs complex mechanism===*/
    public double[] simpleFB;
    public double[] highFiringComplexFB;
    public double[] phaseFB;
    public double fracSimple = 1.0;
    public double simpleAvg = 0.5;
    public double simpleAmpLow = 0.5;
    public double simpleAmpHigh = 0.5;
    public double[] simpleAmp;
    public double deregAmpLow;
    public double deregAmpHigh;
    public double[] deregAmp;
    
    public double simpleSpikeProb = 0.5/1000.0 * 0.1;
    public double highCpxProb = 10.0/1000.0 * 0.1;
    public double lowCpxProb = 0.5/1000.0 * 0.1;

    /*== LGN parameters ==*/
    public double leakLGN = 1.0 / 10.0;
    public double[] phaseON, phaseOFF, vLGNON, vLGNOFF;
    public double spFreqMod;

    /*===LGN Drift grating parameters===*/
    public double I_0 = 100.0 / 1000.0;
    public double shotNoise = 0.08;
    public double noiseRate = 1.0;
    public double threshLGN = 1.15;
    public double contrast = 2.0;
    public double tf = 4.0 / 1000.0;
    public double spFreq = .025;
    public MyPoint orientation = new MyPoint(1, 0); //MyPoint is just an extension of the point2D class
    public double angle = 0.0;

    /*===Extra Drive parameters (NOTE: extra drive conductance will be a part of LGN conductance)===*/
    public double rateEAmb;
    public double rateIAmb;
    public double sAmb = 0.01;

    /*===Quick stats collecting mechanism===*/
    public int statSkipInSteps = 20;
    public int statStep = 0; //will be incremented mod statStkipsInSteps every step
    public double patchOneFRExc = 0;
    public double patchOneFRExcTemp = 0;
    public double patchOneFRInh = 0;
    public double patchOneFRInhTemp = 0;
    public double patchTwoFRExc = 0;
    public double patchTwoFRExcTemp = 0;
    public double patchTwoFRInh = 0;
    public double patchTwoFRInhTemp = 0;
    public int[] membershipExc;
    public int[] membershipInh;

    
    /*===SIE depression mechanism===*/
    public double[] sEIDepLookback;
    public double sEIDepDelta;
    public double[] sEIDepPerSpike = {0, 1, 1.33};
    public double[][] last3Spikes;
    public int[] last3SpikesInd;
    public int[] sEIDepStageNowInh;
    public double sEIDepAmt = .1;
    public double sEIDepTimeConstant = 20.0;
    
    /*===E Fatigue mechanism===*/
    public double eFatSigma;
    public double eFatAmt;
    
    /*===Time Averaging===*/
    public int numAvgs = 40;
    public double[][] timeAvgs;
    public double timeAvgDt = 0.1; //units are seconds
    public int whichAvg = 0;
    public double nextSwitchTime = timeAvgDt / (double)numAvgs;
    public double[] timeAvgWeights;    
    
    /*===Random number generators===*/
    public long seedL4 = 1, seedL6 = 2, seedMisc = 3, seedRun = 4;
    public int seedEngine = 5;
    public RandomEngine engine = new DRand((int) seedEngine);//new DRand(new java.util.Date());
    public Poisson poisson = new Poisson(1.0, engine);
    public Normal normal = new Normal(0.0, 1.0, engine);

    public java.util.Random rndGenL4 = new java.util.Random(  (long) seedL4);
    public java.util.Random rndGenL6 = new java.util.Random(  (long) seedL6);;
    public java.util.Random rndGenMisc = new java.util.Random((long) seedMisc);;
    public java.util.Random rndGenRun = new java.util.Random( (long) seedRun);;

    /*===Stats===*/
    public int eSpikes = 0;
    public int iSpikes = 0;
    public int numSteps = 0;
    public int[] spikesNowExc;
    public int[] spikesNowInh;
    public double[] spikeTimesNowExc;
    public double[] spikeTimesNowInh;
    public int[] spikesNowLGNON;
    public int[] spikesNowLGNOFF;
    public double[] spikeTimesNowLGNON;
    public double[] spikeTimesNowLGNOFF;
    public int numSpikesNowExc;
    public int numSpikesNowInh;
    public int numSpikesNowLGNON;
    public int numSpikesNowLGNOFF;    
    public int[] l6SpikesNow;
    public double[] l6SpikeTimesNow;
    public int numL6SpikesNow;
    
    /*firing rate mechanism*/
    public double[] preRateExc;
    public double[] preRateInh;
    public double[] preRateFB;
    public double rateStartTime;
    /*current record mechanism*/
    public double[] avgVExc;
    public double[] avgVInh;
    public double[][] currentToE;
    public double[][] currentToI;
    public double stepsBetweenCurrentSamples=50;
    public double stepsToNextCurrentSample=50;
    public double numCurrentSamples=0;

    /*===I randomized threshold===*/
    public double[] threshEBounds;// length 2 array: {lower bound, upper bound}, construction prm
    public double[] threshIBounds;// length 2 array: {lower bound, upper bound}
    public double[] threshI;
    public double[] threshE;
    public double[] fatigueExc;
    public double eFatCap = 1000.0;

    /*===I input modifier===*/
    public double sIIOverSEILow = .9;
    public double sIIOverSEIHigh = 1.1;
    public double sIIOverSEIRange;

    /*===Fast spike initial conductance for strength 1===*/
    public double[] fastSpikeInitial = {0.041193481496372, 0.037000133080574, 0.033233652467581, 0.029850586048726};//1.5sec speedup////{0.055017482305651,0.032944600183025,0.019727305498817,0.011812757783723};//1 sec speedup//{0.047552705667974,0.035593342565849,0.026641723477432,0.019941409788497};//1.25 sec speedup////nospeedup:{ 0.093999733184618, 0.005628726538001, 0.000337049493293, 0.000020182604389};//{0.041177,0.036986,0.033221,0.029839}; //exp,lin-exp,quad-exp,cub-exp
    public double tauInitialGABADamper, tauInitialAMPADamper;
    public double[] initialGABADamperExc, initialAMPADamperExc, initialAMPADamperExcFB;
    public double[] initialGABADamperInh, initialAMPADamperInh, initialAMPADamperInhFB;
    
    public double eOffset = dt;
    
    /*===E->E spike delay cylinder (a data structure for quickly implementing E hits with randomized delays)===*/
    //public double maxSpikeDelayMs = 6.0;
    //public int maxSpikeDelaySteps = (int) ceil(maxSpikeDelayMs / dt);
//    public double spikeDelayMs = 3.0;
    public int spikeDelaySteps = 30;//(int) ceil(spikeDelayMs / dt);
    public double[][] delayCylinder;
    public double[][] delayCylinderFB;
    public int spikeDelayIndex = 0;
    public int lgnSpikeDelaySteps;
    public int minSpikeDelaySteps = 0;
    public int spikeDelayActualSteps = 30;

    /*=== PATH variables ===*/
    public String winPATH = "C:\\Users\\Logan\\mfiles\\Sim\\";
    public String linPATH = "/home/charikar/mfiles/Sim/";
    public String macPATH = "/arc/1.3/p3/logan/mfiles/Sim/";

    /*=== Build choices ===*/
    public boolean load4CConnections;
    public String l4CConnectionsFile;
    public boolean loadFBConnections;
    public String fBConnectionsFile;
    public boolean loadMiscRandComponents;
    public String miscRandComponentsFile;

    /*=== for degugging ===*/
    public int debugCount = 0;
    public double highFiringComplexPerc;
    public int MAXPASTAVGCOUPLING = 150;
    public double cornerBoost;
    public double cornerWidth;
    public double lGNSpFreqRightStretch = 1.3;
    
    /*===Empty Constructor===*/
    public Network() {
    }

    /*set 4C layout, allocate state and initialize*/
    public void layout4C() {
        numHC = numHCRows * numHCRows;
        numExcitatory = numExcitatoryPerHC * numHC;
        numInhibitory = numInhibitoryPerHC * numHC;

        HCWidth = 500; //in micrometers
        patchWidth = HCWidth * numHCRows;

        /*===Network state===*/
        vExc = new double[numExcitatory];
        vInh = new double[numInhibitory];
        sleepExc = new double[numExcitatory];
        sleepInh = new double[numInhibitory];

        /*===Network state helpers===*/
        gLGNExc = new double[numExcitatory][3];
        gLGNInh = new double[numInhibitory][3];
        gAmbExc = new double[numExcitatory][3];
        gAmbInh = new double[numInhibitory][3];
        gAMPAFBExc = new double[numExcitatory][3];
        gAMPAFBInh = new double[numInhibitory][3];
        gAMPAExc = new double[numExcitatory][3];
        gAMPAInh = new double[numInhibitory][3];
        gGABAExc = new double[numExcitatory][5];
        gGABAInh = new double[numInhibitory][5];
        gNMDAFBExc = new double[numExcitatory][3];
        gNMDAFBInh = new double[numInhibitory][3];
        gNMDAExc = new double[numExcitatory][3];
        gNMDAInh = new double[numInhibitory][3];
        
        sEIDepStageNowInh = new int[numInhibitory];
        last3Spikes = new double[numInhibitory][3];
        last3SpikesInd = new int[numInhibitory];
        for(int j=0;j<3;j++){
            for(int i=0;i<numInhibitory;i++){
                last3Spikes[i][j]=-1000.0;
            }
        }

        /*=== Place neurons ===*/
        int ePopWidth = (int) ceil(sqrt(numExcitatory)); //in neurons
        int iPopWidth = (int) ceil(sqrt(numInhibitory)); //in neurons
        locExc = new MyPoint[numExcitatory];
        locInh = new MyPoint[numInhibitory];
        for (int i = 0; i < numExcitatory; i++) {
            locExc[i] = new MyPoint((double) (i % ePopWidth) * (patchWidth / (double) ePopWidth),
                    (double) (i / ePopWidth) * (patchWidth / (double) ePopWidth));
        }
        for (int i = 0; i < numInhibitory; i++) {
            locInh[i] = new MyPoint((double) (i % iPopWidth) * (patchWidth / (double) iPopWidth),
                    (double) (i / iPopWidth) * (patchWidth / (double) iPopWidth));
        }
        /*=== Neurons placed ===*/

        /*=== allocate and initialize other neuron state ===*/
        membershipExc = new int[numExcitatory];
        membershipInh = new int[numInhibitory];

        spikesNowExc = new int[numExcitatory];
        spikesNowInh = new int[numInhibitory];
        spikeTimesNowExc = new double[numExcitatory];
        spikeTimesNowInh = new double[numInhibitory];
        
        preRateExc = new double[numExcitatory];
        preRateInh = new double[numInhibitory];
        
        avgVExc = new double[numExcitatory];
        avgVInh = new double[numInhibitory];
        currentToE = new double[numExcitatory][6];//e fb lgn amb i leak
        currentToI = new double[numExcitatory][6];

        threshI = new double[numInhibitory];
        threshE = new double[numExcitatory];
        fatigueExc = new double[numExcitatory];

        sModSeedExc = new double[numExcitatory];
        sModSeedInh = new double[numInhibitory];
        
        last3Spikes = new double[numInhibitory][3];
        last3SpikesInd = new int[numInhibitory];      
        for(int j=0;j<3;j++){
            for(int i=0;i<numInhibitory;i++){
                last3Spikes[i][j]=-1000.0;
            }
        }
        

        delayCylinder = new double[numExcitatory][spikeDelaySteps];
        delayCylinderFB = new double[numExcitatory][spikeDelaySteps];
        
        initialGABADamperExc   = new double[numExcitatory];
        initialAMPADamperExc   = new double[numExcitatory];
        initialAMPADamperExcFB = new double[numExcitatory];
        initialGABADamperInh   = new double[numInhibitory];  
        initialAMPADamperInh   = new double[numInhibitory];  
        initialAMPADamperInhFB = new double[numInhibitory];  

        numInputsExc = new int[numExcitatory];
        numInputsInh = new int[numInhibitory];

        randomPicksExc = new double[numExcitatory];
        randomPicksInh = new double[numInhibitory];
    }

    /* allocate and initialize some LGN state variables */
    public void layoutLGN() {
        numLGNBlocks = numLGNBlockRows * numLGNBlockRows;
        this.numLGNPerBlock = (int) Math.round(numLGNPerRowPerBlock * numLGNPerColPerBlock);
        numLGN = numLGNBlocks * numLGNPerBlock;

        /*==LGN state==*/
        vLGNON = new double[numLGN];
        vLGNOFF = new double[numLGN];
        phaseON = new double[numLGN];
        phaseOFF = new double[numLGN];

        sModSeedLGNON = new double[numLGN];
        sModSeedLGNOFF = new double[numLGN];

        nELGNON = new int[numLGN]; //read: "number connections to E from LGN ON"
        nILGNON = new int[numLGN]; //read: "number connections to E from LGN ON"
        nELGNOFF = new int[numLGN]; //read: "number connections to E from LGN ON"
        nILGNOFF = new int[numLGN]; //read: "number connections to E from LGN ON"

        locLGNON = new MyPoint[numLGN];
        locLGNOFF = new MyPoint[numLGN];

        /*===stats===*/
        spikesNowLGNON  = new int[numLGN];
        spikesNowLGNOFF = new int[numLGN];
        spikeTimesNowLGNON  = new double[numLGN];
        spikeTimesNowLGNOFF = new double[numLGN];
    }

    /*deals just with layout of 4Cand LGN*/
    public void layoutAll() {
        layout4C();
        layoutLGN();
    }

    public void connectLGN() {
        /*== Get LGN inputs ==*/
        try {
            System.out.println("Loading LGN file "+LGNInputFile+"...");
            MatFileReader matfilereader;
            if (System.getProperty("os.name").equals("Linux")) {
                matfilereader = new MatFileReader(linPATH.concat(LGNInputFile));
            } else if (System.getProperty("os.name").equals("Mac OS X")){
                matfilereader = new MatFileReader(macPATH.concat(LGNInputFile));
            } else {
                matfilereader = new MatFileReader(winPATH.concat(LGNInputFile));
            }

            /*==Import out-degrees and out-connectivities
             *of LGN cells (these will need further processing,
             *like changing to java indexing starting at 0)==*/
            double[][] nELGNONPre = ((MLDouble) matfilereader.getMLArray("numLGNONToECortex")).getArray();
            double[][] nELGNOFFPre = ((MLDouble) matfilereader.getMLArray("numLGNOFFToECortex")).getArray();
            double[][] nILGNONPre = ((MLDouble) matfilereader.getMLArray("numLGNONToICortex")).getArray();
            double[][] nILGNOFFPre = ((MLDouble) matfilereader.getMLArray("numLGNOFFToICortex")).getArray();
            double[][] cELGNONPre = ((MLDouble) matfilereader.getMLArray("LGNONToECortex")).getArray();
            double[][] cELGNOFFPre = ((MLDouble) matfilereader.getMLArray("LGNOFFToECortex")).getArray();
            double[][] cILGNONPre = ((MLDouble) matfilereader.getMLArray("LGNONToICortex")).getArray();
            double[][] cILGNOFFPre = ((MLDouble) matfilereader.getMLArray("LGNOFFToICortex")).getArray();
            double[][] maxNumLGNToCortex = ((MLDouble) matfilereader.getMLArray("maxNumLGNToCortex")).getArray();
            double[][] corticalLocLGNON = ((MLDouble) matfilereader.getMLArray("corticalLocLGNON")).getArray();
            double[][] corticalLocLGNOFF = ((MLDouble) matfilereader.getMLArray("corticalLocLGNOFF")).getArray();
            double[][] numInputsExcPre = ((MLDouble) matfilereader.getMLArray("numInputsExc")).getArray();
            double[][] numInputsInhPre = ((MLDouble) matfilereader.getMLArray("numInputsInh")).getArray();
            double[][] randomPicksExcPre = ((MLDouble) matfilereader.getMLArray("randomPicksExc")).getArray();

            /*==Imported, now process them==*/
            int maxOutDegree = (int) maxNumLGNToCortex[0][0];
            cELGNON = new int[numLGN][maxOutDegree];
            cELGNOFF = new int[numLGN][maxOutDegree];
            cILGNON = new int[numLGN][maxOutDegree];
            cILGNOFF = new int[numLGN][maxOutDegree];
            nELGNON = new int[numLGN];
            nELGNOFF = new int[numLGN];
            nILGNON = new int[numLGN];
            nILGNOFF = new int[numLGN];

            for (int i = 0; i < numLGN; i++) {
                nELGNON[i] = (int) nELGNONPre[i][0];
                nELGNOFF[i] = (int) nELGNOFFPre[i][0];
                nILGNON[i] = (int) nILGNONPre[i][0];
                nILGNOFF[i] = (int) nILGNOFFPre[i][0];
                for (int j = 0; j < nELGNON[i]; j++) {
                    cELGNON[i][j] = (int) cELGNONPre[i][j] - 1;
                }
                for (int j = 0; j < nELGNOFF[i]; j++) {
                    cELGNOFF[i][j] = (int) cELGNOFFPre[i][j] - 1;
                }
                for (int j = 0; j < nILGNON[i]; j++) {
                    cILGNON[i][j] = (int) cILGNONPre[i][j] - 1;
                }
                for (int j = 0; j < nILGNOFF[i]; j++) {
                    cILGNOFF[i][j] = (int) cILGNOFFPre[i][j] - 1;
                }
                locLGNON[i] = new MyPoint(corticalLocLGNON[i][0], corticalLocLGNON[i][1]);
                locLGNOFF[i] = new MyPoint(corticalLocLGNOFF[i][0], corticalLocLGNOFF[i][1]);
            }
            for (int i = 0; i < numExcitatory; i++) {
                numInputsExc[i] = (int) numInputsExcPre[i][0];
            }
            for (int i = 0; i < numInhibitory; i++) {
                numInputsInh[i] = (int) numInputsInhPre[i][0];
            }
            for (int i = 0; i < numExcitatory;i++) {
                randomPicksExc[i] = randomPicksExcPre[i][0];
            }
            /*==Processing LGN->cortex connection complete==*/
        } catch (IOException ex) {
            System.out.println("Error in LGN import code");
            System.exit(1);
        }
        /*== Got LGN inputs ==*/
    }

    public void connect4C() {
        if (load4CConnections) {
            System.out.println("loading 4C connections from " + l4CConnectionsFile + "...");
            loadConnectivity(l4CConnectionsFile);
        } else {
            System.out.println("generating new 4C connections network...");
            /*== setup cortical network ==*/
            allocateConnectivity();
            setupCortexNetwork();
            /*== cortical network setup ==*/
        }
    }

    /* allocate space to place the connectivity matrix (and do some other initializations) */
    public void allocateConnectivity() {

        /*===Calculations to improve efficient memory allocation for connectivity matrix===*/
        double EDensity = (double) numExcitatoryPerHC / (double) (HCWidth * HCWidth);
        double IDensity = (double) numInhibitoryPerHC / (double) (HCWidth * HCWidth);
        double EDiscArea = PI * pow(2 * eLenScale, 2);
        double IDiscArea = PI * pow(2 * iLenScale, 2);
        double numEInEDisc = EDiscArea * EDensity;
        double numIInEDisc = EDiscArea * IDensity;
        double numEInIDisc = IDiscArea * EDensity;
        double numIInIDisc = IDiscArea * IDensity;
        int maxEToEOutDegree = (int) floor(egal[0] * numEInEDisc);
        int maxEToIOutDegree = (int) floor(peakCIE * numIInEDisc);
        int maxIToEOutDegree = (int) floor(peakCEI * numEInIDisc);
        int maxIToIOutDegree = (int) floor(peakCII * numIInIDisc);
        int maxEToEInDegree = maxEToEOutDegree;
        int maxEToIInDegree = (int) floor(peakCIE * numEInEDisc);
        int maxIToEInDegree = (int) floor(peakCEI * numIInIDisc);
        int maxIToIInDegree = maxIToIOutDegree;

        System.out.println("Showing allocated in/out degrees (look in Network.java for meaning):");
        System.out.println(maxEToEOutDegree);
        System.out.println(maxEToIOutDegree);
        System.out.println(maxIToEOutDegree);
        System.out.println(maxIToIOutDegree);
        System.out.println(maxEToEInDegree);
        System.out.println(maxEToIInDegree);
        System.out.println(maxIToEInDegree);
        System.out.println(maxIToIInDegree);

        /*===Network Connectivity===*/
        numEPostExc = new int[numExcitatory]; //num E cells postsynaptic to excitatory [array element]
        numIPostExc = new int[numExcitatory];
        numEPostInh = new int[numInhibitory];
        numIPostInh = new int[numInhibitory];
        numEPreExc = new int[numExcitatory];
        numIPreExc = new int[numExcitatory];
        numEPreInh = new int[numInhibitory];
        numIPreInh = new int[numInhibitory];
        cEExc = new int[numExcitatory][maxEToEOutDegree];
        cIExc = new int[numExcitatory][maxEToIOutDegree];
        cEInh = new int[numInhibitory][maxIToEOutDegree];
        cIInh = new int[numInhibitory][maxIToIOutDegree];
        cPreEExc = new int[numExcitatory][maxEToEInDegree];
        cPreIExc = new int[numExcitatory][maxIToEInDegree];
        cPreEInh = new int[numInhibitory][maxEToIInDegree];
        cPreIInh = new int[numInhibitory][maxIToIInDegree];

    }

    public void setMiscRandComponents() {
        if (loadMiscRandComponents) {
            System.out.println("loading misc randomized components from "+miscRandComponentsFile+"...");
            loadMiscRandComponents(miscRandComponentsFile);
        } else {
            System.out.println("generating new set of all misc randomized components");
            /*== Initialize I,E thresholds ==*/
            setThreshERange(threshEBounds[0], threshEBounds[1]);
            setThreshIRange(threshIBounds[0], threshIBounds[1]);

            /*== set strength modifiers ==*/
            setSModSeeds();


        }
    }

    public void setupCortexNetwork() {
        /*=== make sure degree numbers are set to 0 ===*/
        for (int i = 0; i < numExcitatory; i++) {
            numEPreExc[i] = 0;
            numEPostExc[i] = 0;
            numIPreExc[i] = 0;
            numIPostExc[i] = 0;
        }
        for (int i = 0; i < numInhibitory; i++) {
            numEPreInh[i] = 0;
            numEPostInh[i] = 0;
            numIPreInh[i] = 0;
            numIPostInh[i] = 0;
        }
        /*=== degrees reset ===*/

        /* calculate expected number presynapic for each pair type*/
        double EDensity = (double) numExcitatoryPerHC / ((double) HCWidth * (double) HCWidth); // num E per square micrometer
        double IDensity = (double) numInhibitoryPerHC / ((double) HCWidth * (double) HCWidth); // num E per square micrometer

        double EDensityPerESD = EDensity * eLenScale * eLenScale;
        double EDensityPerISD = EDensity * iLenScale * iLenScale;
        double IDensityPerESD = IDensity * eLenScale * eLenScale;
        double IDensityPerISD = IDensity * iLenScale * iLenScale;

        double desiredNumPrePeak1EToE = 2.0 * PI * (1.0 - exp(-2.0)) * EDensityPerESD;
        double desiredNumPrePeak1IToE = 2.0 * PI * (1.0 - exp(-2.0)) * IDensityPerISD;
        double desiredNumPrePeak1EToI = 2.0 * PI * (1.0 - exp(-2.0)) * EDensityPerESD;
        double desiredNumPrePeak1IToI = 2.0 * PI * (1.0 - exp(-2.0)) * IDensityPerISD;

        System.out.println("Expected # presynaptic for connectivities (EE,IE,EI,II)=(.15,.6,.6,.6):");
        System.out.println("E->E: " + (desiredNumPrePeak1EToE * .15));
        System.out.println("I->E: " + (desiredNumPrePeak1IToE * .6));
        System.out.println("E->I: " + (desiredNumPrePeak1EToI * .6));
        System.out.println("I->I: " + (desiredNumPrePeak1IToI * .6));

        /* make peakCEI, peakCIE, peakCII constant arrays, to be compatible with choosePresynaptic function */
        double[] EIegal = {peakCEI, peakCEI, peakCEI, peakCEI, peakCEI, peakCEI, peakCEI};
        double[] IEegal = {peakCIE, peakCIE, peakCIE, peakCIE, peakCIE, peakCIE, peakCIE, peakCIE, peakCIE};
        double[] IIegal = {peakCII, peakCII, peakCII, peakCII, peakCII, peakCII, peakCII, peakCII, peakCII};
        
        /*== E->E and I->E ==*/
        for (int i = 0; i < numExcitatory; i++) {
            //choosePresynaptic(postIndex,postLocation,sizeOfPreGroup,preGroupLocs, lenScale, expNumPre, cOut, numCOut, cIn, numCIn)
            /*== choose E presynaptic ==*/
            /*== first compute the required bump to connectivity to keep expected num presynaptic the same ==*/
            if (i % 3000 == 0) {
                System.out.println("Connecting E neuron " + i);
            }
            choosePresynaptic(i,
                    locExc[i],
                    numExcitatory,
                    locExc,
                    eLenScale,
                    numInputsExc,
                    cEExc,
                    numEPostExc,
                    cPreEExc,
                    numEPreExc,
                    egal,
                    desiredNumPrePeak1EToE,
                    randomPicksExc);

            /*== choose I presynaptic ==*/
            choosePresynaptic(i,
                    locExc[i],
                    numInhibitory,
                    locInh,
                    iLenScale,
                    numInputsExc,
                    cEInh,
                    numEPostInh,
                    cPreIExc,
                    numIPreExc,
                    EIegal,
                    desiredNumPrePeak1IToE,
                    randomPicksExc);
        }

        /*== redraw outliers ==*/
        System.out.println("redrawing E->E outliers...");
        //go through E populations for E->E
        for (int numLGNPre = 0; numLGNPre < 7; numLGNPre++) {
            System.out.println("numLGN=" + numLGNPre);
            //find all E cells with numLGNPre LGN inputs
            //make sure not to pick any random LGN template cells!
            int[] picks = new int[numExcitatory];
            int numPicks = 0;
            for (int i = 0; i < numExcitatory; i++) {
                if (numInputsExc[i] == numLGNPre){// && randomPicksExc[i] < 0.5 ) {
                    picks[numPicks] = i;
                    numPicks++;
                }
            }
            //found all E cells with numLGNPre LGN inputs

            //if none, go on to next numLGNPre
            if (numPicks==0) continue;

            //find average number presynaptic
            double mean = 0;
            for (int i=0;i<numPicks;i++){
                mean += (double) numEPreExc[picks[i]];
            }
            mean /= (double) numPicks;
            //found avg num presynaptic

            //get lower and upper bound for num presynaptic
            double lowerbd      = mean*(1 - numPreDistCutoff);
            double upperbd      = mean*(1 + numPreDistCutoff);
            System.out.println("lowerbd="+lowerbd);
            System.out.println("upperbd="+upperbd);
            //got

            //list all out of bounds of picks
            int[] outOfBounds = new int[numPicks];
            int numOutOfBounds= 0;
            for (int i=0;i<numPicks;i++){
                int numPre = numEPreExc[picks[i]];
                if (numPre < lowerbd || numPre > upperbd){                    
                    outOfBounds[numOutOfBounds] = i;
                    numOutOfBounds++;
                }
            }
            //listed

            //repeatedly redraw until no more out of bounds
            int redrawCount = 0;
            int MAXREDRAW = 15;
            while(numOutOfBounds > 0 && redrawCount < MAXREDRAW){
                redrawCount++;
                System.out.println("E->E, "+numLGNPre+" LGN: redrawing "+numOutOfBounds+"...");
                for (int i=0;i<numOutOfBounds;i++){
                    int pick = picks[outOfBounds[i]];
                    deletePresynaptic(pick,
                            cEExc,
                            numEPostExc,
                            cPreEExc,
                            numEPreExc);
                    choosePresynaptic(pick,
                            locExc[pick],
                            numExcitatory,
                            locExc,
                            eLenScale,
                            numInputsExc,
                            cEExc,
                            numEPostExc,
                            cPreEExc,
                            numEPreExc,
                            egal,
                            desiredNumPrePeak1EToE,
                            randomPicksExc);
                }

                //mark all out of bounds
                numOutOfBounds = 0;
                for (int i = 0; i < numPicks; i++) {
                    int numPre = numEPreExc[picks[i]];
                    if (numPre < lowerbd || numPre > upperbd){
                        outOfBounds[numOutOfBounds] = i;
                        numOutOfBounds++;
                    }
                }
                //marked
            }
            //redrew E->E
        }

        //same for I->E
        System.out.println("redrawing I->E outliers...");
        //find average number presynaptic
        double mean = 0;
        for (int i = 0; i < numExcitatory; i++) {
            mean += (double) numIPreExc[i];
        }
        mean /= (double) numExcitatory;
        //found avg num presynaptic

        //get lower and upper bound for num presynaptic
        double lowerbd = mean * (1 - numPreDistCutoff);
        double upperbd = mean * (1 + numPreDistCutoff);
        //got

        //mark all out of bounds
        int[] outOfBounds = new int[numExcitatory];
        int numOutOfBounds = 0;
        for (int i = 0; i < numExcitatory; i++) {
            int numPre = numIPreExc[i];
            if (numPre < lowerbd || numPre > upperbd) {
                outOfBounds[numOutOfBounds] = i;
                numOutOfBounds++;
            }
        }
        //marked

        //repeatedly redraw until no more out of bounds
        int redrawCount = 0;
        int MAXREDRAW = 15;
        while (numOutOfBounds > 0 && redrawCount < MAXREDRAW) {
            redrawCount++;
            System.out.println("I->E: redrawing " + numOutOfBounds + "...");
            for (int i = 0; i < numOutOfBounds; i++) {
                int pick = outOfBounds[i];
                deletePresynaptic(pick,
                        cEInh,
                        numEPostInh,
                        cPreIExc,
                        numIPreExc);
                choosePresynaptic(pick,
                        locExc[pick],
                        numInhibitory,
                        locInh,
                        iLenScale,
                        numInputsExc,
                        cEInh,
                        numEPostInh,
                        cPreIExc,
                        numIPreExc,
                        EIegal,
                        desiredNumPrePeak1IToE,
                        randomPicksExc);
            }

            //mark all out of bounds
            numOutOfBounds = 0;
            for (int i = 0; i < numExcitatory; i++) {
                int numPre = numIPreExc[i];
                if (numPre < lowerbd || numPre > upperbd) {
                    outOfBounds[numOutOfBounds] = i;
                    numOutOfBounds++;
                }
            }
            //marked
        }
        //redrew I->E

        System.out.println("done.");

        /*== E->I and I->I ==*/
        for (int i = 0; i < numInhibitory; i++) {
            if (i % 3000 == 0) {
                System.out.println("Connecting I neuron " + i);
            }
            /*== choose E presynaptic ==*/
            choosePresynaptic(i,
                    locInh[i],
                    numExcitatory,
                    locExc,
                    eLenScale,
                    numInputsInh,
                    cIExc,
                    numIPostExc,
                    cPreEInh,
                    numEPreInh,
                    IEegal,
                    desiredNumPrePeak1EToI,
                    randomPicksInh);

            /*== choose I presynaptic ==*/
            choosePresynaptic(i,
                    locInh[i],
                    numInhibitory,
                    locInh,
                    iLenScale,
                    numInputsInh,
                    cIInh,
                    numIPostInh,
                    cPreIInh,
                    numIPreInh,
                    IIegal,
                    desiredNumPrePeak1IToI,
                    randomPicksInh);
        }
        /*== Chose all cortical connections ==*/
    }

    public final double quickNormInv(double alpha) {
        return (10 / log(41)) * log(1 - log(-log(alpha) / log(2)) / log(22));
    }

    public final double modifier(double seedi, double seedj) {
        return 1.0 + sModSTD * quickNormInv((seedi * 1024.0 + seedj) % 1.0);
    }
    
    public void registerESpike(int index){
        spikesNowExc[numSpikesNowExc]=index;
        spikeTimesNowExc[numSpikesNowExc]=t;
        numSpikesNowExc++;
        preRateExc[index]++;
    }
    
    public void registerISpike(int index){
        spikesNowInh[numSpikesNowInh]=index;
        spikeTimesNowInh[numSpikesNowInh]=t;
        numSpikesNowInh++;
        preRateInh[index]++;
    }
    
    public void registerFBSpike(int index){
        l6SpikesNow[numL6SpikesNow]=index;
        l6SpikeTimesNow[numL6SpikesNow]=t;
        numL6SpikesNow++;
        preRateFB[index]++;
    }

    public void evenOutL6ToE(int[] numToAdd, int[] numToDel, double SD, double numSDs){
        addPreToAll(excPostFB,
                    numExcPostFB,
                    fBPreExc,
                    numFBPreExc,
                    locExc,
                    locFBAsPoints,
                    SD,
                    numSDs,
                    numToAdd,
                    "L6->E");
        delPreFromAll(excPostFB,
                      numExcPostFB,
                      fBPreExc,
                      numFBPreExc,
                      numToDel,
                      "L6->E");
    }
    public void evenOutL6ToI(int[] numToAdd, int[] numToDel, double SD, double numSDs){
        addPreToAll(inhPostFB,
                    numInhPostFB,
                    fBPreInh,
                    numFBPreInh,
                    locInh,
                    locFBAsPoints,
                    SD,
                    numSDs,
                    numToAdd,
                    "L6->I");
        delPreFromAll(inhPostFB,
                      numInhPostFB,
                      fBPreInh,
                      numFBPreInh,
                      numToDel,
                      "L6->I");
    }
    public void evenOutEToE(int[] numToAdd, int[] numToDel, double SD, double numSDs){
        String connString="E->E";
        addPreToAll(cEExc,
                    numEPostExc,
                    cPreEExc,
                    numEPreExc,
                    locExc,
                    locExc,
                    SD,
                    numSDs,
                    numToAdd,
                    connString);
        delPreFromAll(cEExc,
                      numEPostExc,
                      cPreEExc,
                      numEPreExc,
                      numToDel,
                      connString);
    }
    public void evenOutEToI(int[] numToAdd, int[] numToDel, double SD, double numSDs){
        String connString="E->I";
        addPreToAll(cIExc,
                    numIPostExc,
                    cPreEInh,
                    numEPreInh,
                    locInh,
                    locExc,
                    SD,
                    numSDs,
                    numToAdd,
                    connString);
        delPreFromAll(cIExc,
                      numIPostExc,
                      cPreEInh,
                      numEPreInh,
                      numToDel,
                      connString);
    }
    public void evenOutIToE(int[] numToAdd, int[] numToDel, double SD, double numSDs){
        String connString="I->E";
        addPreToAll(cEInh,
                    numEPostInh,
                    cPreIExc,
                    numIPreExc,
                    locExc,
                    locInh,
                    SD,
                    numSDs,
                    numToAdd,
                    connString);
        delPreFromAll(cEInh,
                      numEPostInh,
                      cPreIExc,
                      numIPreExc,
                      numToDel,
                      connString);
    }
    public void evenOutIToI(int[] numToAdd, int[] numToDel, double SD, double numSDs){
        String connString="I->I";
        addPreToAll(cIInh,
                    numIPostInh,
                    cPreIInh,
                    numIPreInh,
                    locInh,
                    locInh,
                    SD,
                    numSDs,
                    numToAdd,
                    connString);
        delPreFromAll(cIInh,
                      numIPostInh,
                      cPreIInh,
                      numIPreInh,
                      numToDel,
                      connString);
    }
    
    public void addPreToAll(int[][] postFromPre,
                            int[] numPostFromPre,
                            int[][] preToPost,
                            int[] numPreToPost,
                            MyPoint[] locPost,
                            MyPoint[] locPre,
                            double SD,
                            double numSDs,
                            int[] numToAdd,
                            String connectionType){        

        for (int i=0;i<numToAdd.length;i++){
            if (i%1000==0){
                System.out.println("Adding to connection type " + connectionType + ", cell " + i);
            }
            if (numToAdd[i]==0) continue;
            
            double numSDsActual;
            if (  (locPost[i].x < cornerWidth & locPost[i].y < cornerWidth) 
                | (locPost[i].x < cornerWidth & locPost[i].y > 1500.0-cornerWidth) 
                | (locPost[i].x > 1500.0-cornerWidth & locPost[i].y < cornerWidth) 
                | (locPost[i].x > 1500.0-cornerWidth & locPost[i].y > 1500.0-cornerWidth) ){
                numSDsActual = numSDs*cornerBoost;
            } else{
                numSDsActual = numSDs;
            }
            
            addNPresynaptic(i,
                            postFromPre,
                            numPostFromPre,
                            preToPost,
                            numPreToPost,
                            locPost[i],
                            locPre,
                            SD,
                            numSDsActual,
                            numToAdd[i]);
        }
    }
   
    public void delPreFromAll(int[][] postFromPre,
                              int[] numPostFromPre,
                              int[][] preToPost,
                              int[] numPreToPost,
                              int[] numToDel,
                              String connectionType){

        for (int i=0;i<numToDel.length;i++){
            if (i%1000==0){
                System.out.println("Deleting from connection type " + connectionType + ", cell " + i);
            }
            if (numToDel[i]==0) continue;
            
            deleteNPresynaptic(i,
                               postFromPre,
                               numPostFromPre,
                               preToPost,
                               numPreToPost,
                               numToDel[i]);
        }
    }
    
    public void locFBToPoints(){
        locFBAsPoints = new MyPoint[numFB];
        
        for (int i=0;i<numFB;i++){
            locFBAsPoints[i] = new MyPoint(locFB[i][0],locFB[i][1]);
        }
    }
    
    public void addNPresynaptic(int post,
                                int[][] postFromPre,
                                int[] numPostFromPre,
                                int[][] preToPost,
                                int[] numPreToPost,
                                MyPoint locPost,
                                MyPoint[] locPre,
                                double SD,
                                double numSDs,
                                int numToAdd){
        int numPre = numPostFromPre.length;
        double reach = numSDs*SD;
        // create a list of all unchosen presynaptic, along with distances to those neurons
        int[] possiblePre = new int[4096];
        //double[] distPossiblePre = new double[4096];
        int numPossiblePre = 0;
        
        /// create a boolean array for quick determination of prior connectivity
        boolean[] alreadyConnected = new boolean[numPre];
        for (int i=0;i<numPreToPost[post];i++){
            int pre = preToPost[post][i];
            alreadyConnected[pre] = true;
        }        
        /// go throuh all pre, check for distance to post and whether connected--add to list if conditions met
        for (int i=0;i<numPre;i++){    
            double dist = locPost.distance(locPre[i]);
            boolean closeEnough = dist < reach;
            if (closeEnough & !alreadyConnected[i]){
                // then add to list
                possiblePre[numPossiblePre] = i;
                //distPossiblePre[numPossiblePre] = dist;
                numPossiblePre++;
            }
        }
        // list created                
        
        // choose numToAdd neurons from list uniformly
        int[] chosen;
        int actualNumToAdd;
        if (numToAdd > numPossiblePre){
            chosen = new int[numPossiblePre];
            for (int i=0;i<numPossiblePre;i++){
                chosen[i]=i;
            }
            System.out.println("addNPresynaptic: cannot add more neurons.  x="+locPost.x+" y="+locPost.y+". Adding " + numPossiblePre + " instead of " + numToAdd + ".");
        } else {
            chosen = this.chooseKofN(numToAdd, numPossiblePre);
        }
        // add chosen to connectivity
        for (int i=0;i<chosen.length;i++){
            int toAdd = possiblePre[chosen[i]];
            preToPost[post][numPreToPost[post]] = toAdd;
            numPreToPost[post]++;
            postFromPre[toAdd][numPostFromPre[toAdd]] = post;
            numPostFromPre[toAdd]++;
        }
    }
    
    public void deleteNPresynaptic(int post,
                                   int[][] postFromPre,
                                   int[] numPostFromPre,
                                   int[][] preToPost,
                                   int[] numPreToPost,
                                   int numToDel){
        int[] chosen = this.chooseKofN(numToDel, numPreToPost[post]);        
        
        // record which pre to delete
        int[] chosenIndexes = new int[numToDel];
        for (int i=0;i<numToDel;i++){
            chosenIndexes[i] = preToPost[post][chosen[i]];
            preToPost[post][chosen[i]] = -1;
        }
        
        // delete from postFromPre list
        for (int i=0;i<numToDel;i++){
            int pre = chosenIndexes[i];            
            int placeToDelete = 0;
            //find post in the presynaptic cell's list
            while( postFromPre[pre][placeToDelete] != post){
                placeToDelete++;
            }
            //shift list left over deleted spot
            numPostFromPre[pre]--;
            for (int j=placeToDelete;j<numPostFromPre[pre];j++){
                postFromPre[pre][j] = postFromPre[pre][j+1];
            }
        }
        
        // delete from preToPost list
        int ti = 0;
        int si = 0;
        
        while (true){
            while(si < numPreToPost[post] & preToPost[post][si]==-1){
                si++;
            }
            if (si >= numPreToPost[post]){
                break;
            }                
            preToPost[post][ti] = preToPost[post][si];
            ti++;
            si++;
        }
        numPreToPost[post] -= numToDel;
    }

    public void deletePresynaptic(int post,
            int[][] cPostFromPre,
            int[] numPostFromPre,
            int[][] cPreToPost,
            int[] numPreToPost) {
        for (int i = 0; i < numPreToPost[post]; i++) {
            int pre = cPreToPost[post][i];

            ////delete connection
            //find post in pre's connectivity list
            int postsPlace = 0;
            while (cPostFromPre[pre][postsPlace] != post) {
                postsPlace++;
            }

            //shift all others in list left 1
            for (int j = postsPlace; j < numPostFromPre[pre] - 1; j++) {
                cPostFromPre[pre][j] = cPostFromPre[pre][j + 1];
            }
            numPostFromPre[pre]--;

            ////deleted
        }

        numPreToPost[post] = 0;
    }

    public final void choosePresynaptic(int post,
            MyPoint locPost,
            int sizePreGroup,
            MyPoint[] preLocs,
            double lenScale,
            int[] numInputs,
            int[][] cPostFromPre,
            int[] numPostFromPre,
            int[][] cPreToPost,
            int[] numPreToPost,
            double[] egalPlains,
            double desiredNumPrePeak1,
            double[] randomPicks) {
        /*== Choose presynaptic connections ==*/
        double s = (double) numInputs[post] / 8.0;
        double distToBndryX = min(locPost.x, patchWidth - locPost.x);
        double distToBndryY = min(locPost.y, patchWidth - locPost.y);
        double distToBndry = min(distToBndryX, distToBndryY);

        double connectivityWithEgal; //to be calculated in following if-else statement

        if (distToBndry < 1.8 * lenScale) {
            /* then we need to numerically adjust the peak probability to get the expected number presynaptic to be the same */

            /* find all neurons within 2*lenScale */
            double[] allR = new double[3500];
            int numR = 0;
            for (int i = 1; i < sizePreGroup; i++) {
                double sds = preLocs[i].distance(locPost) / lenScale;
                if (sds < 2.0) {
                    allR[numR] = sds;
                    numR++;
                }
            }

            /* map f(x)=e^(-0.5*x^2) over that */
            double[] allP = new double[3500];
            int numP = numR;
            for (int i = 1; i < numR; i++) {
                allP[i] = exp(-0.5 * allR[i] * allR[i]);
            }

            /* use bisection method to find correct connectivity */
            double CLow = egalPlains[numInputs[post]];
            double CHigh = 7.39;
            double expNumPreLow = 0;
            double expNumPreHigh = numR;
            for (int i = 1; i < numP; i++) {
                expNumPreLow += allP[i];
            }
            expNumPreLow *= CLow;

            double threshold = 10; // getting expected number presynaptic within threshold is when the approximation stops

            double desiredExpNumPre;
            if (randomPicks[post] > 0) { //ie if it is a random pick
                desiredExpNumPre = desiredNumPrePeak1 * max(egalPlains[numInputs[post]], randomPicksEEConnectivity);
            } else {
                desiredExpNumPre = desiredNumPrePeak1 * egalPlains[numInputs[post]];
            }

            if (expNumPreHigh < desiredExpNumPre + threshold) {
                /*
                System.out.println("went high");
                System.out.println(expNumPreHigh);
                System.out.println(desiredExpNumPre);
                System.out.println(egalPlains[numInputs[post]]);
                */

                connectivityWithEgal = CHigh;
            } else if (expNumPreLow > desiredExpNumPre - threshold) {
                connectivityWithEgal = CLow;
            } else {
                double CMid = (CHigh + CLow) / 2.0;
                double expNumPreMid = 0;
                for (int i = 1; i < numP; i++) {
                    expNumPreMid += min(1.0, CMid * allP[i]);
                }
                while (expNumPreMid > desiredExpNumPre + threshold
                        || expNumPreMid < desiredExpNumPre - threshold) {
                    if (expNumPreMid < desiredExpNumPre) {
                        CLow = CMid;
                    } else {
                        CHigh = CMid;
                    }
                    CMid = (CLow + CHigh) / 2.0;
                    expNumPreMid = 0;
                    for (int i = 1; i < numP; i++) {
                        expNumPreMid += min(1.0, CMid * allP[i]);
                    }
                }
                connectivityWithEgal = CMid;
                //debugCount++;
                //if (true) System.out.println("expNumPreLow="+expNumPreLow+" desiredExpNumPre="+desiredExpNumPre+" connectivity="+connectivityWithEgal+" expNumPreMid"+expNumPreMid);

            }
        } else {
            if (randomPicks[post] > 0) {//ie it is a random pick
                connectivityWithEgal = max(egalPlains[numInputs[post]], randomPicksEEConnectivity);//(numInputs[post] < egalHelper.length ? egalHelper[numInputs[post]] + connectivity : connectivity);
            } else {
                connectivityWithEgal = egalPlains[numInputs[post]];//(numInputs[post] < egalHelper.length ? egalHelper[numInputs[post]] + connectivity : connectivity);
            }
        }

        /* actually choose connections now that connectivity is established */
        for (int j = 0; j < sizePreGroup; j++) {
            double sds = preLocs[j].distance(locPost) / lenScale;
            if (sds < 2.0 && post != j && rndGenL4.nextDouble() < connectivityWithEgal * exp(-sds * sds / 2)) {
                // then make a connection
                cPostFromPre[j][numPostFromPre[j]] = post; //place on connectivity list
                numPostFromPre[j]++; //update number of E's postsynaptic to j
                cPreToPost[post][numPreToPost[post]] = j; //place on presynaptic connectivity list
                numPreToPost[post]++; //update number of E's presynaptic to i
            }
        }
        /*== chose presynaptic connections ==*/
    }

    public void loadLayer6Recording(){
        try {
            MatFileReader matfilereader;
            if (System.getProperty("os.name").equals("Linux")) {
                matfilereader = new MatFileReader(linPATH.concat(fbInputFile));
            } else if (System.getProperty("os.name").equals("Mac OS X")){
                matfilereader = new MatFileReader(macPATH.concat(fbInputFile));
            } else {
                matfilereader = new MatFileReader(winPATH.concat(fbInputFile));
            }

            /*==Import out-degrees and out-connectivities
             *of LGN cells (these will need further processing,
             *like changing to java indexing starting at 0)==*/
            double[][] fbSpikeTimesPre = ((MLDouble) matfilereader.getMLArray("fbSpikeTimes")).getArray();
            double[][] fbSpikesPre = ((MLDouble) matfilereader.getMLArray("fbSpikes")).getArray();
            double[][] numFBSpikesPre = ((MLDouble) matfilereader.getMLArray("numFBSpikes")).getArray();
            double[][] fbLocPre = ((MLDouble) matfilereader.getMLArray("fbLoc")).getArray();
            double[][] numFBPre = ((MLDouble) matfilereader.getMLArray("numFB")).getArray();

            /* Process */
            numFB = (int) numFBPre[0][0];
            numFBSpikes = (int) numFBSpikesPre[0][0];
            fbSpikeTimes = new double[600000];
            fbSpikes = new int[600000];
            noiseFBSpikeTimes = new double[600000];
            noiseFBSpikes = new int[600000];
            dBFBSpikes = new int[600000];
            dBFBSpikeTimes = new double[600000];
            for (int i = 0; i < numFBSpikes; i++) {
                fbSpikeTimes[i] = fbSpikeTimesPre[i][0];
                fbSpikes[i] = (int) fbSpikesPre[i][0];
            }
            fbSpikeTimes[numFBSpikes] = fbEnd + 1.0;
            locFB = new double[numFB][2];
            for (int i = 0; i < numFB; i++) {
                locFB[i][0] = fbLocPre[i][0];
                locFB[i][1] = fbLocPre[i][1];
            }

            /* setup stats */
            l6SpikesNow = new int[numFB*10];
            l6SpikeTimesNow = new double[numFB*10];
            numL6SpikesNow = 0;
            preRateFB = new double[numFB];            
            
            /* setup the p-value arrays */
            fbPComplex = new double[numFB];
            fbPComplexPast = new double[numFB];
            noiseFBPComplex = new double[numFB];
            timeLastSpike = new double[numFB];
            timeLastNoiseSpike = new double[numFB];
            
            /* setup the time avgs arrays */
            initializeTimeAvgs();
            
            /* setup anti-spike mechanism */
            antiSpikeFB = new double[numFB];

            /* setup spike spread mechanism */
            nextFBSpikeTime = new double[numFB];
            numWaitingFBSpikes = new int[numFB];

        } catch (IOException ex) {
            System.out.println("Error in fb import code");
            System.exit(1);
        }
    }

    public double fBSpFreqDep(double cpd){
        if (cpd <= 2.5){
            return exp(-0.5*sqr((cpd-2.5)/1.33));
        } else {
            return exp(-0.5*sqr((cpd-2.5)/2.25));
        }
    }

    /*==set p-values according to location in cortex==*/
    public void setFBPValues() {
        /* find sf dependence */
        /* for sf dependence, hijack LGN sf dependence and stretch horiz by 1.33 about 2.5c/d*/
        double sfFac = fBSpFreqDep(spFreq);
        //double facPre = exp(-0.5 * sqr((spFreq-2.5)/1.5));
        //double sfFac = facPre;
        //if (facPre < .26) sfFac=.3;
        //if (facPre < .01) sfFac=.1;

        /* find normalizing factor of convolving Gaussian */
        double SD = 75.0;
        double gridWidth = 30;
        double dx = 2.0 * SD / gridWidth;
        double N = 0.0;
        for (double x = -SD; x < SD; x += dx) {
            for (double y = -SD; y < SD; y += dx) {
                double normDSqr = (sqr(x) + sqr(y)) / sqr(SD);
                if (normDSqr < 1.0) {
                    N += exp(-0.5 * normDSqr);
                }
            }
        }

        /* go through each layer 6 neuron and assign a p-value */
        for (int i = 0; i < numFB; i++) {
            /* get the neuron location */
            double[] loc = locFB[i];

            /* do a cheap convolution */
            double conv = 0.0;
            for (double x = loc[0] - SD; x < loc[0] + SD; x += dx) {
                for (double y = loc[1] - SD; y < loc[1] + SD; y += dx) {
                    double dlocx = x - loc[0];
                    double dlocy = y - loc[1];
                    double normDSqr = (sqr(dlocx) + sqr(dlocy)) / sqr(SD);
                    if (normDSqr < 1.0) {
                        double[] newLoc = {x, y};
                        conv += constantPValues(newLoc) * exp(-0.5 * normDSqr);
                    }
                }
            }
            conv /= N;

            /* store the result in the current layer 6 neuron, accounting for sf */
            fbPComplex[i] = sfFac * conv;
        }
    }

    public void connectLayer6(){
        if (loadFBConnections){
            System.out.println("loading feedback connections file "+fBConnectionsFile+"...");
            loadFBConnectivity(fBConnectionsFile);
        } else {
            System.out.println("generating new layer6 to 4Calpha connections....");
            /* allocate fb connectivity matrix */
            numExcPostFB = new int[numFB];
            excPostFB = new int[numFB][800];
            fBPreExc = new int[numExcitatory][250];
            numFBPreExc = new int[numExcitatory];
            fBPreInh = new int[numInhibitory][200];
            numFBPreInh = new int[numInhibitory];
            numSimplePreExc = new int[numExcitatory];

            /* setup simple vs complex layer 6 cells */
            simpleFB = new double[numFB];
            highFiringComplexFB = new double[numFB];
            for (int i = 0; i < numFB; i++) {
                simpleFB[i] = (rndGenL6.nextDouble() < fracSimple ? 1.0 : 0.0);
                highFiringComplexFB[i] = (1.0-simpleFB[i]) * (rndGenL6.nextDouble() < highFiringComplexPerc ? 1.0 : 0.0);
            }
            phaseFB = new double[numFB];
            simpleAmp = new double[numFB];
            deregAmp = new double[numFB];
            for (int i = 0; i < numFB; i++) {
                phaseFB[i] = rndGenL6.nextDouble();
                simpleAmp[i] = simpleAmpLow + (simpleAmpHigh-simpleAmpLow)*rndGenL6.nextDouble();
                deregAmp[i] = deregAmpLow + (deregAmpHigh-deregAmpLow)*rndGenL6.nextDouble();
            }

            final double fbProjectingDensity = (10.0/13.0)*(3000.0/(500.0*500.0))*(1.0/8.0);

            /*== create layer6->4C connectivity matrix */
            for (int i = 0; i < numExcitatory; i++) {
                /* Let user know how much progress is done */
                if (i % 3000 == 0) {
                    System.out.println("Processed " + i + " neurons");
                }

                /* find the peak probability from the egalitarian connectivity */
                //wide
                double expNumFBPreEWide = fbWideEgal[numInputsExc[i]];
                double numFBPreEWideAtPeak1 = fbProjectingDensity*2*PI*(1-exp(-fbWideNumSDsE*fbWideNumSDsE/2)) * fbWideSDE * fbWideSDE;
                double peakWide = expNumFBPreEWide / numFBPreEWideAtPeak1;
                //narrow
                double expNumFBPreENarrow = fbNarrowEgal[numInputsExc[i]];
                double numFBPreENarrowAtPeak1 = fbProjectingDensity*2*PI*(1-exp(-fbNarrowNumSDsE*fbNarrowNumSDsE/2)) * fbNarrowSDE * fbNarrowSDE;
                double peakNarrow = expNumFBPreENarrow / numFBPreENarrowAtPeak1;


                /* find presynaptic layer 6 cells */
                double locX = locExc[i].x;
                double locY = locExc[i].y;
                for (int j = 0; j < numFB; j++) {
                    if (locFB[j][0] == -10000){
                        continue;
                    }
                    double dxSqr = sqr(locX - locFB[j][0]);
                    double dySqr = sqr(locY - locFB[j][1]);
                    double d = sqrt(dxSqr + dySqr);
                    //System.out.println("neuron (" + locX + ", " + locY + ") and fb ("+locFB[j][0]+", "+locFB[j][1]+") dist sqr normalized is " + dSqr);
                    //try{Thread.sleep(250);}catch(InterruptedException e){System.out.println("caught");}
                    if (d < fbWideNumSDsE * fbWideSDE) {
                        double dWideSDs = d / fbWideSDE;
                        double probPreSyn = peakWide * exp(-dWideSDs*dWideSDs / 2);
                        if (d < fbNarrowNumSDsE * fbNarrowSDE){
                            double dNarrowSDs = d / fbNarrowSDE;
                            probPreSyn += peakNarrow * exp(-dNarrowSDs*dNarrowSDs / 2);
                        }
                        //System.out.println("prob="+probPreSyn);
                        if (rndGenL6.nextDouble() < probPreSyn) {
                            /* make a connection from jth FB cell to ith E cell*/
                            //System.out.println("success");
                            excPostFB[j][numExcPostFB[j]] = i;
                            numExcPostFB[j]++;

                            fBPreExc[i][numFBPreExc[i]] = j;
                            numFBPreExc[i]++;
                            if (simpleFB[j] > 0.5) numSimplePreExc[i]++;
                        }
                    }
                }
            } //run through E cells

            //go through each numlgn pop and find outlier info
            double[] mean = new double[7];
            double[] meanOfSquare = new double[7];
            double[] SD = new double[7];

            for (int i=0;i<numExcitatory;i++){
                int nLGNIn = numInputsExc[i];
                mean[nLGNIn] += (double) numFBPreExc[i];
                meanOfSquare[nLGNIn] += sqr( (double)numFBPreExc[i] );
            }
            for (int i=0;i<7;i++){
                mean[i] /= numExcitatory;
                meanOfSquare[i] /= numExcitatory;
                SD[i] = meanOfSquare[i] - mean[i] * mean[i];
            }

            int[] faultyList = new int[numExcitatory];
            int numFaulty = 0;
            //go through excitatory cells, marking anyone who is not with 1/2 sd of mean
            for (int i=0;i<numExcitatory;i++){
                int nLGNIn = numInputsExc[i];
                double simpleSD = Math.sqrt((double)numFBPreExc[i]*fracSimple*(1.0-fracSimple));
                double simpleMean = (double)numFBPreExc[i] * fracSimple;
                if (numFBPreExc[i] < mean[nLGNIn]-0.5*SD[nLGNIn]
                        || numFBPreExc[i] > mean[nLGNIn]+0.5*SD[nLGNIn]
                        || numSimplePreExc[i] < simpleMean-simpleSD
                        || numSimplePreExc[i] > simpleMean+simpleSD){
                    faultyList[numFaulty] = i;
                    numFaulty++;
                }
            }
            /*delete connections*/
            for (int i=0;i<numFaulty;i++){
                /*delete connections*/
                int l4Index = faultyList[i];
                for (int j=0;j<numFBPreExc[l4Index];j++){
                    int fBIndex = fBPreExc[l4Index][j];
                    /*delete l4Index from fBIndex's excPostFB*/
                    int k=0;
                    while(excPostFB[fBIndex][k]!=l4Index && k < numFB) k++;
                    /*k is now at the index where we want to delete*/
                    while(k<numExcPostFB[fBIndex]){
                        excPostFB[fBIndex][k] = excPostFB[fBIndex][k+1];
                        k++;
                    }
                    numExcPostFB[fBIndex]--;
                }
                numFBPreExc[l4Index]=0;
                numSimplePreExc[l4Index]=0;
            }
            /*redraw connections*/
            for (int iFaulty=0;iFaulty<numFaulty;iFaulty++){
                int i = faultyList[iFaulty];
                /* Let user know how much progress is done */
                if (i % 3000 == 0) {
                    System.out.println("Processed " + i + " neurons");
                }

                /* find the peak probability from the egalitarian connectivity */
                double expNumFBPreEWide = fbWideEgal[numInputsExc[i]];
                double numFBPreEWideAtPeak1 = fbProjectingDensity*2*PI*(1-exp(-fbWideNumSDsE*fbWideNumSDsE/2)) * fbWideSDE * fbWideSDE;
                double peakWide = expNumFBPreEWide / numFBPreEWideAtPeak1;
                //narrow
                double expNumFBPreENarrow = fbNarrowEgal[numInputsExc[i]];
                double numFBPreENarrowAtPeak1 = fbProjectingDensity*2*PI*(1-exp(-fbNarrowNumSDsE*fbNarrowNumSDsE/2)) * fbNarrowSDE * fbNarrowSDE;
                double peakNarrow = expNumFBPreENarrow / numFBPreENarrowAtPeak1;

                /* find presynaptic layer 6 cells */
                double locX = locExc[i].x;
                double locY = locExc[i].y;
                for (int j = 0; j < numFB; j++) {
                    if (locFB[j][0] == -10000){
                        continue;
                    }
                    double dxSqr = sqr(locX - locFB[j][0]);
                    double dySqr = sqr(locY - locFB[j][1]);
                    double d = sqrt(dxSqr + dySqr);
                    //System.out.println("neuron (" + locX + ", " + locY + ") and fb ("+locFB[j][0]+", "+locFB[j][1]+") dist sqr normalized is " + dSqr);
                    //try{Thread.sleep(250);}catch(InterruptedException e){System.out.println("caught");}
                    if (d < fbWideNumSDsE * fbWideSDE) {
                        double dWideSDs = d / fbWideSDE;
                        double probPreSyn = peakWide * exp(-dWideSDs*dWideSDs / 2);
                        if (d < fbNarrowNumSDsE * fbNarrowSDE){
                            double dNarrowSDs = d / fbNarrowSDE;
                            probPreSyn += peakNarrow * exp(-dNarrowSDs*dNarrowSDs / 2);
                        }
                        //System.out.println("prob="+probPreSyn);
                        if (rndGenL6.nextDouble() < probPreSyn) {
                            /* make a connection from jth FB cell to ith E cell*/
                            //System.out.println("success");
                            excPostFB[j][numExcPostFB[j]] = i;
                            numExcPostFB[j]++;

                            fBPreExc[i][numFBPreExc[i]] = j;
                            numFBPreExc[i]++;
                            if (simpleFB[j] > 0.5) numSimplePreExc[i]++;
                        }
                    }
                }
            }/*redrew connections*/

            /* allocate fb connectivity matrix to I's */
            numInhPostFB = new int[numFB];
            inhPostFB = new int[numFB][700];
            
            /* find the peak probability */
            double expNumFBPreIWide = fbIExpNumPreWide;
            final double numFBPreIWideAtPeak1 = fbProjectingDensity*2*PI*(1-exp(-fbWideNumSDsI*fbWideNumSDsI/2)) * fbWideSDI * fbWideSDI;
            final double peakWideGuess = expNumFBPreIWide / numFBPreIWideAtPeak1;
            //narrow
            double expNumFBPreINarrow = fbIExpNumPreNarrow;
            double numFBPreINarrowAtPeak1 = fbProjectingDensity*2*PI*(1-exp(-fbNarrowNumSDsI*fbNarrowNumSDsI/2)) * fbNarrowSDI * fbNarrowSDI;
            final double peakNarrowGuess = expNumFBPreINarrow / numFBPreINarrowAtPeak1;
            
            //calculate correction for hitting ceiling
            double dx = locFB[1][0]-locFB[0][0];
            double width = 2.0*fbWideNumSDsI*fbWideSDI;
            final MyPoint center = new MyPoint(width/2,width/2);
            MyPoint[] neurons = makeNeuronArray(width,dx);
            final MyPoint[] neuronsInRange = selectNeuronsInRange(neurons, center, width/2);
            RealFunc F = new RealFunc(){
                @Override
                public double eval(double p){
                    double s=0;
                    for (MyPoint neuron : neuronsInRange){
                        double d = neuron.distance(center);
                        double dw = d / fbWideSDI;
                        if (d < fbNarrowSDI*fbNarrowNumSDsI){
                            double dn = d / fbNarrowSDI;
                            s += min(1.0,
                                 p*( peakWideGuess*exp(-dw*dw/2.0) + peakNarrowGuess*exp(-dn*dn/2.0) )
                                );
                        } else {
                            s += min(1.0,
                                 p*( peakWideGuess*exp(-dw*dw/2.0) )
                                );
                        }
                        
                    }
                    return s;
                }
            };
            double target=expNumFBPreIWide+expNumFBPreINarrow;
            double x1=.9;
            double x2=1.0;
            double y1=F.eval(x1);
            double y2=F.eval(x2);
            double epsilon=1;
            NewtonsMethod nm = new NewtonsMethod(target, x1, x2, y1, y2, epsilon);
            nm.F=F;
            double peakAdjustment = nm.solveForX();
            double peakWide = peakWideGuess*peakAdjustment;
            double peakNarrow = peakNarrowGuess*peakAdjustment;
            

            /* run through I cells */
            for (int i = 0; i < numInhibitory; i++) {
                /* Let user know how much progress is done */
                if (i % 3000 == 0) {
                    System.out.println("Processed " + i + " neurons");
                }
                

                /* find presynaptic layer 6 cells */
                double locX = locInh[i].x;
                double locY = locInh[i].y;
                for (int j = 0; j < numFB; j++) {
                    if (locFB[j][0] == -10000){
                        continue;
                    }
                    double dxSqr = sqr(locX - locFB[j][0]);
                    double dySqr = sqr(locY - locFB[j][1]);
                    double dSDs = sqrt(dxSqr + dySqr) / fbWideSDI;
                    if (dSDs < fbWideNumSDsI ) {
                        double probPreSyn = peakWide * exp(-dSDs*dSDs / 2);
                        if (dSDs < fbNarrowNumSDsI) {
                            probPreSyn += peakNarrow * exp(-dSDs*dSDs / 2);
                        }
                        if (rndGenL6.nextDouble() < probPreSyn) {
                            /* make a connection from jth FB cell to ith E cell*/
                            inhPostFB[j][numInhPostFB[j]] = i;
                            numInhPostFB[j]++;
                            
                            fBPreInh[i][numFBPreInh[i]] = j;
                            numFBPreInh[i]++;
                        }
                    }
                }
            }
        }
    }
    
    public MyPoint[] makeNeuronArray(double width, double dx){
        int n = (int) ceil(width / dx);
        MyPoint[] neuronArray = new MyPoint[n*n];
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                neuronArray[i*n+j] = new MyPoint(dx*i,dx*j);
            }
        }
        
        return neuronArray;
    }
    
    public MyPoint[] selectNeuronsInRange(MyPoint[] neurons, MyPoint center, double range){
        //find distances and count number in range
        double[] distances = new double[neurons.length];
        int numInRange = 0;
        for (int i=0;i<neurons.length;i++){
            distances[i] = center.distance(neurons[i]);
            if (distances[i] < range){
                numInRange++;
            }
        }
        
        //add all in range to new array
        MyPoint[] inRange = new MyPoint[numInRange];
        int numPutIn=0;
        for (int i=0;i<neurons.length;i++){
            if (distances[i] < range){
                inRange[numPutIn]=neurons[i];
                numPutIn++;
            }
        }
        
        return inRange;
    }
            

    public double constantPValues(double[] loc) {
        double constantP;
        /* switch to 1st HC coordinates */
        double[] locHC1 = {NetHelper.movCoordToHC1(loc[0], HCWidth),
            NetHelper.movCoordToHC1(loc[1], HCWidth)};

        double locLen = MyPoint.distance(0, 0, locHC1[0], locHC1[1]);

        /* find constant value based on angle */
        double angleHC1 = NetHelper.pointToAngle(locHC1[0], locHC1[1]);
        //double orientAngle = NetHelper.pointToAngle(kOrient.x,kOrient.y);

        int ind = ((int) floor((angleHC1 + 2.0 * this.angle) / (PI / 6.0))) % 12;

        constantP = pValues[ind];

        return constantP;
    }

    /* set LGN strength */
    public void setSELGNOverSEE(double newS) {
        sELGNOverSEE = newS;
        enforceRatios();
    }

    public void setSILGNOverSEE(double newS) {
        sILGNOverSEE = newS;
        enforceRatios();
    }

    public void setSEE(double newSEE) {
        sEE = newSEE;
        enforceRatios();
    }

    public void enforceRatios(){
        sLGN = sELGNOverSEE * sEE * aChSELGNFac;
        sLGNInh = sILGNOverSEE * sEE;
        sEI = sEIOverSEE * sEE;
        sIE = sIEOverSEE * sEE;
        sII = sEI;
        sIEL6 = sIEOverSEE*sIEL6OverSEEL6FracL4 * (sEEL6High+sEEL6Low)/2.0;
        setSIIOverSEIRange();
    }

    public void engageGrating(){
        setOrientationByAngle(angle); //sets phases and angle variable
        //get sf dependence for lgn
        spFreqMod= lGNSpFreqDep(spFreq);
        //set p-values
        setFBPValues();
        //update peak LGN coupled min L6
        setPeakLGNMinL6();
    }

    public void setAChSELGNFac(double newFac){
        aChSELGNFac = newFac;
        enforceRatios();
    }


    /* set the ratio of SEI to SEE */
    public void setSEIOverSEE(double newSEIOverSEE) {
        sEIOverSEE = newSEIOverSEE;
        enforceRatios();
    }

    public void setSIEOverSEE(double newSIEOverSEE) {
        sIEOverSEE = newSIEOverSEE;
        enforceRatios();
    }

    /* set the ambient drive parameters */
    public void setAmbient(double rateE, double rateI) {
        rateEAmb = rateE;
        rateIAmb = rateI;
    }

    /* set the I threshold range */
    public void setThreshIRange(double low, double high) {
        for (int i = 0; i < numInhibitory; i++) {
            threshI[i] = low + (high - low) * rndGenMisc.nextDouble();
        }
    }

    public void setThreshERange(double low, double high) {
        for (int i = 0; i < numExcitatory; i++) {
            threshE[i] = low + rndGenMisc.nextDouble() * (high-low);
        }
    }

    public void setSModSeeds() {
        for (int i = 0; i < numExcitatory; i++) {
            sModSeedExc[i] = rndGenMisc.nextDouble();
        }
        for (int i = 0; i < numInhibitory; i++) {
            sModSeedInh[i] = rndGenMisc.nextDouble();
        }
        for (int i = 0; i < numLGN; i++) {
            sModSeedLGNON[i] = rndGenMisc.nextDouble();
            sModSeedLGNOFF[i] = rndGenMisc.nextDouble();
        }
    }

    /* set the range of possible I to I strength modifiers */
    public void setSIIOverSEIRange() {
        sIIOverSEIRange = sIIOverSEIHigh - sIIOverSEILow;
    }

    /* set light grating contrast */
    public void setContrast(double newContrast) {
        contrast = newContrast;
    }

    public void setIntensity(double newIntensity) {
        I_0 = newIntensity;
    }

    public void setSpFreq(double newSpFreq) {
        spFreq = newSpFreq;
        setPhases();
        spFreqMod= lGNSpFreqDep(spFreq);
    }

    public void setTF(double newTF){ //TF=temporal frequency
        tf = newTF;

    }

    public void setPhases() {
        double sFRadiansPerMicron = 2*PI*spFreq / 2000.0;
        for (int i = 0; i < numLGN; i++) {
            phaseON[i] = sFRadiansPerMicron * orientation.innerProduct(locLGNON[i]);
            phaseOFF[i] = sFRadiansPerMicron * orientation.innerProduct(locLGNOFF[i]) + PI;
        }
    }

    public double lGNSpFreqDep(double cpdPre){
        double sa = 0.066 * 2.0 * 1.2;
        double sb = 0.093 * 2.0 * 1.2;
        double cpd;
        if (cpdPre > 2.2){
            cpd = 2.2+(cpdPre-2.2)/lGNSpFreqRightStretch;
        } else {
            cpd = cpdPre; 
        }
        return max(
                0.0,
                0.95 * (1 / 0.3359) * (1.0 * exp(-sa * sa * cpd * cpd * PI) - 0.74 * exp(-sb * sb * cpd * cpd * PI))
        );
    }

    public void setOrientation(MyPoint newOrient) {
        orientation = newOrient;
        angle = NetHelper.pointToAngle(orientation.x, orientation.y);
        setPhases();
    }

    public void setOrientationByAngle(double angle) {
        this.angle = angle;
        orientation = new MyPoint(Math.cos(angle),Math.sin(angle));
        setPhases();
    }

    /* temporal frequency of drift grating */
    public void setFrequency(double newW) {
        tf = newW;
    }

    /* resets the LGN inputs if a different LGN input file is given--also rebuilds the network */
    public void setLGNInputsIfNew(String LGNInputFile) {
        System.out.println("Code for setLGNInputsIfNew commented out, not doing anything");
        /*
         if (!this.LGNInputFile.equals(LGNInputFile)){
         this.LGNInputFile = LGNInputFile;

         System.out.println("New LGN Input File detected. Rebuilding Network...");
         getLGNInputs();
         setupCortexNetwork(); //need to also reset the cortical network bc egal can be different now
         setSModSeeds();
         }
         */
    }
    
    public void evolveExponential(double tau, double[] y, double numExponentials){
        for (int i=0; i<numExponentials; i++){
            y[i] *= exp(-dt / tau);
        }
    }
    
    public void evolveDifOfExponentials(double taua, double taub, double[][] g, double numNeurons){
        for (int i=0; i<numNeurons; i++){
            g[i][1] *= exp(-dt / taua);
            g[i][2] *= exp(-dt / taub);
            g[i][0] = g[i][1] - g[i][2];
        }
    }
    
    public void evolveTwoDifOfExponentials(double taua, double taub, double tauc, double taud, double[][] g, double numNeurons){
        for (int i=0; i<numNeurons; i++){
            g[i][1] *= exp(-dt / taua);
            g[i][2] *= exp(-dt / taub);
            g[i][3] *= exp(-dt / tauc);
            g[i][4] *= exp(-dt / taud);
            g[i][0] = (g[i][1] - g[i][2]) + (g[i][3]-g[i][4]);     //weighting happens in synapse
        }
    }
    
    public void synapseDifOfExponentials(double taua, double taub, double[] g, double spaceK, double S, double tStart){
        double area = taua - taub;
        double fac = S * spaceK / area;
        g[1] += fac;
        g[2] += fac;        
    }
    
    public void synapseTwoDifOfExponentials(double taua, double taub, double tauc, double taud, double weightToSecond, double[] g, double spaceK, double S, double tStart){
        double areaFirst = taua - taub;
        double facFirst = (1.0-weightToSecond) * S * spaceK / areaFirst;
        double areaSecond = tauc - taud;
        double facSecond = weightToSecond * S * spaceK / areaSecond;
        g[1] += facFirst;
        g[2] += facFirst;        
        g[3] += facSecond;
        g[4] += facSecond;
    }

    /* calculate the conductances for the next timestep */
    public void evolveAMPAGABA(double tau, double[][] g, double numConductances, double dt) {
        for (int i = 0; i < numConductances; i++) {
            /* main conductance is g[i][0] */
            g[i][0] = exp(-dt / tau) * (g[i][0] + (dt / tau) * (3 * g[i][1] + (dt / tau) * (3 * g[i][2] + (dt / tau) * g[i][3])));
            g[i][1] = exp(-dt / tau) * (g[i][1] + (dt / tau) * (2 * g[i][2] + (dt / tau) * g[i][3]));
            g[i][2] = exp(-dt / tau) * (g[i][2] + (dt / tau) * g[i][3]);
            g[i][3] *= exp(-dt / tau);
        }
    }

    /* set the conductances after being synapsed on */
    public double synapseAMPAGABA(double tau, double[] g, double spaceK, double S, double tStart) {
        double tStartScaled = tStart / tau;
        //double shiftFac = exp(-tStartScaled)*(tStartScaled*(tStartScaled*(tStartScaled+3.0)+6.0)+6.0)/6.0;
        double accum = S * spaceK * exp(-tStart / tau) / (6.0 * tau);// / shiftFac;
        g[3] += accum;
        accum *= tStartScaled;
        g[2] += accum;
        accum *= tStartScaled;
        g[1] += accum;
        accum *= tStartScaled;
        g[0] += accum;
        return accum;
    }

    /* same as above except for NMDA */
    public void evolveNMDA(double[][] g, double numConductances) {
        for (int i = 0; i < numConductances; i++) {
            /* main conductance is g[i][0] */
            g[i][1] *= exp(-dt / tauNMDA1);
            g[i][2] *= exp(-dt / tauNMDA2);
            g[i][0] = g[i][1] - g[i][2];
        }
    }

    public void synapseNMDA(double[] g, double spaceK, double S) {
        double coeff = S * spaceK / (tauNMDA1 - tauNMDA2);
        g[1] += coeff;
        g[2] += coeff;
    }

    public final double sqr(double x) {
        return x * x;
    }

    public final double gaussian2D(double sigma, double normSqr) {
        if (normSqr > 4 * sigma * sigma) {
            return 0.0;
        } else {
            return 1.0 / (2.0 * PI * sqr(sigma)) * exp(-normSqr / (2.0 * sqr(sigma)));
        }
    }

    /* find the voltage for the next timestep */
    public void evolveVExc() {
        for (int i = 0; i < numExcitatory; i++) {
            if (sleepExc[i] >= 0) {
                vExc[i] = vExc[i] * exp(-dt * gLExc) + 
                        dt * (
                               ( gAmbExc[i][0] + 
                                 gLGNExc[i][0] + 
                                 //(dtExc/dt)*(gAMPAExc[i][0]-initialAMPADamperExc[i]) + 
                                 gAMPAExc[i][0] + 
                                 gNMDAExc[i][0] + 
                                 //(dtExc/dt)*(gAMPAFBExc[i][0]-initialAMPADamperExcFB[i]) + 
                                 gAMPAFBExc[i][0] + 
                                 gNMDAFBExc[i][0]
                               ) * (eRev - vExc[i]) 
                               + 
                               //(dtInh/dt)*(gGABAExc[i][0]-initialGABADamperExc[i]) * (iRev - vExc[i])
                               gGABAExc[i][0] * (iRev - vExc[i])
                             );
            }
        }
    }

    public void evolveVInh() {
        for (int i = 0; i < numInhibitory; i++) {
            if (sleepInh[i] >= 0) {
                vInh[i] = vInh[i] * exp(-dt * gLInh) + 
                        dt * (
                               ( gAmbInh[i][0] + 
                                 gLGNInh[i][0] + 
                                 (dtExc/dt)*(gAMPAInh[i][0]-initialAMPADamperInh[i]) + 
                                 //gAMPAInh[i][0] +
                                 gNMDAInh[i][0] + 
                                 (dtExc/dt)*(gAMPAFBInh[i][0]-initialAMPADamperInhFB[i]) + 
                                 //gAMPAFBInh[i][0] +
                                 gNMDAFBInh[i][0]
                               ) * (eRev - vInh[i]) 
                               + 
                               (dtInh/dt)*(gGABAInh[i][0]-initialGABADamperInh[i]) * (iRev - vInh[i])
                               //gGABAInh[i][0] * (iRev - vInh[i])
                             );
            }
        }
    }

    /* move the refractory clocks forward one */
    public void evolveSleep() {
        for (int i = 0; i < numExcitatory; i++) {
            sleepExc[i]++;
        }
        for (int i = 0; i < numInhibitory; i++) {
            sleepInh[i]++;
        }
    }
    
    public double targetModifier(int target){
        double[] loc = { locExc[target].x, locExc[target].y };
        double[] locHC1 = {NetHelper.movCoordToHC1(loc[0], HCWidth),
            NetHelper.movCoordToHC1(loc[1], HCWidth)};

        /* find constant value based on angle */
        double angleHC1 = NetHelper.pointToAngle(locHC1[0], locHC1[1]);
        //double orientAngle = NetHelper.pointToAngle(kOrient.x,kOrient.y);

        int ind = ((int) floor((angleHC1+PI/6.0) / (PI / 3.0))) % 2;
        
        double[] targetMods = { vertTemplateMod, horizTemplateMod };
        
        return targetMods[ind];
    }
    
    public void setTargetModifiers(){
        targetMod = new double[numExcitatory];
        for (int i=0;i<numExcitatory;i++){
            targetMod[i] = targetModifier(i);
        }
    }
    
    public double getFacilitationFac(double l4Rate){
        int i = 1;
        while(facxs[i] < l4Rate){
            i++;
        }
        i--;
        return facys[i] + facslopes[i] * (l4Rate-facxs[i]);
    }

    public void evolveLGN() {
        /*stats*/
        numSpikesNowLGNON = 0;
        numSpikesNowLGNOFF = 0;

        for (int i = 0; i < numLGN; i++) {


            /* update the voltage */

            double kickDir = 2 * round(rndGenRun.nextDouble()) - 1;
            vLGNON[i] = vLGNON[i] * exp(-dt * leakLGN) + dt * I(t, phaseON[i]) + kickDir * shotNoise * (double) poisson.nextInt(dt * noiseRate);

            kickDir = 2 * round(rndGenRun.nextDouble()) - 1;
            vLGNOFF[i] = vLGNOFF[i] * exp(-dt * leakLGN) + dt * I(t, phaseOFF[i]) + kickDir * shotNoise * (double) poisson.nextInt(dt * noiseRate);

            /* check if a spike occured */
            if (vLGNON[i] > threshLGN) {
                /* spike occured, send it out */
                for (int j = 0; j < nELGNON[i]; j++) {
                    //synapseAMPAGABA(tauAMPA, gLGNExc[cELGNON[i][j]], 1.0, sLGN * modifier(sModSeedLGNON[i],sModSeedExc[i]));
                    int target = cELGNON[i][j];
                    //delayCylinder[target][(spikeDelayIndex + (int)floor(random() * lgnSpikeDelaySteps)) % spikeDelaySteps] += 1.2 * sELGNOverSEE * modifier(sModSeedLGNON[i],sModSeedExc[target]);
                    
                    synapseDifOfExponentials(tauAMPAToEa, tauAMPAToEb, gLGNExc[target], 1.0, sEE * ((numInputsExc[target] > 7 ? sELGNOverSEERnd : sELGNOverSEE)) * modifier(sModSeedLGNON[i], sModSeedExc[target]), dt);
                }
                for (int j = 0; j < nILGNON[i]; j++) {
                    synapseDifOfExponentials(tauAMPAToEa,tauAMPAToEb, gLGNInh[cILGNON[i][j]], 1.0, sLGNInh, dt);
                }

                /* record spike */
                spikesNowLGNON[numSpikesNowLGNON] = i;
                spikeTimesNowLGNON[numSpikesNowLGNON] = t;
                numSpikesNowLGNON++;

                /* reset voltage, no refractory */
                vLGNON[i] = 0;
            }

            if (vLGNOFF[i] > threshLGN) {
                for (int j = 0; j < nELGNOFF[i]; j++) {
                    //synapseAMPAGABA(tauAMPA, gLGNExc[cELGNOFF[i][j]], 1.0, sLGN);
                    int target = cELGNOFF[i][j];
                    //delayCylinder[target][(spikeDelayIndex + (int)floor(random() * lgnSpikeDelaySteps)) % spikeDelaySteps] += 1.2 * sELGNOverSEE * modifier(sModSeedLGNOFF[i],sModSeedExc[target]);
                    synapseDifOfExponentials(tauAMPAToEa, tauAMPAToEb, gLGNExc[target], 1.0, sEE * sELGNOverSEE * modifier(sModSeedLGNOFF[i], sModSeedExc[target]), dt);
                }
                for (int j = 0; j < nILGNOFF[i]; j++) {
                    synapseDifOfExponentials(tauAMPAToEa, tauAMPAToEb, gLGNInh[cILGNOFF[i][j]], 1.0, sLGNInh, dt);
                }

                /* record spike */
                spikesNowLGNOFF[numSpikesNowLGNOFF] = i;
                spikeTimesNowLGNOFF[numSpikesNowLGNOFF] = t;
                numSpikesNowLGNOFF++;

                /* reset voltage, no refractory */
                vLGNOFF[i] = 0;
            }
        } /* for all LGN */

    } /* evolve LGN */


    public double I(double t, double phase) {
        double r = 2 * PI * tf * t + phase;
        double s = sin(r);

        if ((r / PI + 0.5) % 2.00 > 1.0) {
            return I_0 * (1 + contrast * spFreqMod * s);
        } else {
            return I_0 * (1 + contrast * spFreqMod * signum(s) * sqrt(abs(s)));
        }
    }

    /* find all conductances for next timestep */
    public void evolveConductances() {
        evolveDifOfExponentials(tauAMPAToEa,tauAMPAToEb, gAMPAExc, numExcitatory);
        evolveDifOfExponentials(tauAMPAToEa,tauAMPAToEb, gAMPAFBExc, numExcitatory);
        evolveTwoDifOfExponentials(tauGABAa,tauGABAb,tauGABAc,tauGABAd, gGABAExc, numExcitatory);
        //evolveExponential(tauInitialGABADamper, initialGABADamperExc, numExcitatory);
        //evolveExponential(tauInitialAMPADamper, initialAMPADamperExc, numExcitatory);
        //evolveExponential(tauInitialAMPADamper, initialAMPADamperExcFB, numExcitatory);
        evolveNMDA(gNMDAExc, numExcitatory);
        evolveNMDA(gNMDAFBExc, numExcitatory);

        evolveDifOfExponentials(tauAMPAToIa,tauAMPAToIb, gAMPAInh, numInhibitory);
        evolveDifOfExponentials(tauAMPAToIa,tauAMPAToIb, gAMPAFBInh, numInhibitory);
        evolveTwoDifOfExponentials(tauGABAa,tauGABAb,tauGABAc,tauGABAd, gGABAInh, numInhibitory);
        //evolveExponential(tauInitialGABADamper, initialGABADamperInh, numInhibitory);
        //evolveExponential(tauInitialAMPADamper, initialAMPADamperInh, numInhibitory);
        //evolveExponential(tauInitialAMPADamper, initialAMPADamperInhFB, numInhibitory);
        evolveNMDA(gNMDAInh, numInhibitory);
        evolveNMDA(gNMDAFBInh, numInhibitory);

        evolveDifOfExponentials(tauAmbToEa,tauAmbToEb, gAmbExc, numExcitatory);
        evolveDifOfExponentials(tauAmbToIa,tauAmbToIb, gAmbInh, numInhibitory);

        evolveDifOfExponentials(tauAMPAToEa,tauAMPAToEb, gLGNExc, numExcitatory);
        evolveDifOfExponentials(tauAMPAToEa,tauAMPAToEb, gLGNInh, numInhibitory);
    }

    /* implement layer 6 recording feedback */
    public void enactFB() {
        numL6SpikesNow = 0;

        /* background portion of fb */
        /* average fbPComplex=___.  Convert whatever that is to the prob desired */
        /* when avg fbPComplex=1, then l6 recording fires normally, at 30 sp/sec */
        /* so avg fbPComplex * 30 = l6 fr per sec and avgfbPComplex*30/l6L4FRRatio=l4fr */
        /* so avgfbPcomplex*30/l6L4FRRatio*l6L4FRRatioPoisson = l6 poisson rate */
        /* now l6L4FRRatioPoisson = desired l6 rate / bgRecRate */
        for (int i = 0 ; i<numFB; i++){
            double meanL4Rate;
                                
            //double meanPresent = 2.0 * fbPComplex[i] * exp((timeLastSpike[i] - t) / sigmaFB);
            double meanFRPre = 2.0 * getL6L4FRRatioPoisson( timeAvgs[i][whichAvg] );
            for (int j=0;j<12;j++){
                int l6Far = longRangeCoupling[i][j];
                //meanPresent += fbPComplex[l6Far] * exp((timeLastSpike[l6Far] - t) / sigmaFB);
                meanFRPre += getL6L4FRRatioPoisson( timeAvgs[l6Far][whichAvg] );
            }
            meanFRPre /= 14.0;
            
            double poissonRate = meanFRPre;
            double spikeProb = dt * poissonRate / 1000.0;
            if (rndGenRun.nextDouble() < spikeProb){
                // record spike
                registerFBSpike(i);
                
                for (int j=0;j<numExcPostFB[i];j++){
                    spikeEEFB(excPostFB[i][j], 1.0);
                }
                for (int j=0;j<numInhPostFB[i];j++){
                    int target = inhPostFB[i][j];                    
                    spikeIEFB(target, 1.0);
                }
            }
        }
        /* bg portion of fb over */
        
        
        while (fbClock > fbSpikeTimes[currentFBSpike]) {
            /* hit every neuron postsynaptic to fbSpike[currentFBSpike] */
            int firingNeuron = fbSpikes[currentFBSpike];

            //update fmax of firing neuron
            double p;
            {
                double meanL4Rate;
                
                double meanFRBothPre = 2.0 * getL6FRBoth( timeAvgs[firingNeuron][whichAvg] );
                double meanFRPoisPre = 2.0 * getL6L4FRRatioPoisson( timeAvgs[firingNeuron][whichAvg] );                
                for (int j=0;j<12;j++){
                    int l6Far = longRangeCoupling[firingNeuron][j];
                    meanFRBothPre += getL6FRBoth( timeAvgs[l6Far][whichAvg] );
                    meanFRPoisPre += getL6L4FRRatioPoisson( timeAvgs[l6Far][whichAvg] );
                }
                meanFRBothPre /= 14.0;
                meanFRPoisPre /= 14.0;

                double l6Rate = meanFRBothPre - meanFRPoisPre;
                double recordingRate = 27.0;
                p = l6Rate / recordingRate;
            }
            //p *= pow(contrast,fbContrastExponent);
            
            int numGuaranteedSpikes = (int) min( floor(p), 8.0 );

            int numSpikes = (rndGenRun.nextDouble() < p - numGuaranteedSpikes ? numGuaranteedSpikes + 1 : numGuaranteedSpikes);
            int actualNumSpikes = 0;
            int s=0;
            for (int i=0;i<numSpikes;i++){
                if (i > fbSpikeSpreadProbChanges[s]){
                    s++;
                }

                if (rndGenRun.nextDouble() < fbSpikeSpreadProbs[s]){
                    actualNumSpikes++;
                }
            }
            numSpikes = actualNumSpikes;

            if (numSpikes > 0){

                // record spike
                registerFBSpike(firingNeuron);
                
                //facilitation and spike code here
                double meanL4Rate;
                                
                double meanPresent = 2.0 * getL6FRBoth( timeAvgs[firingNeuron][whichAvg] );
                for (int j=0;j<12;j++){
                    int l6Far = longRangeCoupling[firingNeuron][j];
                    meanPresent += getL6FRBoth( timeAvgs[l6Far][whichAvg] );
                }
                meanPresent /= 14.0;
                meanL4Rate = meanPresent;
                
                double facilitationFac = getFacilitationFac(meanL4Rate); //min(facilitationCap, 1.0 + facilitationSlope * pow(max(0.0, meanL4Rate - facilitationStart),facilitationExponent) );
                //double facilitationFac = facilitationCap;
                for (int k = 0; k < numExcPostFB[firingNeuron]; k++) {
                    /* spike on kth postsynaptic neuron */
                    int target = excPostFB[firingNeuron][k];
                    
                    delayCylinderFB[target][(int) floor(rndGenRun.nextDouble() * spikeDelaySteps)] += facilitationFac;
                }
                for (int k = 0; k < numInhPostFB[firingNeuron]; k++) {
                    /* spike on kth postsynaptic neuron */
                    int target = inhPostFB[firingNeuron][k];                                        
                    spikeIEFB(target, ((1.0-facIOverFacE) + facilitationFac*facIOverFacE));
                }

                numWaitingFBSpikes[firingNeuron] += numSpikes-1;
                nextFBSpikeTime[firingNeuron] = t + fbSpreadLow + (fbSpreadHigh-fbSpreadLow) * rndGenRun.nextDouble();
            }

            /* update currentFBSpike */
            currentFBSpike++;
        }

        /* go through all l6 cells and detect waiting spikes */
        for (int i=0;i<numFB;i++){
            if (t > nextFBSpikeTime[i] && numWaitingFBSpikes[i] > 0){
                // record spike
                registerFBSpike(i);

                //facilitation and spike code here
                double meanL4Rate;
                                
                double meanPresent = 2.0 * getL6FRBoth( timeAvgs[i][whichAvg] );
                for (int j=0;j<12;j++){
                    int l6Far = longRangeCoupling[i][j];
                    meanPresent += getL6FRBoth( timeAvgs[l6Far][whichAvg] );
                }
                meanPresent /= 14.0;
                meanL4Rate = meanPresent;
                
                double facilitationFac = getFacilitationFac(meanL4Rate); //min(facilitationCap, 1.0 + facilitationSlope * pow(max(0.0, meanL4Rate - facilitationStart),facilitationExponent) );
                //double facilitationFac = facilitationCap;
                for (int k = 0; k < numExcPostFB[i]; k++) {
                    /* spike on kth postsynaptic neuron */
                    int target = excPostFB[i][k];
                    
                    delayCylinderFB[target][(int) floor(rndGenRun.nextDouble() * spikeDelaySteps)] += facilitationFac;

                }
                for (int k = 0; k < numInhPostFB[i]; k++) {
                    /* spike on kth postsynaptic neuron */
                    int target = inhPostFB[i][k];                    
                    spikeIEFB(target, (1.0-facIOverFacE) + facilitationFac*facIOverFacE);
                }
                
                numWaitingFBSpikes[i]--;
                nextFBSpikeTime[i] = t + fbSpreadLow + (fbSpreadHigh-fbSpreadLow) * rndGenRun.nextDouble();
            }
        }

        /* update the clock */
        fbClock += dt;
        noiseFBClock += dt;

        if (fbClock > fbRandEnd) {
            /* set a new random end time */
            fbRandEnd = fbEnd - rndGenRun.nextDouble() * 500.0;

            /* set a new random start time */
            fbRandStart = fbStart + rndGenRun.nextDouble() * 500.0;

            /* seek the first spike after the random feedback start time */
            int left = 0;
            int right= numFBSpikes - 1;
            int mid  = (left + right)/2;
            if (fbSpikeTimes[left] > fbRandStart){
                currentFBSpike = left;
            } else {
                while (fbSpikeTimes[mid] < fbRandStart || fbSpikeTimes[mid-1] >= fbRandStart){
                    if (fbSpikeTimes[mid] < fbRandStart){
                        left = mid;
                        if (left + right % 2 == 0){
                            mid = (left + right) / 2;
                        } else {
                            mid = (left + right) / 2 + 1;
                        }
                    } else {
                        right = mid;
                        mid = (left+right)/2;
                    }
                }
                currentFBSpike = mid;
            }

            /* reset the clock */
            fbClock = fbRandStart;
        }
    }
    
    
    public ArrayList<Integer> selectRegion(double angle, double angularWidth, double minDistFromCenter, double minDistFromBoundary){
        double angle1 = angle - angularWidth/2.0+PI/2;
        double angle2 = angle + angularWidth/2.0-PI/2;
        double[] vec1 = { cos(angle1), sin(angle1) };
        double[] vec2 = { cos(angle2), sin(angle2) };
        // reasoning for vec1 and vec2: a vector will only be within the angular width of angle if the inner product with vec1 and vec2 are both positive
        
        // initialize the list of selected neurons
        ArrayList<Integer> selected = new ArrayList<>(numFB);
        
        /* go through L6 cells */
        for (int i=0;i<numFB;i++){
            /* test if in region */
            boolean inRegion = false;
                        
            double x = locFB[i][0];
            double y = locFB[i][1];
            MyPoint hc1Loc = NetHelper.movLocToHC1(x, y, HCWidth);
            double dx = hc1Loc.x-HCWidth/2.0;
            double dy = hc1Loc.y-HCWidth/2.0;
            double distFromCenter = sqrt(dx*dx+dy*dy);
            if (dx*vec1[0] + dy*vec1[1] > 0
              & dx*vec2[0] + dy*vec2[1] > 0
              & distFromCenter > minDistFromCenter
              & x > minDistFromBoundary
              & y > minDistFromBoundary
              & x < 1500.0 - minDistFromBoundary
              & y < 1500.0 - minDistFromBoundary){
                inRegion = true;
            }
                
            /* done testing */
            
            /* if in, add to list */
            if (inRegion){
                selected.add(i);
            }            
        }
        
        return selected;
    }
    
    public void createL4L6Maps(){
        /*== initialize map ==*/
        l4L6Map = new int[numExcitatory][120];
        numL4L6Map = new int[numExcitatory];
        
        l6L4Map = new int[numFB][120]; //reverse map
        numL6L4Map = new int[numFB];                
        
        double boxHalfWidth=l4L6MapWidth/2.0;
        /*== iterate through all l6 cells (THIS CHOOSES SIMPLE 4Calpha CELLS FOR THE OUTER BOX) ==*/
        System.out.println("Making noise fb map");
        for (int i=0;i<numFB;i++){
            /*== find all l4 cells within boxHalfWidth radius of l6 cell ==*/
            double HCCenterX = floor(locFB[i][0] / HCWidth) * HCWidth + HCWidth/2.0;
            double HCCenterY = floor(locFB[i][1] / HCWidth) * HCWidth + HCWidth/2.0;
            double fx = HCCenterX + l4L6Scaling * (locFB[i][0]-HCCenterX);
            double fy = HCCenterY + l4L6Scaling * (locFB[i][1]-HCCenterY);
                
            for (int j=0;j<numExcitatory;j++){
                double x = locExc[j].x;
                double y = locExc[j].y;
                double dx = x-fx;
                double dy = y-fy;
                double dist = Math.sqrt(dx*dx + dy*dy);
                if (  dist < boxHalfWidth  ) {
                    int l4Cell = j;
                    int l6Cell = i;
                    l4L6Map[l4Cell][numL4L6Map[l4Cell]] = l6Cell;
                    numL4L6Map[l4Cell]++;
                    l6L4Map[l6Cell][numL6L4Map[l6Cell]] = l4Cell;
                    numL6L4Map[l6Cell]++;
                }
            }
        } 

        System.out.println("Done");
    }
    
    public void setupRegionList(){
        regionList = new ArrayList<>(pastAvgAngles.length-1);
        
        for (int setInd=0; setInd<pastAvgAngles.length-1; setInd++){            
            double angle1 = pastAvgAngles[setInd];
            double angle2 = pastAvgAngles[setInd+1];
            double setAngle = 0.5 * (angle1+angle2);
            double angleSweep = angle2-angle1;
            
            ArrayList<Integer> selection = this.selectRegion(setAngle, angleSweep, minDistFromCenter, minDistFromBoundary);
            regionList.add(selection);
        }
    }
    
    public void makePastAvg(){
        pastAvgCoupling = new int[numFB][MAXPASTAVGCOUPLING];
        numPastAvgCoupling = new int[numFB];        
        
        for (int setInd=0; setInd<pastAvgAngles.length-1; setInd++){            
            double angle1 = pastAvgAngles[setInd];
            double angle2 = pastAvgAngles[setInd+1];
            double setAngle = 0.5 * (angle1+angle2);
            double angleSweep = angle2-angle1;
            
            ArrayList<Integer> selectionPre = this.selectRegion(setAngle, angleSweep, 0.0, 0.0);
            ArrayList<Integer> selectionPost= this.selectRegion(setAngle, angleSweep, minDistFromCenter, minDistFromBoundary);
            
            for (int i=0;i<selectionPre.size();i++){
                for (int j=0;j<selectionPost.size();j++){
                    int l6Pre = selectionPre.get(i);
                    int l6Post= selectionPost.get(j);
                    pastAvgCoupling[l6Pre][numPastAvgCoupling[l6Pre]]=l6Post;
                    numPastAvgCoupling[l6Pre]++;
                }
            }
        }
    }
    
    public void makePresentAvg(){
        presentAvgCoupling = new int[numFB][30];
        numPresentAvgCoupling = new int[numFB];
        double presentAvgWidth = 150.0;
        
        for (int l6Pre=0;l6Pre < numFB;l6Pre++){
            double x1 = locFB[l6Pre][0];
            double y1 = locFB[l6Pre][1];
            for (int l6Post=0;l6Post<numFB;l6Post++){
                double x2 = locFB[l6Post][0];
                double y2 = locFB[l6Post][1];
                if (abs(x1-x2)<presentAvgWidth/2.0
                  & abs(y1-y2)<presentAvgWidth/2.0){
                    presentAvgCoupling[l6Pre][numPresentAvgCoupling[l6Pre]]=l6Post;
                    numPresentAvgCoupling[l6Pre]++;
                }
            }
        }
    }
    
    public void makeLongRangeCoupling(){
        longRangeCoupling = new int[numFB][12];
        
        System.out.println("Making L6 long range coupling...");
        
        double[] localDXs = { 50.0, 0.0, -50.0, 0.0 };
        double[] localDYs = { 0.0, 50.0, 0.0, -50.0};
        
        for (int l6Cell=0;l6Cell<numFB;l6Cell++){
            double x = locFB[l6Cell][0];
            double y = locFB[l6Cell][1];
            int HCm = (int)floor( x / HCWidth );
            int HCn = (int)floor( y / HCWidth );
            if (HCm >= 3 || HCn >= 3 || HCm < 0 || HCn < 0){
                for (int i=0;i<8;i++){
                    longRangeCoupling[l6Cell][i] = l6Cell;
                }
                continue;
            }
            double HCCenterX = (HCm + 0.5) * HCWidth;
            double HCCenterY = (HCn + 0.5) * HCWidth;
            double HCX = x - HCCenterX;
            double HCY = y - HCCenterY;
            
            int numLongRangeCoupling = 0;
            for (int HCm2 = 0; HCm2 < 3; HCm2++){
                for (int HCn2 = 0; HCn2 < 3; HCn2++){
                    if (HCm2==HCm && HCn2==HCn){
                        continue;
                    }
                    
                    double ghostHCCenterX = ((double)HCm2 + 0.5)*HCWidth;
                    double ghostHCCenterY = ((double)HCn2 + 0.5)*HCWidth;
                    double xFlips = 1.0 - 2.0 * (((double)(HCm2-HCm) + 2.0) % 2.0);
                    double yFlips = 1.0 - 2.0 * (((double)(HCn2-HCn) + 2.0) % 2.0);
                    double ghostX = ghostHCCenterX + HCX * xFlips;
                    double ghostY = ghostHCCenterY + HCY * yFlips;
                    
                    int l6Far = 0;
                    double l6FarDist = 100000;
                    for (int l6FarTest=0;l6FarTest<numFB;l6FarTest++){
                        double dx = ghostX - locFB[l6FarTest][0];
                        double dy = ghostY - locFB[l6FarTest][1];
                        double distSqr = dx*dx+dy*dy;
                        if (distSqr < l6FarDist){
                            l6Far=l6FarTest;
                            l6FarDist = distSqr;
                        }
                    }
                    
                    longRangeCoupling[l6Cell][numLongRangeCoupling] = l6Far;
                    numLongRangeCoupling++;
                }
            }            
            for (int anglei=0; anglei<4; anglei++){
                double ghostX = x + localDXs[anglei];
                double ghostY = y + localDYs[anglei];
                
                int l6Far = 0;
                double l6FarDist = 100000;
                for (int l6FarTest=0;l6FarTest<numFB;l6FarTest++){
                    double dx = ghostX - locFB[l6FarTest][0];
                    double dy = ghostY - locFB[l6FarTest][1];
                    double distSqr = dx*dx+dy*dy;
                    if (distSqr < l6FarDist){
                        l6Far=l6FarTest;
                        l6FarDist = distSqr;
                    }
                }
                
                longRangeCoupling[l6Cell][numLongRangeCoupling] = l6Far;
                numLongRangeCoupling++;
            }
        }
    }

    /* utility function which uniformly chooses k objects out of n */
    public int[] chooseKofN(int k, int n){
        int[] chosen = new int[k];

        int numChosen = 0;

        int indexToTest = 0;
        while (numChosen < k){
            // test index to test to see if it's in
            double probItsChosen = (double)(k-numChosen) / (double)(n-indexToTest);
            if (rndGenRun.nextDouble() < probItsChosen){
                chosen[numChosen] = indexToTest;
                numChosen++;
            }
            indexToTest++;
        }

        return chosen;
    }

    
    
    /* implement the ambient drive for this timestep */
    public void enactAmbientDrive() {
        for (int i = 0; i < numExcitatory; i++) {
            //ambient drive part here
            double numSpikes = (double) poisson.nextInt(dt * rateEAmb);
            if (numSpikes > 0) {
                synapseDifOfExponentials(tauAmbToEa, tauAmbToEb, gAmbExc[i], numSpikes, sAmb, dt);
            }
        }
        for (int i = 0; i < numInhibitory; i++) {
            //ambient drive part here
            double numSpikes = (double) poisson.nextInt(dt * rateIAmb);
            if (numSpikes > 0) {
                synapseDifOfExponentials(tauAmbToIa, tauAmbToIb, gAmbInh[i], numSpikes, sAmb, dt);
            }

        }
    }

    /* find any cortical spikes this timestep and implement them */
    public void enactCorticalSpikes() {
        /* prepare the stats vector */
        numSpikesNowExc = 0;
        numSpikesNowInh = 0;

        delayCylinderStep();

        /* detect E population spikes */
        for (int i = 0; i < numExcitatory; i++) {
            if (vExc[i] >= threshE[i] + min(eFatCap,fatigueExc[i]) ) {
                /* stats collecting */
                eSpikes++;
                registerESpike(i);
                
                /* tell the patch as well */
                if (membershipExc[i] == 1) {
                    patchOneFRExcTemp++;
                } else if (membershipExc[i] == 2) {
                    patchTwoFRExcTemp++;
                }
                
                /* add to fmax of all associated l6 cells */
                for (int j=0;j<numL4L6Map[i];j++){
                    int l6Cell = l4L6Map[i][j];
                    for (int k=0;k<whichAvg;k++){
                        double weight = timeAvgWeights[whichAvg-k-1];//n*(1.0 - d*d) / (n - (1.0/3.0)/(n*n) * (n+1.0) * ( (n+1.0)*(n-0.5) + 0.5)); 
                        timeAvgs[l6Cell][k] += weight / (timeAvgDt*numL6L4Map[l6Cell]);
                    }
                    for (int k=whichAvg+1;k<numAvgs+1;k++){
                        double weight = timeAvgWeights[numAvgs+whichAvg-k];
                        timeAvgs[l6Cell][k] += weight / (timeAvgDt*numL6L4Map[l6Cell]);
                    }
                }

                delaySpikeE(i);
                
                /* tell E fatigue mechanism */
                fatigueExc[i] += eFatAmt / exp(-refractoryE / eFatSigma);
            }
        }
        /* among I population */
        for (int i = 0; i < numInhibitory; i++) {
            if (vInh[i] >= threshI[i]) {
                /* stats collecting */
                iSpikes++;
                registerISpike(i);

                /* tell the patch as well */
                if (membershipInh[i] == 1) {
                    patchOneFRInhTemp++;
                } else if (membershipInh[i] == 2) {
                    patchTwoFRInhTemp++;
                }

                fastSpikeI(i);
                
                /* record spike for depression mechanism */
                /* update I neuron spike record */
                last3SpikesInd[i]++;
                last3SpikesInd[i]=last3SpikesInd[i] % 3;
                last3Spikes[i][last3SpikesInd[i]] = t;
            }
        }
    }

    /* slightly faster version of above */
    public void enactCorticalSpikesNOSTAT() {
        delayCylinderStep();

        /* detect E population spikes */
        for (int i = 0; i < numExcitatory; i++) {
            if (vExc[i] >= threshE[i]) {
                /* notify feedback */
                //for(int j=0;j<50;j++) timeLocalFR[j]++; //so timeLocalFR is measured in spikes / 50ms

                delaySpikeE(i);
            }
        }
        /* among I population */
        for (int i = 0; i < numInhibitory; i++) {
            if (vInh[i] >= threshI[i]) {
                fastSpikeI(i);
            }
        }
    }

    public void delayCylinderStep() {
        spikeDelayIndex++;
        spikeDelayIndex %= spikeDelaySteps;
        for (int i = 0; i < numExcitatory; i++) {
            if (delayCylinder[i][spikeDelayIndex] > 0) {
                spikeEE(i, delayCylinder[i][spikeDelayIndex]);
                delayCylinder[i][spikeDelayIndex] = 0;
            }
        }

        for (int i = 0; i < numExcitatory; i++) {
            if (delayCylinderFB[i][spikeDelayIndex] > 0) {
                spikeEEFB(i, delayCylinderFB[i][spikeDelayIndex]);
                delayCylinderFB[i][spikeDelayIndex] = 0;
            }
        }
    }

    /* make E cell send delayed spikes */
    public void delaySpikeE(int index) {
        /* begin refractory period */
        vExc[index] = 0;
        sleepExc[index] = -(refractoryE / dt);
        /* enact spike */
        /*===For E->E, put on all target cell delay rings===*/
        for (int j = 0; j < numEPostExc[index]; j++) {
            int target = cEExc[index][j];
            delayCylinder[target][(int) floor((spikeDelayIndex + minSpikeDelaySteps + rndGenRun.nextDouble() * spikeDelayActualSteps) % spikeDelaySteps)] += modifier(sModSeedExc[index], sModSeedExc[target]);
        }
        for (int j = 0; j < numIPostExc[index]; j++) {
            int target = cIExc[index][j];
            spikeIE(target, 1.0);
        }
    }

    /* make I cell send spikes */
    public void fastSpikeI(int index) {
        /* begin refractory period */
        vInh[index] = 0;
        sleepInh[index] = -(refractoryI / dt);
        /* enact spike to E's */
        for (int j = 0; j < numEPostInh[index]; j++) {
            int target = cEInh[index][j];
            fastSpikeEI(target, index);
        }

        /* enact spike to I's */
        for (int j = 0; j < numIPostInh[index]; j++) {
            int target = cIInh[index][j];
            fastSpikeII(target);
        }
    }

    public void spikeEE(int target, double modifier) {
        //double f0 = getF0();        
        //double strength = sEE * (sEELowMultiplier + (sEEHighMultiplier-sEELowMultiplier) * rndGenRun.nextDouble()) / (1.0-f0*tauInitialAMPADamper);
        double strength = sEE * (sEELowMultiplier + (sEEHighMultiplier-sEELowMultiplier) * rndGenRun.nextDouble()); 
        synapseDifOfExponentials(tauAMPAToEa, tauAMPAToEb, gAMPAExc[target], 1.0-fracNMDAEE, strength * modifier, eOffset);
        synapseNMDA(gNMDAExc[target], fracNMDAEE, strength * modifier);
        
        //initialAMPADamperExc[target] += initialRise;
    }

    public void spikeIE(int target, double modifier) {
        //double f0 = getF0();
        //double initialRise = synapseAMPAGABA(tauAMPA, gAMPAInh[target], 2.0 / 3.0, sIE * modifier / (1.0-f0*tauInitialAMPADamper), eOffset);
        synapseDifOfExponentials(tauAMPAToIa, tauAMPAToIb, gAMPAInh[target], 1.0-fracNMDAIE, sIE * modifier, eOffset);
        synapseNMDA(gNMDAInh[target], fracNMDAIE, sIE * modifier);
        
        //initialAMPADamperInh[target] += initialRise;
    }

    public void spikeEEFB(int target, double modifier) {
        //double f0 = getF0();
        //double strength = (sEEL6Low + rndGenRun.nextDouble() * (sEEL6High-sEEL6Low)) / (1.0-f0*tauInitialAMPADamper);
        double strength = (sEEL6Low + rndGenRun.nextDouble() * (sEEL6High-sEEL6Low)) * sEgal[numInputsExc[target]];
        synapseDifOfExponentials(tauAMPAToEa, tauAMPAToEb, gAMPAFBExc[target], 1.0-fracNMDAEE, strength * modifier, eOffset);
        synapseNMDA(gNMDAFBExc[target], fracNMDAEE, strength * modifier);
        
        //initialAMPADamperExcFB[target] += initialRise;
    }

    public void spikeIEFB(int target, double modifier) {
        //double f0 = getF0();
        //double initialRise = synapseAMPAGABA(tauAMPA, gAMPAFBInh[target], 2.0 / 3.0, sIEL6 * modifier / (1.0-f0*tauInitialAMPADamper), eOffset);
        synapseDifOfExponentials(tauAMPAToIa, tauAMPAToIb, gAMPAFBInh[target], 1.0-fracNMDAIE, sIEL6 * modifier, eOffset);
        synapseNMDA(gNMDAFBInh[target], fracNMDAIE, sIEL6 * modifier);
        
        //initialAMPADamperInhFB[target] += initialRise;
    }

    public void fastSpikeEI(int target, int source) {
        double modifier = sEI * getSEIDepFac(source) * modifier(sModSeedExc[target], sModSeedInh[source]);

        synapseTwoDifOfExponentials(tauGABAa, tauGABAb, tauGABAc, tauGABAd, weightToSecondGABA, gGABAExc[target], 1.0, modifier, dt);
    }

    public void fastSpikeII(int target) {
        double modifier = (sIIOverSEILow + sIIOverSEIRange * rndGenRun.nextDouble()) * sII;

        synapseTwoDifOfExponentials(tauGABAa, tauGABAb, tauGABAc, tauGABAd, weightToSecondGABA, gGABAInh[target], 1.0, modifier, dt);
    }
    
    public double getSEIDepFac(int source){        
        int ind = last3SpikesInd[source];
        double lastSpikeTime = last3Spikes[source][ind];
        return 1.0-sEIDepAmt*exp(-(t-lastSpikeTime)/sEIDepTimeConstant);
    }
    
    public void evolveEFatigue(){
        for (int i=0;i<numExcitatory;i++){
            fatigueExc[i] *= exp(-dt/eFatSigma);
        }
    }
    
    public double getL6FRBoth(double fr){
        //[0 8 13.5 18 1000]
        //[5.5 5.5 25.1 48 48+(67.5-48)/(30-18)*(1000-18)]
        //[(5.5-5.5)/(8.0-0.0) (25.1-5.5)/(13.5-8.0) (48.0-25.1)/(18.0-13.5) (67.5-48.0)/(30.0-18.0)]
        int i = 1;
        while (l6L4xs[i] < fr){
            i++;
        }
        i--;
        return max( l6L4ys[i] + l6L4slopes[i] * (fr-l6L4xs[i]), minL6 );
    }
    
    public double getL6L4FRRatioPoisson(double fr){
        //Find x1, x2
        int indr=0;
        for (int i=1;i<fbPoissonxs.length;i++){
            if (fr < fbPoissonxs[i]){
                indr=i;
                break;
            }
        }
        int indl=indr-1;
        return fbPoissonys[indl] + fbPoissonSlopes[indl] * (fr - fbPoissonxs[indl]);
    }
    
    public void setFBPoisson(double[] fbPoissonxs, double[] fbPoissonys){
        this.fbPoissonxs=fbPoissonxs;
        this.fbPoissonys=fbPoissonys;
        for (int i=0;i<fbPoissonxs.length-1;i++){
            fbPoissonSlopes[i] = (fbPoissonys[i+1]-fbPoissonys[i]) / (fbPoissonxs[i+1]-fbPoissonxs[i]);
        }
        fbPoissonSlopes[fbPoissonxs.length-1] = 0.0;
    }
    
    public void enactMinL6(){
        double ms = floor(t);
        if ((int)ms % minL6RefreshPeriod == 0 & t-ms < dt){
            
            double maxAvg = 0;
            
            for (int i=0;i<numFB;i++){
                //double avg = fbPComplex[i] * exp( (timeLastSpike[i] - t) / sigmaFB );
                double avg = 2.0 * timeAvgs[i][whichAvg];
                for (int j=0;j<12;j++){
                    int l6Cell = longRangeCoupling[i][j];
                    //avg += fbPComplex[l6Cell] * exp( (timeLastSpike[l6Cell] - t) / sigmaFB );
                    avg += timeAvgs[l6Cell][whichAvg];
                }
                avg /= 14.0;
                if (avg > maxAvg){
                    maxAvg = avg;
                }

            }
            minL6 = max(peakLGNMinL6, minL6MaxL6Ratio * this.getL6FRBoth(maxAvg));
        }
    }
    
    public void setPeakLGNMinL6(){
        double s =  contrast * spFreqMod;
        peakLGNMinL6 = s * largestMinL6 + (1-s) * fbPoissonxs[0];
    }
    
    public void initializeTimeAvgs(){
        timeAvgs = new double[numFB][numAvgs+1];
    }
    
    public void enactTimeAvgSwitch(){
        if (t > nextSwitchTime){
            // then switch
            for (int i=0;i<numFB;i++){
                timeAvgs[i][whichAvg] = 0.0;
            }
            whichAvg=(whichAvg+1) % (numAvgs+1);
            nextSwitchTime = t + (timeAvgDt*1000.0) / (double)numAvgs;
        }
    }

    public void resetStats() {
        eSpikes = 0;
        iSpikes = 0;
    }

    public void resetNetState() {
        resetStats();

        for (int i = 0; i < numExcitatory; i++) {
            vExc[i] = rndGenRun.nextDouble() * .9;
            sleepExc[i] = 0;
        }
        for (int i = 0; i < numInhibitory; i++) {
            vInh[i] = rndGenRun.nextDouble() * .9;
            sleepInh[i] = 0;
        }
        gLGNExc = new double[numExcitatory][4];
        gAMPAExc = new double[numExcitatory][4];
        gGABAExc = new double[numExcitatory][5];
        gNMDAExc = new double[numExcitatory][3];
        gNMDAInh = new double[numInhibitory][3];
        gGABAInh = new double[numInhibitory][5];
        gLGNInh = new double[numInhibitory][4];
        gAMPAInh = new double[numInhibitory][4];

        /* reset e-delay cylinder */
        for (int i = 0; i < numExcitatory; i++) {
            for (int j = 0; j < spikeDelaySteps; j++) {
                delayCylinder[i][j] = 0.0;
            }
        }
        for (int i = 0; i < numExcitatory; i++) {
            for (int j = 0; j < spikeDelaySteps; j++) {
                delayCylinderFB[i][j] = 0.0;
            }
        }

        t = 0;
        numSteps = 0;
    }
    
    public void recordCurrent(){
        if (stepsToNextCurrentSample > 0){
            stepsToNextCurrentSample--;
        } else {
            stepsToNextCurrentSample = stepsBetweenCurrentSamples;
            
            //grab sample (2nd index: E FB LGN Amb I Leak)            
            for (int i=0;i<numExcitatory;i++){
                avgVExc[i] += vExc[i];
                double mE=(14/3-vExc[i]);
                double mI=(vExc[i]+2/3);
                currentToE[i][0] += mE*(gAMPAExc[i][0]+gNMDAExc[i][0]);
                currentToE[i][1] += mE*(gAMPAFBExc[i][0]+gNMDAFBExc[i][0]);
                currentToE[i][2] += mE*gLGNExc[i][0];
                currentToE[i][3] += mE*gAmbExc[i][0];
                currentToE[i][4] += mI*gGABAExc[i][0];
                currentToE[i][5] += vExc[i]*gLExc;
            }
            for (int i=0;i<numInhibitory;i++){
                avgVInh[i] += vInh[i];
                double mE=(14/3-vInh[i]);
                double mI=(vInh[i]+2/3);
                currentToI[i][0] += mE*(gAMPAInh[i][0]+gNMDAInh[i][0]);
                currentToI[i][1] += mE*(gAMPAFBInh[i][0]+gNMDAFBInh[i][0]);
                currentToI[i][2] += mE*gLGNInh[i][0];
                currentToI[i][3] += mE*gAmbInh[i][0];
                currentToI[i][4] += mI*gGABAInh[i][0];
                currentToI[i][5] += vInh[i]*gLInh;
            }
            numCurrentSamples++;
        }
    }
    
    public void resetRates(){
        /*rate mechanism*/
        Arrays.fill(preRateExc,0.0);
        Arrays.fill(preRateInh,0.0);
        Arrays.fill(preRateFB,0.0);
        rateStartTime=t;
        
        /*current record mechanism*/
        for (int i=0;i<numExcitatory;i++){
            Arrays.fill(currentToE[i],0.0);
        }
        for (int i=0;i<numInhibitory;i++){
            Arrays.fill(currentToI[i],0.0);
        }
        Arrays.fill(avgVExc,0.0);
        Arrays.fill(avgVInh,0.0);
        stepsToNextCurrentSample=stepsBetweenCurrentSamples;
        numCurrentSamples=0;
    }

    /*===Go forward in time one time step===*/
    public void step() {
        evolveVExc();
        evolveVInh();
        evolveConductances();
        evolveSleep();
        evolveLGN();
        evolveEFatigue();

        enactAmbientDrive();
        enactFB();
        enactCorticalSpikes();
        enactMinL6();
        enactTimeAvgSwitch();

        recordCurrent();
        
        /* update time */
        t += dt;
    }

    public void jump(int numTimeSteps) {
        for (int i = 0; i < numTimeSteps; i++) {
            step();
        }
    }

    /* the following is for saving and loading data to disk */
    public void saveFBConnectivity(String fname) {
        try {
            FileOutputStream fos = new FileOutputStream(linPATH.concat(fname));
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeInt(numFB);
            oos.writeObject(excPostFB);
            oos.writeObject(inhPostFB);
            oos.writeObject(numExcPostFB);
            oos.writeObject(numInhPostFB);
            oos.writeObject(locFB);
            oos.writeObject(simpleFB);
            oos.writeObject(phaseFB);
            oos.writeObject(highFiringComplexFB);
            oos.writeDouble(fracSimple);
            oos.writeDouble(simpleAvg);
            oos.writeDouble(simpleAmpLow);
            oos.writeDouble(simpleAmpHigh);
            oos.writeObject(simpleAmp);
            oos.writeDouble(deregAmpLow);
            oos.writeDouble(deregAmpHigh);
            oos.writeObject(deregAmp);
        } catch (Exception e) {
        }
    }

    public void loadFBConnectivity(String fname) {
        try {
            FileInputStream fis = new FileInputStream(linPATH.concat(fname));
            ObjectInputStream ois = new ObjectInputStream(fis);
            numFB = ois.readInt();
            excPostFB = (int[][]) ois.readObject();
            inhPostFB = (int[][]) ois.readObject();
            numExcPostFB = (int[]) ois.readObject();
            numInhPostFB = (int[]) ois.readObject();
            locFB = (double[][]) ois.readObject();
            simpleFB = (double[]) ois.readObject();
            phaseFB = (double[]) ois.readObject();
            highFiringComplexFB = (double[]) ois.readObject();
            fracSimple = ois.readDouble();
            simpleAvg = ois.readDouble();
            simpleAmpLow = ois.readDouble();
            simpleAmpHigh = ois.readDouble();
            simpleAmp = (double[]) ois.readObject();
            deregAmpLow = ois.readDouble();
            deregAmpHigh = ois.readDouble();
            deregAmp = (double[]) ois.readObject();
        } catch (Exception e) {
        }
    }

    public void saveConnectivity(String fname) {
        try {
            FileOutputStream fos = new FileOutputStream(linPATH.concat(fname));
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(numEPostExc);
            oos.writeObject(numIPostExc);
            oos.writeObject(numEPostInh);
            oos.writeObject(numIPostInh);
            oos.writeObject(numEPreExc);
            oos.writeObject(numIPreExc);
            oos.writeObject(numEPreInh);
            oos.writeObject(numIPreInh);
            oos.writeObject(cEExc);
            oos.writeObject(cIExc);
            oos.writeObject(cEInh);
            oos.writeObject(cIInh);
            oos.writeObject(cPreEExc);
            oos.writeObject(cPreIExc);
            oos.writeObject(cPreEInh);
            oos.writeObject(cPreIInh);

        } catch (Exception e) {
        }
    }

    public void loadConnectivity(String fname) {
        try {
            FileInputStream fis = new FileInputStream(linPATH.concat(fname));
            ObjectInputStream ois = new ObjectInputStream(fis);
            numEPostExc = (int[]) ois.readObject();
            numIPostExc = (int[]) ois.readObject();
            numEPostInh = (int[]) ois.readObject();
            numIPostInh = (int[]) ois.readObject();
            numEPreExc = (int[]) ois.readObject();
            numIPreExc = (int[]) ois.readObject();
            numEPreInh = (int[]) ois.readObject();
            numIPreInh = (int[]) ois.readObject();
            cEExc = (int[][]) ois.readObject();
            cIExc = (int[][]) ois.readObject();
            cEInh = (int[][]) ois.readObject();
            cIInh = (int[][]) ois.readObject();
            cPreEExc = (int[][]) ois.readObject();
            cPreIExc = (int[][]) ois.readObject();
            cPreEInh = (int[][]) ois.readObject();
            cPreIInh = (int[][]) ois.readObject();

        } catch (Exception e) {
        }
    }

    /* saves:
     E thresholds
     I thresholds
     synaptic strength mod seeds

     */
    public void saveMiscRandComponents(String fname) {
        try {
            FileOutputStream fos = new FileOutputStream(linPATH.concat(fname));
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(sModSeedLGNON);
            oos.writeObject(sModSeedLGNOFF);
            oos.writeObject(sModSeedExc);
            oos.writeObject(sModSeedInh);
            oos.writeObject(threshE);
            oos.writeObject(threshI);
            oos.writeObject(threshEBounds);
            oos.writeObject(threshIBounds);
        } catch (Exception e) {
        }
    }

    public void loadMiscRandComponents(String fname) {
        try {
            FileInputStream fis = new FileInputStream(linPATH.concat(fname));
            ObjectInputStream ois = new ObjectInputStream(fis);
            sModSeedLGNON = (double[]) ois.readObject();
            sModSeedLGNOFF = (double[]) ois.readObject();
            sModSeedExc = (double[]) ois.readObject();
            sModSeedInh = (double[]) ois.readObject();
            threshE = (double[]) ois.readObject();
            threshI = (double[]) ois.readObject();
            threshEBounds = (double[]) ois.readObject();
            threshIBounds = (double[]) ois.readObject();
        } catch (Exception e) {
        }
    }
}
