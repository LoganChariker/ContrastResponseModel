%===announce which script is running===%
if ~exist('scriptnum','var')
    scriptnum=1;
end
disp(['running script ' num2str(scriptnum) '...']);

%===select model version to run===%
projectname = 'SimJNSVer14-6';

%===durations===%
waitDurationInMs = 5000; %run the model for this duration without recording
runDurationInMs  = 5000; %record for this duration
recording = true; %setting to true causes spikes to be recorded
recordWhat = 'spikesAndPerCellAvgCurrent';
disp(['stat-collecting runtime = ' num2str(runDurationInMs/1000) ' virtual seconds...']);

%===Save file settings===%
savefolder = './runResults/';
filename = ['run' num2str(scriptnum)];
fullsavename = [savefolder filename];
disp(['savename = ' fullsavename]);

% %===make save folder if doesn't exist===%
if ~exist(savefolder)
    mkdir(savefolder);
end

%===set base parameters (some are no longer used)===%
 %===non-randomized, non-constructive parameters===%
prms = containers.Map();
prms('aChSELGNFac') = 1.05;%
prms('lGNThresh') = 1.2;%
prms('sLGNNoise') = 0.075;%
prms('tf') = 4; %Hz
prms('pValues') = [2.3, 1.5, 1.5, 1.15, 1.15, 1.0, 1.0, 1.15, 1.15, 1.5, 1.5, 2.3, 2.3];
prms('simpleAvg') = 9/23;
prms('rhoE') = { 0 };

 %===constructive variables===%
prms('randomPicksEEConnectivity') =.18;
prms('threshIBounds') = [.9 1.1];
prms('threshEBounds') = [.95 1.05];
prms('peakCEI') = .6;
prms('peakCIE') = .6;
prms('peakCII') = .6;
prms('sModSTD') = .1;
prms('simpleAmpLow') = 8/23; %12/23
prms('simpleAmpHigh')= 12/23;%18/23
prms('deregAmpLow')  = 0/20;
prms('deregAmpHigh') = 2/20;
prms('fbStdSqr') = 4.0;

 %===saved file options===%
prms('load4CConnections') = false;
prms('l4CConnectionsFile') = '';
prms('loadFBConnections') = false;
prms('fBConnectionsFile') = '';
prms('loadMiscRandComponents') = false;
prms('miscRandComponentsFile') = '';
prms('fbStart') = 1450;
prms('fbEnd') = 2450;
prms('fbClock') = 1450;
prms('currentFBSpikes') = 0;

%===layout info===%
prms('numHCRows')=3;
prms('numExcitatoryPerHC')=3000;
prms('numInhibitoryPerHC')=1000;
prms('numLGNBlockRows')=6;
prms('numLGNPerRowPerBlock')=2;
prms('numLGNPerColPerBlock')=2.5;
%===layout info(end)===%
%===base parameters set===%

%===create parameter manager and set base parameters===%
disp('setting up parameter manager...')
PM = ParameterManager();
PM.setBasePrms(prms);

%===add parameter sets to run===%
%allsets
p = containers.Map();

%set1 sEE ranged
boost = [ 4.5-1.8 ];
p_0 = 0.45;
p_6 = 0.7;
p('fracSimple') = { 0.15 };

p('sEELowMultiplier') = { 0.9 };
p('sEEHighMultiplier') = { 1.08 };
p('fbWideNumSDsE') = { 3.0 };
p('fbWideNumSDsI') = { 3.0 };
p('fbNarrowNumSDsE') = { 1.5 };
p('fbNarrowNumSDsI') = { 1.5 };
p('LGNInputFile') = { 'LGNToCortexAug29.mat'... 
                    };

p('numPreDistCutoff') = { 0.05 };
p('numPreDistCutoffFB') = { 0.10 };
p('sEICpxRatio') = { 1 };

% ====================new strength egal code=================%
p('sEgal') = {  [1.03  1.027  1.024  1 1 1 1 1 1 1 1] };

% =====================new fb code isolated here==============%
p('fbSpreadLow') =   { 2 };
p('fbSpreadHigh') =  { 6 };
p('fbSpikeSpreadProbChanges') = { [1 2 3 1000] };
p('fbSpikeSpreadProbs') = { [1 2/3 1/3 0] };

% fb curve
fbcurvex1=[0	 	4.0485	 	4.742	 	5.1083	 	5.5826	 	6.0018	 	6.5151	 	7.0053	 	7.5329	 	7.9947	 	8.492	 	9.0107	 	9.9982	 	10.7928	 	11.3864	 	12.19           13.0515         14.6273         1044.43];
fbcurvey1=[4	 	4           6.0497	 	7.4292	 	9.3914	 	11.8048	 	18.468	 	25.7083	 	33.7356	 	42.0252	 	51.0493	 	61.5425	 	79.9056	 	91.1938	 	97.4283	 	103.348         107.5044	 	112.6054	 	2553.1697];
fbcurvex2=[0	 	4.0485	 	4.742	 	5.1083	 	5.5826	 	6.0018	 	6.5151	 	7.0053	 	7.5329	 	7.9947	 	8.492	 	9.0107	 	9.9982	 	10.7928	 	11.4205	 	12.2717         13.3983         14.7565         22.6073         1044.43];
fbcurvey2=[4	 	4           6.0497	 	7.4292	 	9.3914	 	11.8048	 	18.468	 	25.7083	 	33.7356	 	42.0252	 	51.0493	 	61.5425	 	79.9056	 	91.1938	 	97.8195	 	103.6134	 	108.0623	 	111.0627	 	121.4727	 	2553.1697];

l6L4xs = { fbcurvex2 ...
         };
l6L4ys = { fbcurvey2 ...
         };
     
for curveInd = 1 : length(l6L4xs)
p('curveInd') = {curveInd};
[xs, ys ,sls] = getL6PiecewiseParams(l6L4xs{curveInd},l6L4ys{curveInd});
p('l6L4xs') = { xs };
p('l6L4ys') = { ys };
p('l6L4slopes') = { sls };
p('fbPoissonxs') = { [0 5 1000] };
p('fbPoissonys') = { [1 1 1]*ys(1) };
p('minL6MaxL6Ratio') = { .22 };
p('xL4Start') = { 6.5 };
p('xL4End')   = { 10.5 };
p('fbBGRate') = { 4.5 };

% fac curve
faccurvex=[0	 	29.8	 	37.48	 	40.739	 	43.9809	 	46.0072	 	48.0095	 	50.0119	 	55.0179	 	1142.7];
faccurvey=[1	 	1           1.0182	 	1.0297	 	1.054	 	1.0732	 	1.0856	 	1.0903	 	1.0899	 	1.3045];

facxs = { faccurvex...
        };
facys = { faccurvey...
        };
for faccurveInd = 1 : length(facxs)
p('faccurveInd') = { faccurveInd };
[xs, ys, sls] = getL6PiecewiseParams(facxs{faccurveInd},facys{faccurveInd});
p('facxs') = { xs };
p('facys') = { ys };
p('facslopes') = { sls };
p('facIOverFacE') = { 1 };


% l4 to l6 coupling code
p('l4L6MapWidth') = {75};

% =====================old fb code isolated here==============%
p('fbWideSDE') = { 110  };
p('fbWideSDI') = { 110  };
p('fbNarrowSDE') = { 110  };
p('fbNarrowSDI') = { 110  };
p('fbInputFile') = { 'fbInputFileOct16-Plus00PoisExt2.mat' };
p('threshType') = { 2 };
p('seedL4') = { 1 };
p('seedL6') = { 2 };
p('seedMisc') = { 3 };
p('seedRun') = { 4 };
p('seedEngine') = { 5 };     
p('fracNMDAEE') = { .20 };
p('fracNMDAIE') = { .20 };
p('sEIDepPerSpike') = { [0 1 1.5 1.7] };
p('sEIDepAmt') = { .12 };

sEIs = { 2.07 };
for sEIInd=1:length(sEIs)
p('sEIOverSEE') = {sEIs{sEIInd}};
p('sIIOverSEIBounds') = { [.65 .85]/.75*1.3965/sEIs{sEIInd} };
p('sEIDepTimeConstant') = { 20 };
p('sEIDepLookback') = { [10 19 24] };
p('simpleSpikeProb') = { 0.625/1000.0 * 0.1 };
p('highCpxProb') = { 12.0/1000.0 * 0.1 };
p('lowCpxProb') = { 0.625/1000.0 * 0.1 };
p('spikeDelaySteps') = { 15 };
p('weightToSecondGABA') = { 0 };
p('sIEL6OverSEEL6FracL4') = { .97 };

tauGABAas    = { 5      };
tauGABAbs    = { .5     };
tauGABAcs    = { 3      };
tauGABAds    = { .5     };
tauAMPAToEas = { 2      };
tauAMPAToIas = { 2      };
tauAMPAToEbs = { .6    };
tauAMPAToIbs = { .5     };
sEIDepDeltas = { .051   };

for tauInd=1:length(tauGABAas)
p('tauInd') = { tauInd };
p('tauGABAa') = { tauGABAas{ tauInd } };
p('tauGABAb') = { tauGABAbs{ tauInd } };
p('tauGABAc') = { tauGABAcs{ tauInd } };
p('tauGABAd') = { tauGABAds{ tauInd } };
p('tauAMPAToEa') = { tauAMPAToEas{ tauInd } };
p('tauAMPAToIa') = { tauAMPAToIas{ tauInd } };
p('tauAMPAToEb') = { tauAMPAToEbs{ tauInd } };
p('tauAMPAToIb') = { tauAMPAToIbs{ tauInd } };
p('tauAmbToEa') = { 10 };
p('tauAmbToEb') = { .5 };
p('tauAmbToIa') = { 10 };
p('tauAmbToIb') = { .5 };
p('sEIDepDelta')    = { .051 };

sELGNOverSEEs = { 2.2    };
sILGNOverSEEs = { 3.5    };
fbEgals = { [70 64 58 51 47.5 43.5 39]-1   ...
          };
for setnum = 1:length(sELGNOverSEEs)
p('setnum') = { setnum };
p('sELGNOverSEE') = { sELGNOverSEEs{setnum} };
p('sILGNOverSEE') = { sILGNOverSEEs{setnum} };

fbWideEgal =   [(1-p_0)*fbEgals{setnum}(1:3) (1-p_6)*fbEgals{setnum}(4:end)];
fbNarrowEgal = [    p_0*fbEgals{setnum}(1:3)     p_6*fbEgals{setnum}(4:end)];
p('fbWideEgal')   = { fbWideEgal };
p('fbNarrowEgal') = { fbNarrowEgal };

sEEs = { .023 };
sEEL6Lows = { .35 * .023*2/3 };
sEEL6Highs= { .35 * .023*4/3 };
for sEEInd = 1 : length(sEEs)
p('sEEInd') = { sEEInd };
p('sEE') = { sEEs{sEEInd} };
p('sEEL6Low') = { sEEL6Lows{sEEInd} };
p('sEEL6High') = { sEEL6Highs{sEEInd} };
p('rateEAmb') = { .62 };
p('spikeDelayActualSteps') = { 10 };
p('minSpikeDelaySteps') = { 5 };

% =====================new e fatigue code=====================%
p('eFatAmt') = { .2  };
p('eFatSigma') = { 10/log(10) };
p('eFatCap') = { 100 };
% ============================================================%

egals = {   [18.65    18.40    17.70    16.5 15.1 13.800 12.800]/100 ...             
        }; 
for egalInd = 1:length(egals)
p('egalInd') = { egalInd };
p('egal') = { egals{egalInd} };

p('angle') = { 0*pi/8 };
p('sf') = { 2.5 };
p('contrast') = { 1 };
p('tkHalfHeight') = { 100 };%{ 40 };
p('l4L6CpxToRemove') = { 1 };
p('numPreToRemove') = { 4 };
p('sIEOverSEE') = { .256 };
p('rateIAmb')   = { .50 };

fbIExpNumPres = { 104.6 };
for fbIExpNumPreInd=1:length(fbIExpNumPres)
fbIExpNumPre=fbIExpNumPres{fbIExpNumPreInd};
q_0 = 2/3;
p('fbIExpNumPreWide')   = { fbIExpNumPre*(1-q_0) };
p('fbIExpNumPreNarrow') = { fbIExpNumPre*q_0 };
p('EToIBoost') = { 10 };
p('135EToIDec') = { 12 };
p('135L6ToEBoost') = { .5 };
p('largestMinL6') = { 17 };
PM.addAllCombinations(p);
end
end
end
end
end
end
end
end

disp(['Number of runs = ' num2str(length(PM.prmList))])

%%
%===save parameter manager in run folder===%
if scriptnum==1 %~exist(PMFile)
    disp('saving parameters...');
    PMFile = [savefolder 'PM.mat'];
    try
        save([savefolder 'PM.mat'], 'PM');
        disp('saved.')
    catch err
        disp('could not save.');
    end    
end

%===create java object===%
disp('creating Network object...')
disp('loading jar file...')
javaaddpath( [pwd '/' projectname '.jar'] );
a = sim.Network();

%===update path variables===%
archstr = computer('arch');
if strcmp('win64',archstr)
    a.winPATH = [pwd '\'];        
elseif strcmp('glnxa64',archstr)
    a.linPATH = [pwd '/'];
elseif strcmp('maci64',archstr)
    a.macPATH = [pwd '/'];
else
    error('Error: architecture not recognized.  Cannot set paths.');
end

%===fill in java object parameters===%
PM.copyParameterSetToJava(scriptnum, a);
%===================fill in new fb code======================%
% ====================new strength egal code=================%
a.sEgal = PM.prmList{scriptnum}.parameters('sEgal');

% =====================new e fatigue code=====================%
a.eFatAmt = PM.prmList{scriptnum}.parameters('eFatAmt');        
a.eFatSigma = PM.prmList{scriptnum}.parameters('eFatSigma');
a.eFatCap = PM.prmList{scriptnum}.parameters('eFatCap');

% =====================new fb code isolated here==============%
a.fbSpreadLow = PM.prmList{scriptnum}.parameters('fbSpreadLow');        
a.fbSpreadHigh = PM.prmList{scriptnum}.parameters('fbSpreadHigh');        
a.fbSpikeSpreadProbChanges = PM.prmList{scriptnum}.parameters('fbSpikeSpreadProbChanges');        
a.fbSpikeSpreadProbs = PM.prmList{scriptnum}.parameters('fbSpikeSpreadProbs');        

% fb curve
a.l6L4xs = PM.prmList{scriptnum}.parameters('l6L4xs');        
a.l6L4ys = PM.prmList{scriptnum}.parameters('l6L4ys');        
a.l6L4slopes = PM.prmList{scriptnum}.parameters('l6L4slopes');
a.setFBPoisson( PM.prmList{scriptnum}.parameters('fbPoissonxs'), ...
                PM.prmList{scriptnum}.parameters('fbPoissonys') );
a.minL6MaxL6Ratio = PM.prmList{scriptnum}.parameters('minL6MaxL6Ratio');        
a.xL4Start = PM.prmList{scriptnum}.parameters('xL4Start');        
a.xL4End = PM.prmList{scriptnum}.parameters('xL4End');        
a.fbBGRate = PM.prmList{scriptnum}.parameters('fbBGRate');        

% fac curve
a.facxs = PM.prmList{scriptnum}.parameters('facxs');        
a.facys = PM.prmList{scriptnum}.parameters('facys');        
a.facslopes = PM.prmList{scriptnum}.parameters('facslopes');        
a.facIOverFacE = PM.prmList{scriptnum}.parameters('facIOverFacE');        

% l4 to l6 coupling code
a.l4L6MapWidth = PM.prmList{scriptnum}.parameters('l4L6MapWidth');        

% time averaging
tkHalfHeight = PM.prmList{scriptnum}.parameters('tkHalfHeight');
f = @(x) max(0,min(1,1./(1+exp((x-tkHalfHeight)/7))));
ts = 2.5:2.5:100;
timeAvgWeights = f(ts);
timeAvgWeights = timeAvgWeights / ( sum(timeAvgWeights) / length(timeAvgWeights) );
a.numAvgs = length(timeAvgWeights);
a.timeAvgDt = ts(end)/1000;
a.timeAvgWeights = timeAvgWeights;
a.whichAvg = 0;
a.nextSwitchTime = a.t + (a.timeAvgDt*1000.0) / a.numAvgs;
a.initializeTimeAvgs
% =====================old fb code isolated here==============%
% =====================new sei depression=====================%
sEIDepAmt = PM.prmList{scriptnum}.parameters('sEIDepAmt');
sEIDepTimeConstant = PM.prmList{scriptnum}.parameters('sEIDepTimeConstant');
a.sEIDepAmt = sEIDepAmt;
a.sEIDepTimeConstant = sEIDepTimeConstant;

% =====================RNG SEED CHANGES HERE (FOR TESTING--DELETE THIS)===%
seedL4 = 10;
a.rndGenL4 = java.util.Random(  seedL4);
%a.rndGenL6 = new java.util.Random(  (long) seedL6);;

% =========================================================
% ===============EE delay change==============%
a.spikeDelayActualSteps = PM.prmList{scriptnum}.parameters('spikeDelayActualSteps');
a.minSpikeDelaySteps    = PM.prmList{scriptnum}.parameters('minSpikeDelaySteps');
% ===========================================%
% ===============Ambient conductance timescale============%
a.tauAmbToEa = PM.prmList{scriptnum}.parameters('tauAmbToEa');
a.tauAmbToEb = PM.prmList{scriptnum}.parameters('tauAmbToEb');
a.tauAmbToIa = PM.prmList{scriptnum}.parameters('tauAmbToIa');
a.tauAmbToIb = PM.prmList{scriptnum}.parameters('tauAmbToIb');
%============================================%
%================peak LGN coupled min L6============================%
a.largestMinL6 = PM.prmList{scriptnum}.parameters('largestMinL6');
a.setPeakLGNMinL6();

% ================================================%
%===================filled in=========================%
disp(' ')
disp('---------PARAMETERS---------')
disp(' ')
PM.printParameterSet(scriptnum);
disp(' ')
disp('-----END OF PARAMETERS------');
disp(' ')

%===BUILD NETWORK===%
%===layout the 4C and LGN cells===%
a.layout4C();
a.layoutLGN();
%===create LGN connections===%
a.connectLGN();
%===create 4C connections===%
a.connect4C();
%===load layer 6 recording and set p-values===%
a.loadLayer6Recording();
%===create layer 6 to 4C connections===%
a.connectLayer6();
%===set misc randomized components===%
a.setMiscRandComponents();
%===engage grating===%
a.engageGrating();
%===enforce parameter ratios (like sIEOverSEE)===%
a.enforceRatios();

%==================do new fb build code================%
a.createL4L6Maps;
a.makeLongRangeCoupling;
a.initializeTimeAvgs;
%==================new fb build code=================%

%== copy matfiles into folder ==%
fbInputFile=char(a.fbInputFile);
lgnInputFile=char(a.LGNInputFile);
if ~exist([savefolder fbInputFile])
    try
        copyfile(fbInputFile,savefolder);
    catch err
        disp('error in copyfile line');
    end
end
if ~exist([savefolder lgnInputFile])
    try
        copyfile(lgnInputFile,savefolder);
    catch err
        disp('error in copyfile line');
    end
end


%%

% %===Update l6->l4 connectivity===%
%% get usual cortical locations
numHCRows = 3;
numHC = 3*3;
numExcitatory = 3000 * numHC;
numInhibitory = 1000 * numHC;

locExc = zeros(numExcitatory,2);
locInh = zeros(numInhibitory,2);

%==cortex dimensions==%
cortexHCWidth = 500;
cortexWidth = cortexHCWidth * numHCRows;
cortexNumPerRowExc = ceil(sqrt(numExcitatory));
cortexSpacingExc = cortexWidth / cortexNumPerRowExc;
cortexNumPerRowInh = ceil(sqrt(numInhibitory));
cortexSpacingInh = cortexWidth / cortexNumPerRowInh;

%==fill in cortex locations==%
for i = 1 : numExcitatory
    locExc(i,1) = mod(i-1,cortexNumPerRowExc) * cortexSpacingExc;
    locExc(i,2) = floor((i-1) / cortexNumPerRowExc) * cortexSpacingExc; 
end

for i = 1 : numInhibitory
    locInh(i,1) = mod(i-1,cortexNumPerRowInh) * cortexSpacingInh;
    locInh(i,2) = floor((i-1) / cortexNumPerRowInh) * cortexSpacingInh;
end

%eveningSeed=1; %ORIGINAL
eveningSeed=2; %NEW TEST VALUE
rng(eveningSeed);

loadconnectivities=false;
if loadconnectivities
    disp('Loading connectivities...');
    %load 'connectivitiesMLIndexing.mat';
    load 'connectivitiesMLIndexing5.mat';
    disp('Done.')
else
    a.locFBToPoints();
    a.cornerBoost = 1.0; a.cornerWidth = 0.0;
    [numToAdd, numToDel] = evenOutConnectivityVer2(a.numFBPreExc,locExc,0,'L6->E');
    a.evenOutL6ToE(numToAdd,numToDel,110.0,2.5);
    
    a.cornerBoost = 2.0; a.cornerWidth = 175.0;
    [numToAdd, numToDel] = evenOutConnectivityVer2(a.numFBPreInh,locInh,0,'L6->I');
    a.evenOutL6ToI(numToAdd,numToDel,110.0,2.2);
    
    a.cornerBoost = 1.0; a.cornerWidth = 0.0;
    [numToAdd, numToDel] = evenOutConnectivityVer2(a.numEPreExc,locExc,0,'E->E');
    a.evenOutEToE(numToAdd,numToDel,a.eLenScale,1.0);
    a.cornerBoost = 2.0; a.cornerWidth = 175.0;
    [numToAdd, numToDel] = evenOutConnectivityVer2(a.numIPreExc,locExc,0,'I->E');
    a.evenOutIToE(numToAdd,numToDel,a.iLenScale,1.0);    
    [numToAdd, numToDel] = evenOutConnectivityVer2(a.numEPreInh,locInh,0,'E->I');
    a.evenOutEToI(numToAdd,numToDel,a.eLenScale,1.0);    
    [numToAdd, numToDel] = evenOutConnectivityVer2(a.numIPreInh,locInh,0,'I->I');
    a.evenOutIToI(numToAdd,numToDel,a.iLenScale,1.0);
end

%% remove cpx cells from l6->l4 coupling
disp('Removing cpx cells from l6-l4 coupling...')
l4L6CpxToRemove = PM.prmList{scriptnum}.parameters('l4L6CpxToRemove');

l6L4Map = a.l6L4Map + 1; %convert to MATLAB indexing
numL6L4Map = a.numL6L4Map;

simpleCells = a.numInputsExc >= 3;
cpxCells = a.numInputsExc < 3;

survivalProb = double(simpleCells) + (1-l4L6CpxToRemove).*double(cpxCells);
survivalTest = rand(size(l6L4Map));
survivors = survivalTest < survivalProb(l6L4Map);

disp('deleting...')
rowLen = size(l6L4Map,2);
for i=1:size(l6L4Map,1)
    rowSurvivors = l6L4Map(i,survivors(i,1:numL6L4Map(i)));
    nsurvivors = length(rowSurvivors);
    l6L4Map(i,:) = [rowSurvivors, zeros(1,rowLen-nsurvivors)];
    numL6L4Map(i) = nsurvivors;
end
disp('done deleting.')
disp('reversing map')
[l4L6Map, numL4L6Map]=getReverseConnectivity(l6L4Map,numL6L4Map,a.numExcitatory,120);

disp('implementing...')
a.l6L4Map = l6L4Map-1;
a.numL6L4Map = numL6L4Map;
a.l4L6Map = l4L6Map-1;
a.numL4L6Map = numL4L6Map;
disp('done.')

% %% convert some cpx to simple cells in 135 deg optimal region
% 
% seed = 10;
% rng(seed);
% ballRadius = 100;
% ballCenter = [750 500];
% numToRemove = PM.prmList{scriptnum}.parameters('numPreToRemove');
% 
% distances = rowmap(@(loc) norm(loc-ballCenter),locExc);
% 
% cpxInBall    = find(distances < ballRadius & a.numInputsExc==2 & ~a.randomPicksExc);
% simpleInBall = find(distances < ballRadius & a.numInputsExc==4 & ~a.randomPicksExc);
% 
% %change L4 egal
% numToDel = numToRemove;
% delVector = zeros(1,a.numExcitatory);
% addVector = zeros(1,a.numExcitatory);
% delVector(cpxInBall) = numToDel;
% a.evenOutEToE(addVector,delVector,100.0,1.0); %last two inputs unused
% 
% %change L6 egal    
% numToDel = numToRemove;
% delVector = zeros(1,a.numExcitatory);
% addVector = zeros(1,a.numExcitatory);
% delVector(cpxInBall) = numToDel;
% a.evenOutL6ToE(addVector,delVector,100.0,1.0); %last two inputs unused

%% add some E->I in vertical region

%select I cells in all vertical patches (stay 10 away frm bdy)
type='I';
angleRad=0;
sweepRad=pi/4;
distFmBdy=10;
sel = selectSweep(a, 'I', angleRad, sweepRad, distFmBdy);

%get number to add
numToAdd = PM.prmList{scriptnum}.parameters('EToIBoost');
delVector = zeros(1,a.numInhibitory);
addVector = zeros(1,a.numInhibitory);
addVector(sel) = randi(numToAdd+[-1 1], size(sel));
SD=100;
numSD=1.0;
a.evenOutEToI(addVector,delVector,SD,numSD); %last two inputs unused

%% add some L6->E in 135 deg region

%select E cells in all 135 patches (stay 10 away frm bdy)
type='E';
angleRad=1*pi/2;
sweepRad=pi/4;
distFmBdy=10;
sel = selectSweep(a, type, angleRad, sweepRad, distFmBdy);

%get number to add
numToAdd = PM.prmList{scriptnum}.parameters('135L6ToEBoost');
delVector = zeros(1,a.numInhibitory);
addVector = zeros(1,a.numInhibitory);
addVector(sel) = numToAdd;
addVector(sel) = randRound(addVector(sel));
SD=100;
numSD=1.0;
a.evenOutL6ToE(addVector,delVector,SD,numSD); %last two inputs unused

%% dec some E->I in 135 deg region

%select E cells in all 135 patches (stay 10 away frm bdy)
type='I';
angleRad=1*pi/2;
sweepRad=pi/4;
distFmBdy=10;
sel = selectSweep(a, type, angleRad, sweepRad, distFmBdy);

%get number to add
numToDel = PM.prmList{scriptnum}.parameters('135EToIDec');
delVector = zeros(1,a.numInhibitory);
addVector = zeros(1,a.numInhibitory);
delVector(sel) = numToDel;
delVector(sel) = randRound(delVector(sel));
SD=100;
numSD=1.0;
a.evenOutEToI(addVector,delVector,SD,numSD); %last two inputs unused

%%

%===RUN===%
if recording    
    runAndRecordJavaSim(a, scriptnum, fullsavename, PM, waitDurationInMs, runDurationInMs, recordWhat)
end    

