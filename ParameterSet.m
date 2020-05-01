classdef ParameterSet < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        parameters %map of parameter name-parameter value pairs
    end
    
    methods
        function PS = ParameterSet(prms)          
            PS.parameters = copyMap(prms);
        end
        
        function setKeyVal(PS,key,val)
            PS.parameters(key) = val;            
        end
        
        function areEqual = equals(PS,otherPS)
            %note: this only sees if the keys of the first equal the keys
            %of the second, however second can have more keys
            
            areEqual = true; %until it's not
            
            allkeys = PS.parameters.keys;
            for i = 1 : length(allkeys)
                key = allkeys{i};
                                
                if otherPS.parameters.isKey(key)
                    valueHere = PS.parameters(key);
                    valueThere = otherPS.parameters(key);                    
                else
                    areEqual = false;
                    break;
                end
                
                if ~isequal(valueHere,valueThere)
                        areEqual = false;
                    break;
                end
            end
        end
        
        function areEqual = equalsMod(PS,otherPS,modKeyList)
            %note: this only sees if the keys of the first equal the keys
            %of the second, however second can have more keys
            
            areEqual = true; %until it's not
            
            allkeys = PS.parameters.keys;
            for i = 1 : length(allkeys)
                key = allkeys{i};
                
                onModKeyList = false;
                for j = 1 : length(modKeyList)
                    if isequal(key,modKeyList{j})
                        onModKeyList = true;
                        break;
                    end
                end
                
                if onModKeyList
                    continue;
                end
                
                valueHere = PS.parameters(key);
                if otherPS.parameters.isKey(key)
                    valueThere = otherPS.parameters(key);
                else
                    areEqual = false;
                    break;
                end
                if ~isequal(valueHere, valueThere)
                    areEqual = false;
                    break;
                end                
            end
        end
        
%         function retval = getIndexes(obj, valStrings, values)
%             indexAccum = obj.parameters(:,getIndex(valStrings{1})) == values{1};
%             for i = 2 : length(valStrings)
%                 indexAccum = indexAccum & obj.parameters(:,getIndex(valStrings{i})) == values{i};
%             end
%             retval = find(indexAccum);
%         end


        
        function printParameters(PS)
            prms = PS.parameters;
            ks = prms.keys;
            for i = 1 : length(ks)
                val = prms(ks{i});
                if isnumeric(val)
                    disp([ks{i} ': ' num2str(val)]);
                elseif islogical(val)
                    if val
                        disp([ks{i} ': true'])
                    else
                        disp([ks{i} ': false'])
                    end
                else
                    disp([ks{i} ': ' val]);
                end
            end
        end
        
        function printCertainParameters(PS,keyList) %keyList is a cell containing all the keys to show
            prms = PS.parameters;           
            for i = 1 : length(keyList)
                key = keyList{i};
                val = prms(key);
                if isnumeric(val)
                    disp([key ': ' num2str(val)]);
                elseif islogical(val)
                    if val
                        disp([key ': true'])
                    else
                        disp([key ': false'])
                    end
                else
                    disp([key ': ' val]);
                end
            end
        end
        
        function copyAllParams(PS, javaObj)
            prms = PS.parameters;
            javaObj.sEE=prms('sEE');
            javaObj.sEIOverSEE=prms('sEIOverSEE');
            javaObj.sIEOverSEE=prms('sIEOverSEE');
            javaObj.sELGNOverSEE=prms('sELGNOverSEE');
            javaObj.sILGNOverSEE=prms('sILGNOverSEE');
            javaObj.aChSELGNFac=prms('aChSELGNFac');            
            v=prms('sIIOverSEIBounds'); %'v=' is correct here; check next two lines
            javaObj.sIIOverSEILow=v(1);
            javaObj.sIIOverSEIHigh=v(2);            
            
            javaObj.rateEAmb=prms('rateEAmb');
            javaObj.rateIAmb=prms('rateIAmb');
            javaObj.threshLGN=prms('lGNThresh');
            javaObj.shotNoise=prms('sLGNNoise');
            javaObj.contrast=prms('contrast');
            javaObj.spFreq=prms('sf');
            javaObj.tf=prms('tf') / 1000; %convert Hz to cyc/ms
            javaObj.angle=prms('angle');
            javaObj.pValues=prms('pValues');
            
            javaObj.randomPicksEEConnectivity=prms('randomPicksEEConnectivity');
            javaObj.threshIBounds=prms('threshIBounds');
            javaObj.threshEBounds=prms('threshEBounds');
            javaObj.egal=prms('egal');
            javaObj.peakCEI=prms('peakCEI');
            javaObj.peakCIE=prms('peakCIE');
            javaObj.peakCII=prms('peakCII');
            javaObj.sModSTD=prms('sModSTD');
            javaObj.fracSimple=prms('fracSimple');
            javaObj.simpleAvg=prms('simpleAvg');
            javaObj.simpleAmpLow=prms('simpleAmpLow');
            javaObj.simpleAmpHigh=prms('simpleAmpHigh');
            javaObj.deregAmpLow=prms('deregAmpLow');
            javaObj.deregAmpHigh=prms('deregAmpHigh');
            %javaObj.fbEgal=prms('fbEgal'); 
%             javaObj.fbWideSDE=prms('fbWideSDE');
%             javaObj.fbWideNumSDsE=prms('fbWideNumSDsE');
%             javaObj.fbNarrowSDE=prms('fbNarrowSDE');
%             javaObj.fbNarrowNumSDsE=prms('fbNarrowNumSDsE');
%             javaObj.fbWidePeakCutoffE=prms('fbWidePeakCutoffE');
%             javaObj.fbWideSDI=prms('fbWideSDI');
%             javaObj.fbWideNumSDsI=prms('fbWideNumSDsI');
%             javaObj.fbNarrowSDI=prms('fbNarrowSDI');
%             javaObj.fbNarrowNumSDsI=prms('fbNarrowNumSDsI');
%             javaObj.fbWidePeakCutoffI=prms('fbWidePeakCutoffI');
            javaObj.fbWideSDE=prms('fbWideSDE');
            javaObj.fbWideNumSDsE=prms('fbWideNumSDsE');
            javaObj.fbNarrowSDE=prms('fbNarrowSDE');
            javaObj.fbNarrowNumSDsE=prms('fbNarrowNumSDsE');
            javaObj.fbWideSDI=prms('fbWideSDI');
            javaObj.fbWideNumSDsI=prms('fbWideNumSDsI');
            javaObj.fbNarrowSDI=prms('fbNarrowSDI');
            javaObj.fbNarrowNumSDsI=prms('fbNarrowNumSDsI');
            javaObj.sEEL6Low=prms('sEEL6Low');
            javaObj.sEEL6High=prms('sEEL6High');
            javaObj.sIEL6OverSEEL6FracL4=prms('sIEL6OverSEEL6FracL4');
            javaObj.fbWideEgal=prms('fbWideEgal');
            javaObj.fbNarrowEgal=prms('fbNarrowEgal');
%            javaObj.fbIExpNumPre=prms('fbIExpNumPre');
            javaObj.fbIExpNumPreWide=prms('fbIExpNumPreWide');
            javaObj.fbIExpNumPreNarrow=prms('fbIExpNumPreNarrow');     
            javaObj.numPreDistCutoff=prms('numPreDistCutoff');
            javaObj.numPreDistCutoffFB=prms('numPreDistCutoffFB');
%            javaObj.numPreDistCutoffInner=prms('numPreDistCutoffInner');
%            javaObj.numPreDistCutoffOuter=prms('numPreDistCutoffOuter');
                        
            javaObj.sEEHighMultiplier=prms('sEEHighMultiplier');
            javaObj.sEELowMultiplier=prms('sEELowMultiplier');
            
%             javaObj.dtExc=prms('dtExc');
%             javaObj.dtInh=prms('dtInh');

            javaObj.seedL4=prms('seedL4');
            javaObj.seedL6=prms('seedL6');
            javaObj.seedMisc=prms('seedMisc');
            javaObj.seedRun=prms('seedRun');
            javaObj.seedEngine=prms('seedEngine');
            
            javaObj.sEIDepPerSpike=prms('sEIDepPerSpike');
            javaObj.sEIDepDelta=prms('sEIDepDelta');
            javaObj.sEIDepLookback=prms('sEIDepLookback');
            
            javaObj.tauAMPAToEa=prms('tauAMPAToEa');
            javaObj.tauAMPAToIa=prms('tauAMPAToIa');
            javaObj.tauAMPAToEb=prms('tauAMPAToEb');
            javaObj.tauAMPAToIb=prms('tauAMPAToIb');
            javaObj.tauGABAa=prms('tauGABAa');
            javaObj.tauGABAb=prms('tauGABAb');
            javaObj.tauGABAc=prms('tauGABAc');
            javaObj.tauGABAd=prms('tauGABAd');
            javaObj.weightToSecondGABA=prms('weightToSecondGABA');
            
            javaObj.spikeDelaySteps=prms('spikeDelaySteps');

            javaObj.simpleSpikeProb=prms('simpleSpikeProb');
            javaObj.highCpxProb=prms('highCpxProb');
            javaObj.lowCpxProb=prms('lowCpxProb');

            javaObj.fracNMDAEE=prms('fracNMDAEE');
            javaObj.fracNMDAIE=prms('fracNMDAIE');
            %javaObj.fracAMPAEE=prms('fracAMPAEE');
            %javaObj.fracAMPAIE=prms('fracAMPAIE');
            
% %             javaObj.minSpikeDelayMs=prms('minSpikeDelayMs');
            
            %javaObj.fastSpikeInitial=prms('fastSpikeInitial');
            
            %javaObj.recordingSmearSD=prms('recordingSmearSD');
% %             javaObj.l6L4FRRatio=prms('l6L4FRRatio');
% %             javaObj.recordingBigSmearSD=prms('recordingBigSmearSD');
            %javaObj.recordingBigSmearPercent=prms('recordingBigSmearPercent');
% %             javaObj.sigmaFB=prms('sigmaFB');
% %             javaObj.sigmaFBPast=prms('sigmaFBPast');
% %             javaObj.sigmaNoiseFB=prms('sigmaNoiseFB');
% %             javaObj.sEgal=prms('sEgal');
% %             javaObj.highFiringComplexPerc=prms('highFiringComplexPerc');
            %javaObj.l6PoissonSpPerSec=prms('l6PoissonSpPerSec');
% %             javaObj.fbSpreadLow=prms('fbSpreadLow');
% %             javaObj.fbSpreadHigh=prms('fbSpreadHigh');
% %             javaObj.fbSpikeSpreadProbChanges=prms('fbSpikeSpreadProbChanges');
% %             javaObj.fbSpikeSpreadProbs=prms('fbSpikeSpreadProbs');
% %             javaObj.l6HearsFromSimpleOnly=prms('l6HearsFromSimpleOnly');            
% %             javaObj.constFBFrac=prms('constFBFrac');
% %             javaObj.constL6L4FRRatio=prms('constL6L4FRRatio');
% %             javaObj.puffLow=prms('puffLow');
% %             javaObj.puffHigh=prms('puffHigh');
% %             javaObj.puffToIOverPuffToE=prms('puffToIOverPuffToE');
% %             javaObj.rateENoiseAmb=prms('rateENoiseAmb');
% %             javaObj.rateINoiseAmb=prms('rateINoiseAmb');
% %             javaObj.puffDurationLow=prms('puffDurationLow');
% %             javaObj.puffDurationHigh=prms('puffDurationHigh');
            
% %             javaObj.puffStateDurations=prms('puffStateDurations');
% %             javaObj.puffStateHeights=prms('puffStateHeights');

% %             javaObj.tauAmb=prms('tauAmb');
% %             javaObj.tauNoiseAmb=prms('tauNoiseAmb');
% %             javaObj.sAmb=prms('sAmb');
% %             javaObj.sNoiseAmb=prms('sNoiseAmb');
% %             javaObj.takeAwayPuff=prms('takeAwayPuff');
            
% %             javaObj.fracOriginalL4L6Coupling=prms('fracOriginalL4L6Coupling');
% %             javaObj.numL4L6GaussianSDs=prms('numL4L6GaussianSDs');
% %             javaObj.l4L6GaussianSD=prms('l4L6GaussianSD');
% %             javaObj.fullContrastAntiSpikeFrac=prms('fullContrastAntiSpikeFrac');
% %             javaObj.noContrastAntiSpikeFrac=prms('noContrastAntiSpikeFrac');
% %             
% %             javaObj.fullContrastLowCpxRate=prms('fullContrastLowCpxRate');
% %             javaObj.fullContrastHighCpxRate=prms('fullContrastHighCpxRate');
% %             
% %             javaObj.fbContrastExponent=prms('fbContrastExponent');    
% %             
% %             javaObj.facxs=prms('facxs');
% %             javaObj.facys=prms('facys');
% %             javaObj.facslopes=prms('facslopes');
% %             
% %             javaObj.constFacilitationFac=prms('constFacilitationFac');            
% %             
% %             javaObj.numL4InL4L6InnerMap=prms('numL4InL4L6InnerMap');
% %             javaObj.numL4InL4L6OuterMap=prms('numL4InL4L6OuterMap');
% %             javaObj.l4L6MapWidth=prms('l4L6MapWidth');
            
% %             javaObj.vertTemplateMod=prms('vertTemplateMod');
% %             javaObj.horizTemplateMod=prms('horizTemplateMod');
% %             
% %             javaObj.facIOverFacE=prms('facIOverFacE');
% %             javaObj.xL4End=prms('xL4End');
% %             javaObj.xL4Start=prms('xL4Start');
% %             javaObj.fbBGRate=prms('fbBGRate');
% %             
% %             javaObj.sEIDepLookback=prms('sEIDepLookback');            
% %             javaObj.sEIDepDelta=prms('sEIDepDelta');
% %             
% %             javaObj.l6L4FRRatioStart=prms('l6L4FRRatioStart');
% %             
% %             javaObj.l6L4xs=prms('l6L4xs');
% %             javaObj.l6L4ys=prms('l6L4ys');
% %             javaObj.l6L4slopes=prms('l6L4slopes');
% %             
% %             javaObj.sEIDepPerSpike=prms('sEIDepPerSpike');
% %             
% %             javaObj.fbPastFrac=prms('fbPastFrac');
% %             
% %             javaObj.minL6MaxL6Ratio=prms('minL6MaxL6Ratio');
% %             
% %             javaObj.eFatAmt=prms('eFatAmt');
% %             javaObj.eFatLookback=prms('eFatLookback');
% %             
            %javaObj.strengthEgal=prms('strengthEgal');
            %javaObj.strengthEgalFB=prms('strengthEgalFB');
%             javaObj.sELGNOverSEERnd=prms('sELGNOverSEERnd');
            
%             javaObj.fbFlatE=prms('fbFlatE');
%             javaObj.fbFlatI=prms('fbFlatI');
%             javaObj.fbESlopeSD=prms('fbESlopeSD');
%             javaObj.fbISlopeSD=prms('fbISlopeSD');
%             javaObj.fbESlopeCutoffSDs=prms('fbESlopeCutoffSDs');
%             javaObj.fbISlopeCutoffSDs=prms('fbISlopeCutoffSDs');
            
            %javaObj.fbStdSqr=prms('fbStdSqr');
            %javaObj.fbESD=prms('fbESD');
            %javaObj.fbISD=prms('fbISD');
            
            %javaObj.rhoE=prms('rhoE');
            %javaObj.rhoI=prms('rhoI');
            %javaObj.flashDurationMs=prms('flashDurationMs');
            %javaObj.flashMultiplierFB=prms('flashMultiplierFB');
            %javaObj.gLInh = javaObj.gLExc * prms('leakIOverLeakE');

            javaObj.load4CConnections=prms('load4CConnections');
            javaObj.l4CConnectionsFile=prms('l4CConnectionsFile');
            javaObj.loadFBConnections=prms('loadFBConnections');
            javaObj.fBConnectionsFile=prms('fBConnectionsFile');
            javaObj.loadMiscRandComponents=prms('loadMiscRandComponents');
            javaObj.miscRandComponentsFile=prms('miscRandComponentsFile');
            javaObj.LGNInputFile=prms('LGNInputFile');
            javaObj.fbInputFile=prms('fbInputFile');
            javaObj.fbStart=prms('fbStart');
            javaObj.fbEnd=prms('fbEnd');            
            javaObj.fbRandEnd=prms('fbEnd') - 50 * 0.814723686393179;
            javaObj.fbClock=prms('fbClock');
            javaObj.currentFBSpike=prms('currentFBSpikes');
            
            javaObj.numHCRows=prms('numHCRows');
            javaObj.numExcitatoryPerHC=prms('numExcitatoryPerHC');
            javaObj.numInhibitoryPerHC=prms('numInhibitoryPerHC');
            javaObj.numLGNBlockRows=prms('numLGNBlockRows');
            javaObj.numLGNPerRowPerBlock=prms('numLGNPerRowPerBlock');
            javaObj.numLGNPerColPerBlock=prms('numLGNPerColPerBlock');
        end
        
        %The format of keyValList is {{key1 val1} {key2 val2} {...} ... }
        function huh = match(PS,keyValList)
            huh = true;
            for j = 1 : length(keyValList)
                key = keyValList{j}{1};
                val = keyValList{j}{2};
                if ~isequal(single(PS.parameters(key)),single(val))
                    huh = false;
                end
            end
        end
        
%         function copyParam(PS, paramKey, javaObj)
%             prms = PS.parameters;
%             ks = prms.keys;
%             switch(lower(paramKey))
%                 case 'see'
%                     val = prms(paramKey);
%                     javaObj.setSEE(val);
%                 case 'sieoversee'
%                     val = prms(paramKey);
%                     javaObj.setSIEOverSEE(val);
%                 case 'seioversee'
%                     val = prms(paramKey);
%                     javaObj.setSEIOverSEE(val);
%                 case 'siioverseibounds'
%                     val=prms('sIIOverSEIBounds');                    
%                     javaObj.setModIIRange(val(1),val(2));
%                 case 'selgnoversee'
%                     val=prms('sELGNOverSEE');
%                     javaObj.setSELGNOverSEE(val);
%                 case 'silgnoversee' 
%                     val=prms('sILGNOverSEE');
%                     javaObj.setSILGNOverSEE(val);
%                 case 'tauampa'
%                     val=prms('tauAMPA');
%                     javaObj.tauAMPA = val;
%                 case 'taugaba'
%                     val=prms('tauGABA');
%                     javaObj.tauGABA = val;                    
%                 case 'rateeamb'
%                     val1 = prms('rateEAmb');
%                     val2 = prms('rateIAmb');
%                     javaObj.setAmbient(val1,val2);
%                 case 'rateiamb'
%                     val1 = prms('ambe');
%                     val2 = prms('ambi');
%                     javaObj.setAmbient(val1,val2);
%                 case 'contrast'
%                     val = prms('contrast');
%                     javaObj.setContrast(val);
%                 case 'iThreshLow'
%                     val1 = prms('iThreshLow');
%                     val2 = prms('iThreshHigh');
%                     javaObj.setThreshIRange(val1, val2);
%                 case 'iThreshHigh'
%                     val1 = prms('iThreshLow');
%                     val2 = prms('iThreshHigh');
%                     javaObj.setThreshIRange(val1, val2);
%                 case 'lgnthresh'
%                     val = prms('lGNThresh');
%                     javaObj.threshLGN = val;
%                 case 'slgnnoise'
%                     val = prms('sLGNNoise');
%                     javaObj.shotNoise = val;
%                 case 'sf'
%                     val = prms('sf');
%                     javaObj.setSpFreq(val);
%                 case 'tf'
%                     val = prms('tf');
%                     javaObj.setTF(val);
%                 case 'LGNInputFile'
%                     val = prms('LGNInputFile');
%                     load(val,'randomPicksExc','randomPicksInh');
%                     javaObj.randomPicksExc = randomPicksExc;
%                     javaObj.randomPicksInh = zeros(numInhibitory,1);
%                     javaObj.LGNInputFile = val;                    
%                 case 'cEE'
%                     javaObj.
%                 case 'cIE'
%                     javaObj.
%                 case 'cEI'
%                     javaObj.
%                 case 'cII'
%                     javaObj.
%                 case 'egal'
%                     javaObj.
%                 case 'LGNBaseRate'
%                     javaObj.
%                 case 'bgRate'
%                     javaObj.
%                 case 'fbE'
%                     javaObj.
%                 case 'fbI'
%                     javaObj.
%                 case 'orientationTheta'
%                     javaObj.
%                 case 'spatialFreq'
%                     javaObj.
%                 default
%                     throw MException('ParameterSet:ArgError','parameter given is not a known name');
%             end
%         end
    end % methods
    
    methods(Static)
        function index = getIndex( request )
            %UNTITLED Summary of this function goes here
            %   Detailed explanation goes here
            dictionary = {...
                'ambE',...
                'ambSmall',...
                'sEE',...
                'sEI',...
                'sIE',...
                'drive',...
                'iThreshLow',...
                'iThreshHigh',...
                'sIILow',...
                'sIIHigh',...
                'LGNInputFile',...
                'cEE',...
                'cIE',...
                'cEI',...
                'cII',...
                'egal',...
                'LGNBaseRate',...
                'bgRate',...
                'fbE',...
                'fbI',...
                'orientationTheta',...
                'spatialFreq',...
                };
            index = 0;
            for i=1 : length(dictionary)
                if strcmpi(dictionary{i},request)
                    index = i;
                end
            end
            
        end
    end % static methods
    
end
