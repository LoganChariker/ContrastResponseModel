classdef ParameterManager < handle
    %ParameterManager Keeps sets of parameters organized
    %   
    
    properties
        stdPrms %standard parameters, from which modifications and parameter sweeps are made
        prmList %list of parameters
    end
    
    methods
        function pm = ParameterManager(varargin)
            pm.prmList = {};
        end
        
        function printParameterSet(PM,scriptnum)
            prmList = PM.prmList;
            PS = prmList{scriptnum};
            PS.printParameters();
        end
        
        function printParameterSetCertainParameters(PM,scriptnum,keyList)
            prmList = PM.prmList;
            PS = prmList{scriptnum};
            PS.printCertainParameters(keyList);
        end
        
        function printAllParameterSets(PM)
            for i = 1 : length(PM.prmList)
                disp(['parameter ' num2str(i) ':'])
                PM.printParameterSet(i);
                disp(' ');
            end
        end
        
        function printAllParameterSetsCertainParameters(PM,keyList)
            for i = 1 : length(PM.prmList)
                disp(['parameter set ' num2str(i) ':'])
                PM.printParameterSetCertainParameters(i,keyList);
                disp(' ');
            end
        end
        
        function copyParameterSetToJava(PM,scriptnum,a)
            prmList = PM.prmList;
            PS = prmList{scriptnum};
            PS.copyAllParams(a);
        end
        
        function setBasePrms(PM,prms)
            PM.stdPrms = ParameterSet(prms);
        end
        
        function printBasePrms(PM)
            PM.stdPrms.printParameters();
        end
        
        function prms = changingParameters(PM)
            allkeys = keys(PM.prmList{1}.parameters);
            
            numValueChanging = 0;
            prms = cell(1,length(allkeys));
            for i = 1 : length(allkeys)
                ky = allkeys{i};
                
                %see if ky changes in value
                changesValue = false;
                firstVal = PM.prmList{1}.parameters(ky);
                for j = 2 : length(PM.prmList)
                    if ~isequal(PM.prmList{j}.parameters(ky), firstVal)
                        %then parameter changes
                        changesValue = true;
                        break;
                    end
                end
                
                %we now know if ky changes value
                if changesValue
                    numValueChanging = numValueChanging + 1;
                    prms{numValueChanging} = ky;
                end
            end
            
            prms = prms(1 : numValueChanging);
        end
        
        %print all values a key has in a cell array
        function vals = parameterValues(PM, key)
            numVals = 0;
            vals = cell(1,length(PM.prmList));
            for i = 1 : length(PM.prmList)
                val = PM.prmList{i}.parameters(key);
                
                %find whether val is new to the list
                valnew = true;
                for j = 1 : numVals
                    if isequal(val,vals{j})
                        valnew = false;
                        break;
                    end
                end
                
                %known if val is new to list here
                if valnew
                    numVals = numVals + 1;
                    vals{numVals} = val;
                end
            end
            
            vals = vals(1:numVals);
        end
        
        function easyViewRun(PM)
            chprms = PM.changingParameters();
            for i = 1 : length(chprms)
                strVals = '';
                pmVals = PM.parameterValues(chprms{i});
                for j = 1 : length(pmVals)
                    if isnumeric(pmVals{j}) || islogical(pmVals{j})
                        strVals = [strVals  num2str(pmVals{j}) ' '];
                    else
                        strVals = [strVals '''' pmVals{j} ''' '];
                    end
                end
                disp(['{''' chprms{i} '''  } %%' strVals]);
            end
        end
                
        
        %The format of keyValList is {{key1 val1} {key2 val2} {...} ... }
        %findScriptnum will create a list of all parameteres with the 
        %listed key-value pairs
        %it will also print the out.
        function scriptnums = findScriptnum(PM,keyValList,varargin)
            nVarargs = length(varargin);
            if nVarargs > 0
                keysToShow = varargin{1};
            end
            
            prmList = PM.prmList;
            
            scriptnums = zeros(length(prmList),1);
            
            for i = 1 : length(prmList)
                PS = prmList{i};
                if PS.match(keyValList)
                    %print the parameter set
                    disp(['Parameter set ' num2str(i)]);
                    if nVarargs > 0
                        PS.printCertainParameters(keysToShow);
                        disp(' ');
                    end
                                        
                    scriptnums(i)=1;
                end
            end         
        end
            
        %In the following function addAllCombinations,
        %the form of prms is a map of key-value pairs where a value
        %is a cell, and addAllCombinations adds to the parameter list prmList
        %all combinations of the cell elements, where each cell element is
        %attached to the key originally attached to its cell.  If any
        %base parameter stdPrms keys are not in prms, then they are added
        %to each combination. Also every combination is made a 
        %ParameterSet.
        function addAllCombinations(PM,prms)
            ks = keys(prms);
                        
            prmListLen = length(PM.prmList);
            startingIndex = prmListLen + 1;
            PM.prmList{startingIndex} = ParameterSet(PM.stdPrms.parameters);
            endingIndex = startingIndex;
            
            for i = 1 : length(ks)
                key = ks{i};
                vals = prms(key);
                currentIndex=startingIndex;
                for j = 1 : length(vals)                    
                    for cyclingIndex = startingIndex : endingIndex
                        val = vals{j};
                        cyclingPrmSet = PM.prmList{cyclingIndex};
                        newPrmSet = ParameterSet(cyclingPrmSet.parameters);
                        newPrmSet.setKeyVal(key,val);
                        PM.prmList{currentIndex} = newPrmSet;
                        
                        currentIndex=currentIndex+1;
                    end
                end
                if ~isempty(vals)
                    endingIndex = currentIndex-1;
                end
            end
        end %addAllCombinations
        
        % give each run a unique group index and output as a vector size
        % #runs x 1
        function groupIndexes = getGroupIndexes(PM, modOutPrms)
            numRuns = length(PM.prmList);
            groupIndexes = zeros(numRuns,1);
            
            %find all changing parameters and delete anything from
            %modOutPrms            
            prmsPre = PM.changingParameters();
            prmsToKeep = ~ismember(prmsPre,modOutPrms);
            prms = prmsPre(prmsToKeep);
            if isempty(prms)
                %no sets, just return with group indexes all 0
                return
            end
                        
            % give each run a parameter index
            prmIndexes = zeros(numRuns,length(prms));
            for i=1:length(prms)
                prm=prms{i};
                
                vals=PM.parameterValues(prm);
                                
                for j=1:numRuns
                    val=PM.prmList{j}.parameters(prm);
                    prmIndex=1;
                    while ~isequal( val, vals{prmIndex} )
                        prmIndex=prmIndex+1;
                    end
                    prmIndexes(j,i)=prmIndex;
                end
            end
            
            % give each run a group index based on unique sets of prmIndexes
            groupIndexes = prmIndexes(:,1)-1;
            
            maxCols = max(prmIndexes);            
            scale = maxCols(1);
            
            for i=2:length(maxCols)
                groupIndexes=groupIndexes + scale * (prmIndexes(:,i)-1);
                scale = scale * maxCols(i);
            end
        end
            
        
    end %methods    
end

