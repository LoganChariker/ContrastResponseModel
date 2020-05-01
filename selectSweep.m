function selection = selectSweep( a, type, angleRad, sweepRad, distFromBdy )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% get all locations as nx2 vector
locs = getLocs( a, type );

% convert locs to their HC 1 equivalent locations
hc1Locs = convertToHC1( locs );

% select everything within the sweep
angleR = angleRad-sweepRad/2+pi/2;
angleL = angleRad+sweepRad/2-pi/2;
vecR = [cos(angleR); sin(angleR)];
vecL = [cos(angleL); sin(angleL)];
selection = find( hc1Locs * vecR > distFromBdy & ...
                  hc1Locs * vecL > distFromBdy );

end

function ar = NeuronLoc2Vec( NL )
    ar = [NL.x NL.y];
end

function ar = NeuronLocs2Mat( NLAr )
    ar = rowmap( @NeuronLoc2Vec, NLAr(:) );
end

function locs = getLocs( a, type )
    switch type
        case 'E'
            locs = NeuronLocs2Mat(a.locExc);
        case 'I'
            locs = NeuronLocs2Mat(a.locInh);
        case 'FB'
            locs = a.locFB;
    end
end

function hc1Locs = convertToHC1( locs )
    HCWidth = 500;
    
    %HC coordinate
    HCCoord = floor( locs / HCWidth );
    
    %HC centers
    HCCenters = [250,250] + HCCoord*500;
    
    %coordinate relative to pinwheel
    pinwheelCoord = locs-HCCenters;
    
    %flip the pinwheel coordinates for every HC away from 1st
    hc1Locs = pinwheelCoord .* (-1) .^ HCCoord; 
end    