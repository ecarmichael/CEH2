function behav_out = MS_behav_append_Coord(behav, fieldname)
%% MS_behav_append_Coord: lets the user draw the idealized trajectory and appends the coordinates to the behav structure. It also appends the linearized position. 
%
%
%
%
%     Inputs:
%       - behav: [struct]  behavioural tracking data in the Miniscope
%       format. 
%
%       - fieldname: [string]  what you would like the appended fieldname
%       to be called
%
%
%     Outputs:
%       - behav_out: [struct] same as the behav input but with an
%       additional field (fieldname) containing the coord struct. 
%
%% initialize

if nargin < 2
    fieldname = 'CoorD';
end

behav_out = behav; 

% make the coord
behav_out.(fieldname) = MakeCoord(MS_behav2tsd(behav), 'titl', fieldname);

% behav_out.([fieldname '_lin']) = LinearizePos([], MS_behav2tsd(behav),behav_out.(fieldname));
