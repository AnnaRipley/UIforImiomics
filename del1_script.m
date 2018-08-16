% Prepare data needed for the UI part 1
% Thesis project, Anna Larsson 2018. See Readme for instructions.

load('uidata');

% Zoom radius
radius = 20;

% Custom colormap BlueWhiteRed
blue = [0 0 1];
white = [1 1 1];
red = [1 0 0];
blueinc = [1 1 0];
redinc = [0 1 1];
BWRmap = [ blue
    blue + 0.1* blueinc
    blue + 0.2* blueinc
    blue + 0.3* blueinc
    blue + 0.4* blueinc
    blue + 0.5* blueinc
    blue + 0.6* blueinc
    blue + 0.7* blueinc
    blue + 0.8* blueinc
    blue + 0.9* blueinc
    white
    red + 0.9* redinc
    red + 0.8* redinc
    red + 0.7* redinc
    red + 0.6* redinc
    red + 0.5* redinc
    red + 0.4* redinc
    red + 0.3* redinc
    red + 0.2* redinc
    red + 0.1* redinc
    red];