% create data for del 2
% Thesis project, Anna Larsson 2018. See Readme for instructions.

% Path definitions.
% When using a new setup, change these if needed.
pathDeformed = 'E:/imiomics/POEM data/';
pathUndeformed = 'E:/imiomics/POEM data undeformed/';

load('IDselected.mat');
load('POEMdxa.mat');
load('refFemaleID.mat');
load('refMaleID.mat');

U = dir(pathUndeformed);
sizeU = length(U);

for i = 1 : sizeU
    strr = U(i).name;

    if contains(strr, refMaleID)
        info = vtk_read_header(strr);
        refvolMale = vtk_read_volume(info);
    end
    if contains(strr, refFemaleID)
        info = vtk_read_header(strr);
        refvolFemale = vtk_read_volume(info);
    end
end


dims = info.Dimensions; %size(deformvol);

xmax = dims(1);
ymax = dims(2);
zmax = dims(3);

xinit = round(xmax/2);
yinit = round(ymax/2);
zinit = round(zmax/2);


dimarray = info.PixelDimensions;
normdims = dimarray / min(dimarray);
xscale = normdims(1);
yscale = normdims(2);
zscale = normdims(3);

corscale = [yscale xscale 1];
sagscale = [yscale zscale 1];
axiscale = [xscale zscale 1];


poemsize = size(POEMdxa);
poemlength = poemsize(1);

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

