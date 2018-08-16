% Run to read volume data files and create a (very) large .mat-file
% Runtime: 10 min
% Thesis project, Anna Larsson 2018. See Readme for instructions.

% Path definitions.
% When using a new setup, change these if needed.
pathCorrelation = 'E:/imiomics/correlation maps/';
pathDeformed = 'E:/imiomics/POEM data/';
pathUndeformed = 'E:/imiomics/POEM data undeformed/';

saveDataAs = 'E:\imiomics\uidata.mat';

load('IDs.mat');
load('POEMdxa.mat');


idlength = 0;
idtablelength = size(IDs);
for row = 1 : idtablelength(1)
    if IDs{row, 2}
        idlength = idlength + 1;
    end
end
idlength;
IDarrayLonger = zeros(idlength, 1);
i = 1;
for row = 1 : idtablelength(1)
    if IDs{row, 2}
        IDarrayLonger(i) = IDs{row, 1};
        i = i + 1;
    end
end

poemlength = 0;
poemtablelength = size(POEMdxa);

for row = 1 : poemtablelength(1)
    if( ~isnan(POEMdxa{row,4}) && ismember(POEMdxa{row,1}, IDarrayLonger) )
        poemlength = poemlength + 1;
    end
end
% use poemlength to pre-allocate, to save time
dxavals = zeros(poemlength, 1);
i = 1;
% Shorter IDarray: only IDs that are in the ID list AND have dxa data 
IDarray = zeros(poemlength, 1);
% Initialize imiFatVolarray.
imiFatVolarray(poemlength) = imiFatVol();
for row = 1 : poemtablelength(1)
    if( ~isnan(POEMdxa{row,4}) && ismember(POEMdxa{row,1}, IDarrayLonger) )
        thisval = POEMdxa{row,4};
        dxavals(i) = thisval;
        thisID = POEMdxa{row,1};
        IDarray(i) = thisID;
        thisFem = POEMdxa{row,2};
        imiFatVolarray(i).poemid = thisID;
        imiFatVolarray(i).female = thisFem;
        imiFatVolarray(i).dxaval = thisval;
        i = i + 1;
    end
end

% prepare data for gui del1


D = dir(pathDeformed);
sizeDir = length(D);

noFatvols = 0;
noJacDets = 0;

for i = 1 : sizeDir
    strr = D(i).name;
    if(length(strr) > 5)
        if(contains(strr, 'fat_'))
            noFatvols = noFatvols + 1;
        else
            if(contains(strr, 'JacDet_'))
                noJacDets = noJacDets + 1;
                info = vtk_read_header(strr);
            end
                
        end
    end
end
if noFatvols == noJacDets
    sprintf("All is well.")
else
    sprintf("Numbers of 'fat' and 'JacDet' are not equal!")
end

dims = info.Dimensions;

xmax = dims(1);
ymax = dims(2);
zmax = dims(3);

dimarray = info.PixelDimensions;
normdims = dimarray / min(dimarray);
xscale = normdims(1);
yscale = normdims(2);
zscale = normdims(3);

corscale = [yscale xscale 1];
sagscale = [yscale zscale 1];
axiscale = [xscale zscale 1];

fatscale = 100;
jdscale = 20;


for i = 1 : sizeDir
    strr = D(i).name;
    if(length(strr) > 9)
        thispoemid = str2double( strr(end-9 : end-4 ) );
        if ismember(thispoemid, IDarray) 

            if contains(strr, 'fat_')
                idnotfound = 1;
                di = 1;
                while(idnotfound && di <= poemlength)
                    check = imiFatVolarray(di).poemid;
                    if thispoemid == check
                        sprintf(thispoemid + " (in filename) equals? \n" + check + " (in struct)")
                        idnotfound = 0;
                        info = vtk_read_header([pathDeformed, strr]);
                        imiFatVolarray(di).DefVol = int8(fatscale * vtk_read_volume(info));
                    else
                        di = di + 1;
                    end
                    
                end

            else
                if contains(strr, 'JacDet_')
                    idnotfound = 1;
                    ji = 1;
                    while(idnotfound && ji <= poemlength)
                        check = imiFatVolarray(ji).poemid;
                        if thispoemid == check
                            sprintf(thispoemid + " (in filename) equals? \n" + check + " (in struct)")
                            idnotfound = 0;                    
                            info = vtk_read_header([pathDeformed, strr]);
                            imiFatVolarray(ji).JacDet = int8(jdscale * vtk_read_volume(info));
                        else
                            ji = ji + 1;
                        end
                    end
                end
            end
                
        end
    end
end
% reads all file names in alphabetic order, first comes files 
% '.' and '..' then JacDet_ ...      
    
xinit = round(xmax/2);
yinit = round(ymax/2);
zinit = round(zmax/2);
qinit = 1;

noFemales = 0;
noMales = 0;
for i = 1:poemlength
    if imiFatVolarray(i).female
        noFemales = noFemales + 1;
    else
        noMales = noMales + 1;
    end
end

fi = 1;
mi = 1;
startrowFemale = zeros(noFemales, 1, 'int8');
startrowMale = zeros(noMales, 1, 'int8');
dxavalsFemale = zeros(noFemales, 1);
dxavalsMale = zeros(noMales, 1);
poemIDsFemale = zeros(noFemales, 2);
poemIDsMale = zeros(noMales, 2);
for i = 1:poemlength
    if imiFatVolarray(i).female
        startrowFemale(fi) = imiFatVolarray(i).DefVol(xinit, yinit, zinit);
        dxavalsFemale(fi) = imiFatVolarray(i).dxaval;
        poemIDsFemale(fi,1) = imiFatVolarray(i).poemid;
        poemIDsFemale(fi,2) = i;
        fi = fi + 1;
    else
        startrowMale(mi) = imiFatVolarray(i).DefVol(xinit, yinit, zinit); 
        dxavalsMale(mi) = imiFatVolarray(i).dxaval;
        poemIDsMale(mi,1) = imiFatVolarray(i).poemid;
        poemIDsMale(mi,2) = i;
        mi = mi + 1;
    end
end

A = dir(pathCorrelation);
sizeDir = length(A);

noFiles = 0;

for i = 1 : sizeDir
    strr = A(i).name;
    if(length(strr) > 5)
        noFiles = noFiles + 1;
    end
end


% Make data matrix with colums
% female/male   pval/rval   fat/JacDet
datamatrix = zeros(noFiles, 3);

% Make array of file names
filenames(noFiles) = "";

corrVol = zeros(xmax, ymax, zmax, noFiles); %, 'int8');

ni = 1;
for i = 1 : sizeDir
    strr = A(i).name;
    if(length(strr) > 5)
        filenames(ni) = strr;
        if(contains(strr, 'female'))
            datamatrix(ni, 1) = 1;
        end
        if(contains(strr, 'pval'))
            datamatrix(ni, 2) = 1;
        end
        if(contains(strr, 'fat'))
            datamatrix(ni, 3) = 1;
        end
        
        info = vtk_read_header(strr);
        corrVol(:,:,:,ni) = vtk_read_volume(info);
        
        ni = ni + 1;
    end
end

imshow(corrVol(:,:,zinit,4 ));


% Save all data, path specified in the beginning of this script.
save(saveDataAs);
