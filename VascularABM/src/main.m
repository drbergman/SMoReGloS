function main(INPUT)

% addpath(".")

[loop,si] = ind2sub([81,6],INPUT);
if ~exist("Data","dir")
    mkdir('Data')
end

n_loops = 300;
print_every = Inf;

exType = 'Binary'; %Gradated or binary
[AA,AB,BB,CC] = ndgrid([0.05,0.125,0.245], [0.01,0.05,0.1], [1,2,3],[8,12,16]); % AA is pdiv; AB is sdiv; and BB is tip migration; CC is division limit
RandStream.setGlobalStream(RandStream('mt19937ar','seed','shuffle')); % Make sure random is truly random

AA_str = regexprep(num2str(AA(loop)),'\.','_');
AB_str = regexprep(num2str(AB(loop)),'\.','_');
BB_str = regexprep(num2str(BB(loop)),'\.','_');
CC_str = regexprep(num2str(CC(loop)),'\.','_');

cd Data %everything occurs in data to avoid clutter of folder BSRI19
if ~exist(exType,"dir")
    mkdir(exType);
end
cd(exType)
experimentfolder = strcat('AA_',AA_str,'__AB_',AB_str);
if ~exist(experimentfolder,"dir")
    mkdir(experimentfolder)
end
cd(experimentfolder)
subfolder = strcat('BB_',BB_str,'__CC_',CC_str,'__Sample',num2str(si));%subfolder to keep track of each combination of variables
if ~exist(subfolder,"dir")
    mkdir(subfolder)
end
cd(subfolder)
% folder_names = ["OverTime","Run_Images","Vasculature_Images"];
% for i = 1:length(folder_names)
%     if ~exist(folder_names(i),"dir")
%         mkdir(folder_names(i))
%     end
% end

namedataA =strcat('NumberofCells','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(n_loops),'.txt');

if exist(sprintf("./fig_data/%s",namedataA),"file")
    cd ../../../..%back to main folder
    return;
end

cd ../../../..%back to main folder
fprintf('in\n')
%global symRate;
cartRandDeath=0.051;
symchange = 0.05;
% runnum = 2; % AB(loop);
seedrate = 0;
sdiv = AB(loop);%.05;
pdiv = AA(loop);
dlim =  CC(loop); %AB(loop);
macronum = 0;
fibronum = 0;% BSRI 19 set to 0
%cartnum = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For now at 20- infiltration later
promultm = 1;%BB(loop);
migmultm = 1;%AA(loop);
cartnum = 10;
cartnum2 = cartnum;

%% The Parameters - Perhaps should be a separate file
%migclock = 3; %Migration occurs every x+1 iterations
migmove = 1; %How many voxels the cell can move
DoMigration = 1; %Do you have migration?
sendeath = 0.1;%because time step is smaller%0.1; %percent  of senescent cells that die each it / Sai - 3/4/14- changed from .1

cartLifespan = 56;
CARTDelay=1500;
cartProlifDuration = 48; %ammount of lifespan that car t cells can prolif

%assigning variables needed for respective experiment types.
if strcmp(exType,'Gradated')
    binaryHeterogeniety = false;
    cartHetMean = 0.9;%AA(loop)/100;
    cartHetStd = .1;
elseif strcmp(exType, 'Binary')
    binaryHeterogeniety = true;
    AntigenHetPercentExpressed = 0.9;%AA(loop)/100;%.9; % the fraction of Cancer cells that express the CAA for car T cells
end

%binaryHeterogeniety = false;

CartBreakpt = 0;

%% Setup Initial Niche
%There will be only one metastatic stem cell in initial metastatic Niche

%Sensitivity analysis storage
cellsovertime = zeros(300,1);
stemsovertime = zeros(300,1);
CCR5overtime = zeros(300,1);
hypoxiaovertime = zeros(300,1);
deathovertime = zeros(300,1);
macroovertime = zeros(300,1);
cartovertime = zeros(300,1);
killsovertime = zeros(300,1);
prolifovertime = zeros(300,1);
vascovertime = zeros(300,1);

%Inputs
%Proliferation case: "allcells" - all cells have the possibility to proliferate each it
%       "onecell" - only one cell proliferate each iteration
% ProlifMethod = 'allcells';
% Seedmode = 'OneDirectionDist';%'EmptyRandom'; %location/placement of the seeds

%% Setup Initial Cells
%Setup CA Grid with Cells
%SetupTumor

%Setup CA Grid
gridsizeC = [50,50,50];

%Initialize the grid
%Guide: 0- empty, 1 - tumor cell, 2 - macrophage, 3 - fibroblast, 4 -hypoxic tumor cell
Agentmat = zeros(gridsizeC(1), gridsizeC(2), gridsizeC(3));

%Make 100 cells
%theCells = zeros(100,1);
%numberofstems = 0;
CellState = zeros(100,1); %1 is stem and 2 is progen and
QState = ones(100,1); %Whether the cell is alive/active (1), quiescent (2), or senescent (0)
%CellSymRate = zeros(500000,1);
CCR5level = zeros(100,1);
MigSpeed = ones(100,1);
MigProb = 0.25*ones(100,1);
XYZ = zeros(100,3);
HypCell = zeros(100,1);
HypCounter = zeros(100,1);
%MigCycle = zeros(500000,1);
DivLimit = dlim*ones(100,1); %Used to be 12 - still is because not based on time
DivRate = pdiv*ones(100,1);
DivCycle = zeros(100,1);
AntigenHet = zeros(100,1);
if binaryHeterogeniety
    r = randperm(100);
    AntigenHet(r(1:AntigenHetPercentExpressed*length(AntigenHet)))=1;
else
    pd = makedist('Normal','mu',cartHetMean,'sigma',cartHetStd);
    AntigenHet = random(pd,length(AntigenHet),1);
    AntigenHet(AntigenHet>1)=1;
    AntigenHet(AntigenHet<0)=0;
end
MacroInfProb = 0;%0.6;%used to be 0.5 %Probability a macrophage will infiltrate a tumor


%Make 100 stems with a 20% chance of being CCR5+
%Define initial metastatic stem cells

numberofCCR5 = 0;

for jj = 1:20
    CellState(jj) = 1;
    DivRate(jj) = sdiv;
    %CellSymRate(jj) = 0.05;
    if rand < 0.2
        CCR5level(jj) = 1;
        MigProb(jj) = 1;
        numberofCCR5= numberofCCR5 + 1;
    end
end
numberofstems = 20; %There is one stem cell in first it


%Make 80 progenitor cells
for jj = 21:100
    CellState(jj) = 2;
    if jj < 26
        CCR5level(jj) = 1;
        %MigSpeed(jj) = migmove;
        numberofCCR5 = numberofCCR5+1;
        MigProb(jj) = 1;
    end
    DivCycle(jj) = ceil(12*rand(1)); %7_16_16 Make Cycle Random
end

numberofcells = 100;
%Place Cells randomly in grid
%Choose cells in random order
rorder = randperm(numberofcells);
%Need 1:5, 95:100, and 1:4
%Setup Meshgrid
[AAx,ABx,BBx] = meshgrid(linspace(1,5,5), fliplr(linspace(1,5,5)), linspace(1,4,4)); %resolution of data here
for loopg = 1:numel(ABx)
    %parfor loop = 1:length(stemloop)
    %global symRate;
    xax = AAx(loopg);
    yax = ABx(loopg);
    zax = BBx(loopg);
    pos = rorder(loopg);
    XYZ(loopg,:) = [xax, yax, zax];
end
%cell1.Position = [1,100,1];%[100,100,100]; %Initial position of cell
Agentmat(1:5,1:5,1:4) = 1; %Place cell on grid V6

%Setup Macrophages
%There are 14450 macrophages randomly spread out
%macronum = macs;
%Setup Macrophages
%There are 14450 macrophages randomly spread out
XYZm = [];
Statem = [];
numberofcellsm = 0;
if macronum > 0
    sizegrid = numel(Agentmat);
    rorder2 = randperm(sizegrid);
    themacro = rorder2(1:macronum);
    Agentmat(themacro) = 2;
    [Is,Js,Ks] = ind2sub(size(Agentmat),themacro);
    theinds = [Is;Js;Ks];
    XYZm = theinds';
    Statem = 2*ones(macronum,1);%unknown if needed
end

%Setup Fibroblasts
XYZf = [];
if fibronum > 0
    sizegrid = numel(Agentmat);
    rorder2 = randperm(sizegrid);dfb bv
    thefibro = rorder2(1:fibronum);
    Agentmat(thefibro) = 3;
    [Ifs,Jfs,Kfs] = ind2sub(size(Agentmat),thefibro);
    theinds = [Ifs;Jfs;Kfs];
    XYZf = theinds';
    Statef = 3*ones(fibronum,1);
end

%Setup Killer T-Cells
%This would place them initially if not delaying therapy
XYZcart = [];
%     if cartnum > 0
%         sizegrid = numel(Agentmat);
%         rorder3 = randperm(sizegrid);
%         thecart = rorder3(1:cartnum);
%         Agentmat(thecart) = 4;
%         [Icart,Jcart,Kcart] = ind2sub(size(Agentmat),thecart);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         theinds = [Icart;Jcart;Kcart];
%         XYZcart = theinds';
%         Statecart = 4*ones(cartnum,1);%unknown if needed
%     end

%Neighborhood is 3by3by3
[row,col,len] = ind2sub([3 3 3] ,1:27);
indmat = [row',col',len'];

%Setup Angio Model

%Choose VEGF mode
VEGFmode = 'cells';%'constant';%'gradientperp';%'gradient';%'tumorsphere';%'constant';
c = clock;

%Choose capillary mode %V1.2
%'twocap' or 'aminas' or 'shortcap'
Capmode = 'endcap';

%Choose parameters
sigma = pi/6;%(pi/4)*(tipVEGF.^3)./(.53^3+tipVEGF.^3); %Sensitivity to persistence (smaller is more persistent)
sigma_VEGF = .01; %St dev of local gradient: smaller = smaller changes are more significant, ie segment more sensitive to perturbations
gamma = 1; %Weight of local factors in score
sigma_g = pi/12; %St dev of the global gradient: smaller = more strictly enforced
C_gradient = 0; %Weight of global gradient in score %1 means the cells move toward corner
global vegf_constant;% Made global for use in SegmentS2 MV HF
vegf_constant = .1; %the standard concentration of vegf HF MV

%Minimum length
min_length = 2; %Was 5 should be 2

%Choose time length
maxtime = 240;

%Do VEGF degradation?
VEGFdeg = 0; %1 = Yes, 0 = No

tempseglist = [];
tempseglist1=tempseglist;
proloop1=12;
migloop1=1; %should be between 0.5 and 2
runloop1=1;
Tiplog1 =[];
p1=loop;
prolifval = 12; %= proloop1(p1);
migmulti= BB(loop);%2;%migloop1(p1);
%runnum = runloop1(p1);
dist_vec={};
dist_vec1=dist_vec;

%Version 18
%Last revision: Will Date 2/10/12
%
%               Kerri Oct 21 2013
%               Kerri Sept 21 2015

%This one includes salt and pepper patterning for branching V5
%Fixed the hit2 section so it is more clear - No Bug
%Uses fillAgentGridTortuosity
%Includes a new anastomosis
%Includes migration that fills grid each time
%This includes the changes to correspond to Will's code
%The stalk cells proliferate here
%Tip cells can only proliferate if they are the first sprout
%removing activelist and fixing checkAnast
%Here there is a new capillary setup with the inital setup moved
%I also made branching limited to the quiescent cells (not stalk or tip)
%Changing migration so that is always adds to previous length
%Tip cell splits in half when stalk proliferates
% Add regression to the mix
%Need to make vessels that hit the wall be mature

%Data variables
branchpossible = 0;
branchcount = 0;
%TipsperCap = 0;
pixelsize = 2; %in microns BSRI 19 made it 2 instead of 1 for hardware limitations
%numspstalk = 1; %The number of stalks that can proliferate after the first = supporting stalks

%min_length = 2;
%Storage6
SegmentMatrix = zeros(1000,9); % [N1x, N1y, N1z, N2x, N2y, N2z, capnum, segnum, mature]

%Setup grid
%This sets up the initial Grid that consists of the voxels.  This is a
%structure of arrays with the different species and (possibly agents).
%Each row will contain one 'variable'.
%I got rid of Will's setup segment agents
%Changing how we set up the matrices - logical etc.

%Version 3
%Last revision: Kerri - April 16, 2012
%               Will - Feb 2, 2012
%               Kerri 3/13/13

%Params
initO2 = .02;
maxsegments = 500;
gridsize = [500,500,500];% BSRI 19 set to 500 times 2 for pixel count due to hardware limitations
% for tumorsphere
VEGFmax = 30;
VEGFdenom = 1.67;

%Initialize the grid
voxelgrid = [];
VEGFgrad = {};

%Setup the Agentgrid %ADDED BY KERRI
voxelgrid.Agent = false(gridsize(1),gridsize(2),gridsize(3));
voxelgrid.CapInd = zeros(gridsize(1),gridsize(2),gridsize(3), 'int16');
voxelgrid.PosInd = zeros(gridsize(1),gridsize(2),gridsize(3),'int16');
%voxelgrid.O2 = false(gridsize(1),gridsize(2),gridsize(3));
vector = zeros(5*5*5,3);
count = 1;
for x = -2:2
    for y=-2:2
        for z = -2:2
            vector(count,:) = [x y z];
            count = count+1;
        end
    end
end
vector2 = zeros(3*3*3,3);
count = 1;
for x=-1:1
    for y=-1:1
        for z=-1:1
            vector2(count,:) = [x y z];
            count = count+1;
        end
    end
end

voxelgrid.Check = setdiff(vector,vector2,'rows');  %This finds the values in vector not in vector2

ActiveList = zeros(maxsegments,1);

%% Vegf initialization CHOOSE MODE 1=Amina's 0=flat gradient w/ our test parameters
%VEGFmode = 'constant'; %'gradient','aminas', 'constant' %Kerri Now in
%Control
MatureNodes = ones(1,3);
%Setup Capillary
C1 = capillary;
C2 = capillary;
C3 = capillary;
C4 = capillary;
C5 = capillary;
C6 = capillary;
C7 = capillary;
C8 = capillary;
%Input Initial Capillaries
%THESE VALUES ARE ALL DIVIDED BY 2 FROM ORIGIONAL, due to a change in
%pixelsize from 1 to 2.
%We chose to make the vasculature fit the grid rather than chop off the
%parts that go outside the boundaries BSRI 19
capillary1 = {[115,15,10],[113,14,25];[113,14,25],[108,13,43];[108,13,43],[99,13,68];[99,13,68],[90,15,90];[90,15,90],[76,14,123];[76,14,123],[65,15,145];[65,15,145],[49,15,173];[49,15,173],[36,15,195];[36,15,195],[48,15,230];[48,15,230],[64,15,250];[64,15,250],[93,21,310];[93,21,310],[97,21,380];[97,21,380],[112,21,410];[112,21,410],[117,21,490];[117,21,490],[133,22,650];[133,22,650],[138,22,710];[138,22,710],[134,22,815];[134,22,815],[138,22,875];[138,22,875],[140,22,910];[140,24,910],[143,24,980]};
capillary2 = {[483,20,10],[493,21,20];[493,21,20],[497,21,40];[497,21,40],[512,21,60];[512,21,60],[517,21,80];[517,21,80],[533,22,100];[533,22,100],[538,22,120];[538,22,120],[534,22,190];[534,22,190],[538,22,210];[538,22,210],[540,23,280];[540,23,280],[543,24,350];[543,24,350],[543,21,420];[543,21,420],[547,21,490];[547,21,490],[562,21,575];[562,21,575],[567,21,630];[567,21,630],[583,22,700];[583,22,700],[588,22,770];[588,22,770],[584,22,800];[584,22,800],[588,22,865];[588,22,865],[590,22,930];[590,24,930],[593,24,970]};
capillary3 = {[965,15,10],[958,14,40];[958,14,40],[952,13,55];[952,13,55],[941,13,90];[941,13,90],[938,15,120];[938,15,120],[935,14,145];[935,14,145],[930,15,170];[930,15,170],[941,15,240];[941,15,240],[949,15,310];[949,15,310],[958,15,380];[958,15,380],[964,15,415];[964,15,415],[973,14,490];[973,14,490],[988,13,520];[988,13,520],[999,13,560];[999,13,560],[990,15,630];[990,15,630],[976,14,650];[976,14,650],[965,15,730];[965,15,730],[984,15,760];[984,15,760],[986,15,845];[986,15,845],[998,15,930];[998,15,930],[991,15,990]};
capillary4 = {[97,21,380],[88,21,370];[88,21,370],[95,21,317];[95,21,317],[103,21,260];[103,21,260],[121,21,253];[121,21,253],[134,22,235];[134,22,235],[152,24,213];[152,24,213],[159,27,190];[159,27,190],[163,25,161];[163,25,161],[176,22,130];[176,22,130],[183,24,127];[183,24,127],[188,21,120];[188,21,120],[195,21,117];[195,21,117],[203,21,110];[203,21,110],[221,21,103];[221,21,103],[234,22,90];[234,22,90],[252,24,78];[252,24,78],[259,27,40];[259,27,40],[263,25,31];[263,25,31],[276,22,28];[276,22,28],[283,24,13]};
capillary5 = {[973,14,490],[983,14,482];[983,14,482],[960,13,480];[960,13,480],[936,15,445];[936,15,445],[913,15,418];[913,15,418],[886,19,400];[886,19,400],[870,20,327];[870,20,327],[864,24,290];[864,24,290],[856,27,256];[856,27,256],[843,25,240];[843,25,240],[831,22,217];[831,22,217],[793,19,197];[793,19,197],[770,18,180];[770,18,180],[751,15,165];[751,15,165],[743,10,128];[743,10,128],[726,9,120];[726,9,120],[715,5,112];[715,5,112],[714,9,95];[714,9,95],[701,12,86];[701,12,86],[693,15,72];[693,15,72],[665,17,67]};
%capillary6 = {[100,300,50],[100,300,100];[100,300,50],[100,300,100];};
capillary6 = {[15,115,10],[14,113,25];[14,113,25],[13,108,43];[13,108,43],[13,99,68];[13,99,68],[15,90,90];[15,90,90],[14,76,123];[14,76,123],[15,65,145];[15,65,145],[15,49,173];[15,49,173],[15,36,195];[15,36,195],[15,48,230];[15,48,230],[15,64,250];[15,64,250],[21,93,310];[21,93,310],[21,97,380];[21,97,380],[21,112,410];[21,112,410],[21,117,490];[21,117,490],[22,133,650];[22,133,650],[22,138,710];[22,138,710],[22,134,815];[22,134,815],[22,138,875];[22,138,875],[22,140,910];[24,140,910],[24,143,980]};
capillary7 = {[20,483,10],[21,493,20];[21,493,20],[21,497,40];[21,497,40],[21,512,60];[21,512,60],[21,517,80];[21,517,80],[22,533,100];[22,533,100],[22,538,120];[22,538,120],[22,534,190];[22,534,190],[22,538,210];[22,538,210],[23,540,280];[23,540,280],[24,543,350];[24,543,350],[21,543,420];[21,543,420],[21,547,490];[21,547,490],[21,562,575];[21,562,575],[21,567,630];[21,567,630],[22,583,700];[22,583,700],[22,588,770];[22,588,770],[22,584,800];[22,584,800],[22,588,865];[22,588,865],[22,590,930];[24,590,930],[24,593,970]};
capillary8 = {[15,965,10],[14,958,40];[14,958,40],[13,952,55];[13,952,55],[13,941,90];[13,941,90],[15,938,120];[15,938,120],[14,935,145];[14,935,145],[15,930,170];[15,930,170],[15,941,240];[15,941,240],[15,949,310];[15,949,310],[15,958,380];[15,958,380],[15,964,415];[15,964,415],[14,973,490];[14,973,490],[13,988,520];[13,988,520],[13,999,560];[13,999,560],[15,990,630];[15,990,630],[14,976,650];[14,976,650],[15,965,730];[15,965,730],[15,984,760];[15,984,760],[15,986,845];[15,986,845],[15,998,930];[15,998,930],[15,991,990]};

for ik=1:length(capillary1)
    newSeg = segmentS2;
    newSeg.Node1 = capillary1{ik,1};
    newSeg.Node2 = capillary1{ik,2};
    MatureNodes = cat(1,MatureNodes,newSeg.Node1);
    MatureNodes = cat(1,MatureNodes,newSeg.Node2);
    newSeg.Radius = 5;
    newSeg.Listpos = ik;
    newSeg.CapillaryNum = 1;
    newSeg.updateDirection();

    %fill in SegmentMatrix
    SegmentMatrix(ik,1:3) =  capillary1{ik,1};
    SegmentMatrix(ik,4:6) =  capillary1{ik,2};
    SegmentMatrix(ik,7) =  1;
    SegmentMatrix(ik,8) =  ik;
    SegmentMatrix(ik,9) =  1;
    if ik==length(capillary1)% To prevent spidering BSRI 19
        newSeg.Branched = 1;
    end

    C1.SegmentList{ik} = newSeg;

    %Add new cell location to Agent List %ADDED BY KERRI V1.2
    %Find new cells nodes and radius
    Node1pt = newSeg.Node1;
    Node2pt = newSeg.Node2;
    segradius = newSeg.Radius;
    %Use FillAgentGrid %NEED TO FIX THIS
    [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);
end

for j=1:length(capillary2)
    newSeg = segmentS2;
    newSeg.Node1 = capillary2{j,1};
    newSeg.Node2 = capillary2{j,2};
    MatureNodes = cat(1,MatureNodes,newSeg.Node1);
    MatureNodes = cat(1,MatureNodes,newSeg.Node2);
    newSeg.Radius = 5;
    newSeg.Listpos = j;
    newSeg.CapillaryNum = 2;
    newSeg.updateDirection();
    jj2 = ik+j;
    %fill in SegmentMatrix
    SegmentMatrix(jj2,1:3) =  capillary2{j,1};
    SegmentMatrix(jj2,4:6) =  capillary2{j,2};
    SegmentMatrix(jj2,7) =  1;
    SegmentMatrix(jj2,8) =  j;
    SegmentMatrix(jj2,9) =  1;
    if j==length(capillary2)
        newSeg.Branched = 1;
    end
    %Setup Capillary
    C2.SegmentList{j} = newSeg;

    % Add new cell location to Agent List %ADDED BY KERRI V1.2
    %Find new cells nodes and radius
    Node1pt = newSeg.Node1;
    Node2pt = newSeg.Node2;
    segradius = newSeg.Radius;
    %Use FillAgentGrid %NEED TO FIX THIS
    [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);
end

for jk=1:length(capillary3)
    newSeg = segmentS2;
    newSeg.Node1 = capillary3{jk,1};
    newSeg.Node2 = capillary3{jk,2};
    MatureNodes = cat(1,MatureNodes,newSeg.Node1);
    MatureNodes = cat(1,MatureNodes,newSeg.Node2);
    newSeg.Radius = 5;
    newSeg.Listpos = jk;
    newSeg.CapillaryNum = 3;
    newSeg.updateDirection();
    jj3 = jk+jj2;
    %fill in SegmentMatrix
    SegmentMatrix(jj3,1:3) =  capillary3{jk,1};
    SegmentMatrix(jj3,4:6) =  capillary3{jk,2};
    SegmentMatrix(jj3,7) =  1;
    SegmentMatrix(jj3,8) =  jk;
    SegmentMatrix(jj3,9) =  1;
    if jk==length(capillary3)
        newSeg.Branched = 1;
    end
    %Setup Capillary
    C3.SegmentList{jk} = newSeg;


    % Add new cell location to Agent List %ADDED BY KERRI V1.2
    %Find new cells nodes and radius
    Node1pt = newSeg.Node1;
    Node2pt = newSeg.Node2;
    segradius = newSeg.Radius;
    %Use FillAgentGrid %NEED TO FIX THIS
    [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);
end

for jh=1:length(capillary4)
    newSeg = segmentS2;
    newSeg.Node1 = capillary4{jh,1};
    newSeg.Node2 = capillary4{jh,2};
    MatureNodes = cat(1,MatureNodes,newSeg.Node1);
    MatureNodes = cat(1,MatureNodes,newSeg.Node2);
    newSeg.Radius = 5;
    newSeg.Listpos = jh;
    newSeg.CapillaryNum = 4;
    newSeg.updateDirection();
    %seg_list{i} = newSeg;
    jj4 = jh+jj3;
    %fill in SegmentMatrix
    SegmentMatrix(jj4,1:3) =  capillary4{jh,1};
    SegmentMatrix(jj4,4:6) =  capillary4{jh,2};
    SegmentMatrix(jj4,7) =  1;
    SegmentMatrix(jj4,8) =  jh;
    SegmentMatrix(jj4,9) =  1;
    if jh==length(capillary4)
        newSeg.Branched = 1;
    end
    %Setup Capillary
    %C2 = capillary;
    C4.SegmentList{jh} = newSeg;

    % Add new cell location to Agent List %ADDED BY KERRI V1.2
    %Find new cells nodes and radius
    Node1pt = newSeg.Node1;
    Node2pt = newSeg.Node2;
    segradius = newSeg.Radius;
    %Use FillAgentGrid %NEED TO FIX THIS
    [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);
    %Fill O2 grid
    %[voxelgrid] = fillAgentGridO2(voxelgrid, newSeg, pixelsize);


    %SegmentList{i} = newSeg;
    %obj.Agents.AssignSpecies(newSeg.Node1(1),newSeg.Node1(2),newSeg.Node1(3),'Segment',newSeg);
    %obj.Agents.AssignSpecies(newSeg.Node2(1),newSeg.Node2(2),newSeg.Node2(3),'Segment',newSeg);
end

for hh=1:length(capillary5)
    newSeg = segmentS2;
    newSeg.Node1 = capillary5{hh,1};
    newSeg.Node2 = capillary5{hh,2};
    MatureNodes = cat(1,MatureNodes,newSeg.Node1);
    MatureNodes = cat(1,MatureNodes,newSeg.Node2);

    newSeg.Radius = 5;
    newSeg.Listpos = hh;
    newSeg.CapillaryNum = 5;
    newSeg.updateDirection();
    %seg_list{i} = newSeg;
    jj5 = hh+jj4;
    %fill in SegmentMatrix
    SegmentMatrix(jj5,1:3) =  capillary5{hh,1};
    SegmentMatrix(jj5,4:6) =  capillary5{hh,2};
    SegmentMatrix(jj5,7) =  1;
    SegmentMatrix(jj5,8) =  hh;
    SegmentMatrix(jj5,9) =  1;
    if hh==length(capillary5)
        newSeg.Branched = 1;
    end
    %Setup Capillary
    %C2 = capillary;
    C5.SegmentList{hh} = newSeg;

    % Add new cell location to Agent List %ADDED BY KERRI V1.2
    %Find new cells nodes and radius
    Node1pt = newSeg.Node1;
    Node2pt = newSeg.Node2;
    segradius = newSeg.Radius;
    %Use FillAgentGrid %NEED TO FIX THIS
    [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);
    %Fill O2 grid
    % [voxelgrid] = fillAgentGridO2(voxelgrid, newSeg, pixelsize);


    %SegmentList{i} = newSeg;
    %obj.Agents.AssignSpecies(newSeg.Node1(1),newSeg.Node1(2),newSeg.Node1(3),'Segment',newSeg);
    %obj.Agents.AssignSpecies(newSeg.Node2(1),newSeg.Node2(2),newSeg.Node2(3),'Segment',newSeg);
end

for hj=1:length(capillary6)
    newSeg = segmentS2;
    newSeg.Node1 = capillary6{hj,1};
    newSeg.Node2 = capillary6{hj,2};
    MatureNodes = cat(1,MatureNodes,newSeg.Node1);
    MatureNodes = cat(1,MatureNodes,newSeg.Node2);

    newSeg.Radius = 5;
    newSeg.Listpos = hj;
    newSeg.CapillaryNum = 6;
    newSeg.updateDirection();
    %seg_list{i} = newSeg;
    jj6 = jj5+hj;
    %fill in SegmentMatrix
    SegmentMatrix(jj6,1:3) =  capillary6{hj,1};
    SegmentMatrix(jj6,4:6) =  capillary6{hj,2};
    SegmentMatrix(jj6,7) =  1;
    SegmentMatrix(jj6,8) =  hj;
    SegmentMatrix(jj6,9) =  1;
    if hj==length(capillary6)
        newSeg.Branched = 1;
    end
    %Setup Capillary
    %C2 = capillary;
    C6.SegmentList{hj} = newSeg;

    % Add new cell location to Agent List %ADDED BY KERRI V1.2
    %Find new cells nodes and radius
    Node1pt = newSeg.Node1;
    Node2pt = newSeg.Node2;
    segradius = newSeg.Radius;
    %Use FillAgentGrid %NEED TO FIX THIS
    [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);
    %Fill O2 grid
    % [voxelgrid] = fillAgentGridO2(voxelgrid, newSeg, pixelsize);


    %SegmentList{i} = newSeg;
    %obj.Agents.AssignSpecies(newSeg.Node1(1),newSeg.Node1(2),newSeg.Node1(3),'Segment',newSeg);
    %obj.Agents.AssignSpecies(newSeg.Node2(1),newSeg.Node2(2),newSeg.Node2(3),'Segment',newSeg);
end

for hhh=1:length(capillary7)
    newSeg = segmentS2;
    newSeg.Node1 = capillary7{hhh,1};
    newSeg.Node2 = capillary7{hhh,2};
    MatureNodes = cat(1,MatureNodes,newSeg.Node1);
    MatureNodes = cat(1,MatureNodes,newSeg.Node2);

    newSeg.Radius = 5;
    newSeg.Listpos = hhh;
    newSeg.CapillaryNum = 7;
    newSeg.updateDirection();
    %seg_list{i} = newSeg;
    jj7 = jj6+hhh;
    %fill in SegmentMatrix
    SegmentMatrix(jj7,1:3) =  capillary7{hhh,1};
    SegmentMatrix(jj7,4:6) =  capillary7{hhh,2};
    SegmentMatrix(jj7,7) =  1;
    SegmentMatrix(jj7,8) =  hhh;
    SegmentMatrix(jj7,9) =  1;
    if hhh==length(capillary7)
        newSeg.Branched = 1;
    end
    %Setup Capillary
    %C2 = capillary;
    C7.SegmentList{hhh} = newSeg;

    % Add new cell location to Agent List %ADDED BY KERRI V1.2
    %Find new cells nodes and radius
    Node1pt = newSeg.Node1;
    Node2pt = newSeg.Node2;
    segradius = newSeg.Radius;
    %Use FillAgentGrid %NEED TO FIX THIS
    [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);
    %Fill O2 grid
    % [voxelgrid] = fillAgentGridO2(voxelgrid, newSeg, pixelsize);


    %SegmentList{i} = newSeg;
    %obj.Agents.AssignSpecies(newSeg.Node1(1),newSeg.Node1(2),newSeg.Node1(3),'Segment',newSeg);
    %obj.Agents.AssignSpecies(newSeg.Node2(1),newSeg.Node2(2),newSeg.Node2(3),'Segment',newSeg);
end

for hw=1:length(capillary8)
    newSeg = segmentS2;
    newSeg.Node1 = capillary8{hw,1};
    newSeg.Node2 = capillary8{hw,2};
    MatureNodes = cat(1,MatureNodes,newSeg.Node1);
    MatureNodes = cat(1,MatureNodes,newSeg.Node2);

    newSeg.Radius = 5;
    newSeg.Listpos = hw;
    newSeg.CapillaryNum = 8;
    newSeg.updateDirection();
    %seg_list{i} = newSeg;
    jj8 = hw+jj7;
    %fill in SegmentMatrix
    SegmentMatrix(jj8,1:3) =  capillary8{hw,1};
    SegmentMatrix(jj8,4:6) =  capillary8{hw,2};
    SegmentMatrix(jj8,7) =  1;
    SegmentMatrix(jj8,8) =  hw;
    SegmentMatrix(jj8,9) =  1;
    if hw==length(capillary8)
        newSeg.Branched = 1;
    end
    %Setup Capillary
    %C2 = capillary;
    C8.SegmentList{hw} = newSeg;

    % Add new cell location to Agent List %ADDED BY KERRI V1.2
    %Find new cells nodes and radius
    Node1pt = newSeg.Node1;
    Node2pt = newSeg.Node2;
    segradius = newSeg.Radius;
    %Use FillAgentGrid %NEED TO FIX THIS
    [voxelgrid] = fillAgentGrid2(voxelgrid, newSeg, pixelsize);
    %Fill O2 grid
    % [voxelgrid] = fillAgentGridO2(voxelgrid, newSeg, pixelsize);


    %SegmentList{i} = newSeg;
    %obj.Agents.AssignSpecies(newSeg.Node1(1),newSeg.Node1(2),newSeg.Node1(3),'Segment',newSeg);
    %obj.Agents.AssignSpecies(newSeg.Node2(1),newSeg.Node2(2),newSeg.Node2(3),'Segment',newSeg);
end

%Find nonempty entries of CapillaryList

CapillaryList = cell(maxsegments,1);
CapillaryList{1} = C1;
CapillaryList{2} = C2;
CapillaryList{3} = C3;
CapillaryList{4} = C4;
CapillaryList{5} = C5;
CapillaryList{6} = C6;
CapillaryList{7} = C7;
CapillaryList{8} = C8;
timecap = 1;

%List for all the endothelial cell start positions
VVDlist = zeros(2000);

% Active a node on each capillary
fprintf('initial!\n');
caps = CapillaryList(~cellfun('isempty',CapillaryList));

for ii=1:8%Must be 1: size of CapillaryList

    nonempty = find(~cellfun('isempty',CapillaryList{ii}.SegmentList));
    %Make 2 initial branches
    %initialTip=randi(length(nonempty));
    %All segments in capillary network branch
    for cc = 1:length(nonempty)
        S1 = CapillaryList{ii}.SegmentList{cc};

        %% BRANCH    WILL CHANGED v2.3--BRANCHES NOW
        %display('First branches')
        Hypbin = HypCell == 1;
        XYZh = XYZ(Hypbin,:);
        [CapillaryList, SegmentMatrix,voxelgrid] = branch2(CapillaryList,SegmentMatrix,S1,Agentmat,voxelgrid, XYZh, prolifval, pixelsize);
        CapillaryList{ii}.SegmentList{cc}.Branched = 1;
    end
end %end if available

%% Start Iterations
fprintf('It starts\n');

for time = 1:n_loops                                                                                                                                                                 %For now 50 iteration, 3600 about 10 years

    %  hyploc = zeros(3600,3);
    hypcount = 1;
    fprintf("Starting loop %d with %d cells.\n",time,numberofcells)

    % numberofcells = length(theCells);
    possiblydead = []; %Keep track of the possibly dead
    possiblydead2 = [];

    %% Fibroblast module - they are quicker so they go first
    if size(XYZf) > 1
        numberofcellsf =size(XYZf,1);
        %randperm
        randcellf = randperm(numberofcellsf);
        jf = 0;

        %Go through all of the current cells
        for jj = 1:numberofcellsf %proliferate all cells so far but not the new ones
            %Don't update the numberofcells
            %randomly choose the cells
            jf = randcellf(jj);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%             Migration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Version 2
            %Check if there is empty space
            %% Check if there is space
            %Note: jm is the random fibroblast position
            [emptypos1] = FindEmptySpace(Agentmat,XYZf(jf,:), indmat);

            if ~isempty(emptypos1)

                Fmigrate = 1;
                %Check how many cell spaces to move
                migmoves = Fmigrate - 1;

                %Record old position
                oldpos = XYZf(jf,:);

                %Update new position
                XYZf(jf,:) = emptypos1;

                %for loop for number of moves
                if migmoves > 0
                    for mi = 1:migmoves
                        %find next empty spot
                        emptypos2 = FindEmptySpace(Agentmat,XYZf(jf,:), indmat);
                        if ~isempty(emptypos2)
                            %Update new position
                            XYZf(jf,:) = emptypos2;
                        end
                    end %mi = 1:migmoves
                end %if migmoves >0

                %Remove previous position from Agentmat
                Agentmat(oldpos(1), oldpos(2), oldpos(3)) = 0;

                %Put new position in Agentmat - macrophage = 2
                newpos = XYZf(jf,:);
                Agentmat(newpos(1), newpos(2), newpos(3)) = 3;
            end
        end
    end


    %% Macrophage module - they are not random so they go last
    if length(XYZm) > 1
        numberofcellsm =length(XYZm);
        %randperm
        randcellm = randperm(numberofcellsm);
        jm = 0;

        %Go through all of the current cells
        for jj = 1:numberofcellsm %proliferate all cells so far but not the new ones
            %Don't update the numberofcells
            %randomly choose the cells
            jm = randcellm(jj);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%             Migration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Version 2
            %Check if there is empty space
            %% Check if there is space
            %Note: jm is the random macrophage position
            [emptypos] = FindEmptySpace(Agentmat,XYZm(jm,:), indmat);

            if ~isempty(emptypos)

                %% Macrophage Migration

                %Check migration cycles are the same - now migration happens
                %every iteration
                %Test distance from tumor cells
                mpos =  XYZm(jm,:);
                thedists = pdist2(mpos, XYZ);
                [mindist, minpos] = min(thedists);

                %If close to cell Mmigrate = 8 otherwise = 1

                if mindist < 10 %Should really be 10
                    Mmigrate = 8;
                    Mrun = 8;

                    %xm = sqrt((mpos(1))^2-(thedists(minpos,1))^2);
                    r1 = mpos;
                    r2 = XYZ(minpos,:);
                    onfirst = r1 == 1;
                    onlast = size(Agentmat) == r1;
                    if sum(onfirst >0) || sum(onlast > 0)
                        %On the edge - don't move
                    elseif r1 ~= r2
                        %not on edge
                        %find the unit vector
                        v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector; = unit vector
                        %if free - assign emptypos
                        testpos = r1 + round(v);
                        newpos = r1; %just in case it can't move towards nearest cell
                        %Test if the position is free
                        isfree = Agentmat(testpos(1), testpos(2), testpos(3));

                        while isfree == 0 && Mrun > 0 %If the cell can go towards closest Tumor cell
                            newpos = testpos;
                            onfirst = newpos == 1;
                            onlast = size(Agentmat) == newpos;
                            if sum(onfirst >0) || sum(onlast > 0)
                                %On the edge - don't move
                                Mrun = 0;
                            else
                                Mrun = Mrun -1;
                                testpos = newpos + round(v);
                                %Test if the position is free
                                isfree = Agentmat(testpos(1), testpos(2), testpos(3));
                            end
                        end

                        if Mrun > 0 %If cell can't go toward the tumor cell
                            for mm = 1:Mrun
                                %Randomly move
                                [emptypos2] = FindEmptySpace(Agentmat,newpos, indmat);
                                newpos = emptypos2;
                                XYZm(jm,:) = emptypos2;
                            end
                        else
                            emptypos2 = testpos;
                            XYZm(jm,:) = emptypos2;
                        end

                        %Reassign new cell
                        Agentmat(emptypos2(1),emptypos2(2),emptypos2(3)) = 2;
                        %Remove previous position from Agentmat
                        Agentmat(mpos(1), mpos(2), mpos(3)) = 0;
                    end
                else
                    Mmigrate = 1;
                    %Check how many cell spaces to move
                    migmoves = Mmigrate - 1;

                    %Record old position
                    oldpos = XYZm(jm,:);

                    %Update new position
                    XYZm(jm,:) = emptypos;

                    %for loop for number of moves
                    if migmoves > 0
                        for mi = 1:migmoves
                            %find next empty spot
                            emptypos2 = FindEmptySpace(Agentmat,XYZm(jm,:), indmat);
                            if ~isempty(emptypos2)
                                %Update new position
                                XYZm(jm,:) = emptypos2;
                            end
                        end %mi = 1:migmoves
                    end %if migmoves >0

                    %Remove previous position from Agentmat
                    Agentmat(oldpos(1), oldpos(2), oldpos(3)) = 0;

                    %Put new position in Agentmat - macrophage = 2
                    newpos = XYZm(jm,:);
                    Agentmat(newpos(1), newpos(2), newpos(3)) = 2;
                end
            end
        end
    end
    %% Angiogenesis Module

    %Run Angio
    %Test if there is a hypoxic cell
    hyptest = sum(HypCell);
    %         if 1<= hyptest && hyptest <=5
    %             disp('first hypoxic cell')
    %         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Make this happen 6 times
    if hyptest > 0
        %             disp('angiogenesis')
        %Each iteration is 6 hours so run 6 times
        % Run for many more (10) times and see whether it still grows
        for hrsix = 1:3
            %Housekeeping
            %clear SaveAddress p SaveFig
            close all
            %Find nonempty entries of CapillaryList
            caps = CapillaryList(~cellfun('isempty',CapillaryList));

            %Loop through all capillaries in the list - except initial ones
            for ii=9:length(caps)
                %Find each segment
                thesegs = CapillaryList{ii}.SegmentList(~cellfun('isempty',CapillaryList{ii}.SegmentList));
                tipnum = length(thesegs);

                %% Tip cell Actions

                %Check if it is still active
                [binarya] = CapillaryList{ii}.SegmentList{tipnum}.canActivate(Agentmat, pixelsize);
                [binarym] = CapillaryList{ii}.SegmentList{tipnum}.Mature;
                %Only do if there are active ones and not mature
                if binarya == 1 && binarym == 0
                    % error('After check active')
                    %If the tip cell has reached

                    if CapillaryList{ii}.SegmentList{tipnum}.CellClock >=  CapillaryList{ii}.SegmentList{tipnum}.CycleLength && CapillaryList{ii}.SegmentList{tipnum}.Listpos == 1
                        vector = round(pixelsize*CapillaryList{ii}.SegmentList{tipnum}.direction);
                        %Test whether the direction is valid -ADDED
                        %V4.2 KERRI
                        tipPosition = CapillaryList{ii}.SegmentList{tipnum}.Node2;
                        tipPosition = round(tipPosition./pixelsize);
                        startPoint = [tipPosition(1)-1, tipPosition(2)-1, tipPosition(3)-1]; %start search at bottom left corner of "cube" around tip

                        %Need to check that the startpoint and startpoint + 2 are not
                        %outside of the grid
                        gridsize = size(voxelgrid.Agent);
                        binary1 = checkBoundaries3(startPoint, gridsize);
                        binary2 = checkBoundaries3(startPoint+2, gridsize);

                        if binary1 == 0 && binary2 == 0
                            %ADDED KERRI V4.2 %The search must have been allowed (node not on edge)

                            %Make the segments longer - V4 DELETED
                            %KERRI

                            %Check whether the new segment will leave the grid
                            newpoint = CapillaryList{ii}.SegmentList{tipnum}.Node2 + vector;
                            sizegrid = size(voxelgrid.Agent);
                            binary = checkBoundaries3(round(newpoint/pixelsize), sizegrid);

                            %Do Prolif if it doesn't leave grid
                            if binary == 0 %Prolif if it doesn't leave the grid

                                %Check whether it can prolif %ADDED BY KERRI V3.1
                                [binaryP] = CapillaryList{ii}.SegmentList{tipnum}.canProliferateTip;%Added Tip


                                if binaryP == 1 %Check it can prolif
                                    % display('proliferating tip!') %REMOVED vector V4.2 KERRI
                                    Listpos = tipnum+1;

                                    %Still need to input tip cell here[newSeg] = proliferate(obj, Listposval, prolifval, pixelsize)
                                    [newcell,newstalk,evaldum] = CapillaryList{ii}.SegmentList{tipnum}.proliferate(prolifval, pixelsize);%KAN 3/12/13
                                    %CapillaryList{i}.ActiveList = ActiveList;

                                    if newcell.Node1 == newcell.Node2
                                        error('Nodes are the same')
                                    end

                                    newstalk.CellClock = 0;
                                    CapillaryList{ii}.SegmentList{tipnum} = newstalk;

                                    uniquepos = find(SegmentMatrix(:,7) == ii & SegmentMatrix(:,8) == tipnum);

                                    %fill in SegmentMatrix
                                    SegmentMatrix(uniquepos,1:3) =  newstalk.Node1;
                                    SegmentMatrix(uniquepos,4:6) =  newstalk.Node2;
                                    SegmentMatrix(uniquepos,7) =  newstalk.CapillaryNum;
                                    SegmentMatrix(uniquepos,8) = newstalk.Listpos;
                                    SegmentMatrix(uniquepos,9) =  newstalk.Mature;

                                    %Add newcell to SegmentList
                                    CapillaryList{ii}.SegmentList{Listpos} = newcell;

                                    %fill in SegmentMatrix
                                    epos = find(SegmentMatrix(:,1) == 0, 1 );
                                    if isempty(epos)
                                        epos = length(SegmentMatrix(:,1)) +1;
                                    end
                                    SegmentMatrix(epos,1:3) =  newcell.Node1;
                                    SegmentMatrix(epos,4:6) =  newcell.Node2;
                                    SegmentMatrix(epos,7) =  newcell.CapillaryNum;
                                    SegmentMatrix(epos,8) = newcell.Listpos;
                                    SegmentMatrix(epos,9) =  newcell.Mature;
                                    if newcell.Mature == 1
                                        MatureNodes = cat(1,MatureNodes,newcell.Node1);
                                        MatureNodes = cat(1,MatureNodes,newcell.Node2);
                                    end
                                    %                            display('prolif len')

                                    %% Add new cell location to Agent List %ADDED BY KERRI V1.2
                                    %Find new cells nodes and radius
                                    %                                              Node1pt = newcell.Node1;
                                    %                                              Node2pt = newcell.Node2;
                                    %                                              segradius = newcell.Radius;

                                    %% Test Anastomosis %Added V8
                                    %Check Anast fills in the grid
                                    if newcell.Listpos == 1
                                        %This is the first segment
                                        [doanast, gridindex] = CheckAnastBranch(voxelgrid, newcell, CapillaryList, pixelsize);
                                    else
                                        %Previous segment
                                        oldseg =  CapillaryList{ii}.SegmentList{tipnum};
                                        [doanast, gridindex] = CheckAnast(voxelgrid, newcell, oldseg, CapillaryList, pixelsize);
                                    end

                                    if doanast == 1
                                        %Do anastomosis
                                        [voxelgrid, SegmentMatrix, newseg, CapillaryList, CapillaryList{ii}.ActiveList,MatureNodes] = doAnastSM(newcell, voxelgrid, SegmentMatrix, gridindex, CapillaryList, CapillaryList{ii}.ActiveList, pixelsize,MatureNodes);
                                        %Update new cell
                                        newseg.updateDirection();
                                        CapillaryList{ii}.SegmentList{Listpos} = newseg;
                                        %Update new positions
                                        SegmentMatrix(epos,4:6) =  newseg.Node2;

                                    else
                                        %Fill in the grid
                                        [voxelgrid,CapillaryList{ii}.SegmentList{Listpos}] = fillAgentGridTortNodes(voxelgrid, CapillaryList{ii}.SegmentList{Listpos}, pixelsize);
                                    end
                                end % if binaryP == 1
                            else
                                %It would have left the grid
                                %deactivate the segment
                                %Make boundary capillaries mature
                                Pos = tipnum;
                                for cpos = 1:Pos
                                    CapillaryList{ii}.SegmentList{cpos}.Mature = 1;
                                    %Record mature nodes
                                    MatureNodes = cat(1,MatureNodes,CapillaryList{ii}.SegmentList{cpos}.Node1);
                                    MatureNodes = cat(1,MatureNodes,CapillaryList{ii}.SegmentList{cpos}.Node2);
                                end
                            end
                        else %if a node is on the edge of the grid
                            %deactivate the segment
                            %                      ActiveList = CapillaryList{i}.SegmentList{theactive(j)}.deactivate(CapillaryList{i}.ActiveList);
                            %                      CapillaryList{i}.ActiveList = ActiveList;
                            Pos = tipnum;
                            for cpos = 1:Pos
                                CapillaryList{ii}.SegmentList{cpos}.Mature = 1;
                                % SegmentMatrix(uniquepos,9) =  newstalk.Mature;
                                %Record mature nodes
                                MatureNodes = cat(1,MatureNodes,CapillaryList{ii}.SegmentList{cpos}.Node1);
                                MatureNodes = cat(1,MatureNodes,CapillaryList{ii}.SegmentList{cpos}.Node2);
                            end
                        end
                    else
                        %Add to the cell cycle
                        newcycle = CapillaryList{ii}.SegmentList{tipnum}.CellClock + 1;
                        CapillaryList{ii}.SegmentList{tipnum}.CellClock = newcycle;

                        vector = CapillaryList{ii}.SegmentList{tipnum}.findDirection(Agentmat,pixelsize);

                        %Test whether the segment is NOT on the edge
                        if ~isempty(vector) %The search must have been allowed (node not on edge)
                            tempSeg = CapillaryList{ii}.SegmentList{tipnum};
                            %
                            l = norm(tempSeg.Node2 - tempSeg.Node1); %distance from Node1 to 2
                            l = max(l,1); %The length must be one micron at least
                            max_migration = 1.5*l;%From Amina's Paper: M = T2[VEGF] + migNoVegfMatrix
                            migNoVEGFMatrix = 6.2;
                            N1position = ceil((tempSeg.Node2./20));
                            theneighbs = Agentmat(max(N1position(1)-1,1):min(N1position(1)+1, size(Agentmat,1)),max(N1position(2)-1,1):min(N1position(2)+1, size(Agentmat,2)),max(N1position(3)-1,1):min(N1position(3)+1, size(Agentmat,3)))~= 0;
                            sumcells = sum(sum(sum(theneighbs)))/numel(theneighbs); %now use 30 because 2/3 of spaces are empty in tumor 20*(3/2)
                            mig_distN2 = migmulti*2*(.4*vegf_constant*(3/2)*sumcells + migNoVEGFMatrix);%6.2; %From Amina's Paper: M = dt(T2[VEGF] + migNoVegfMatrix)
                            mig_distN3 = min(min(mig_distN2, max_migration), 60);
                            mig_dist = mig_distN3 + l; %new elongation is the distance plus the old length of the segment
                            elongation = max(l, mig_dist); % Randomly chooses length of elongation for now, 1 + (.2 + 1.3*rand))

                            newpoint = tempSeg.Node2 + (vector);
                            elongation_vec =(newpoint - tempSeg.Node1)/norm(newpoint - tempSeg.Node1);
                            if isnan(elongation_vec)
                                %persistance - go the way it was originall
                                %going
                                vector =  round( (tempSeg.Node2 - tempSeg.Node1)/l);
                                newpoint = tempSeg.Node2 + (vector);
                                elongation_vec =(newpoint - tempSeg.Node1)/norm(newpoint - tempSeg.Node1);
                            end
                            sizegrid = size(voxelgrid.Agent);
                            %Changed V6.2 check boundaries of new
                            %point
                            binary = checkBoundaries3(round((tempSeg.Node1 + round(elongation_vec*elongation))/pixelsize), sizegrid);%checkBoundaries2(newpoint, sizegrid);

                            if binary == 0 %Migrate if it doesn't leave the grid
                                %didmig = 1;
                                %Remove the previous prolif segment
                                [voxelgrid] = emptyAgentGrid3(voxelgrid,CapillaryList{ii}.SegmentList{tipnum}, pixelsize);

                                sameSeg = CapillaryList{ii}.SegmentList{tipnum}.migrate(elongation_vec,elongation);
                                sameSeg.updateDirection(); %Update direction vector
                                %                                         disp('oldseg')
                                %                                           disp(CapillaryList{i}.SegmentList{theactive(j)}.Node2)
                                %                                         disp('newseg')
                                if sameSeg.Node1 == sameSeg.Node2
                                    error('Nodes are the same')
                                end

                                CapillaryList{ii}.SegmentList{tipnum} = sameSeg;
                                uniquepos = find(SegmentMatrix(:,7) == ii & SegmentMatrix(:,8) == tipnum);
                                SegmentMatrix(uniquepos,4:6) =  sameSeg.Node2;

                                %% Test Anastomosis %Added V8
                                if sameSeg.Listpos == 1
                                    %This is the first segment
                                    [doanast, gridindex] = CheckAnastBranch(voxelgrid, sameSeg, CapillaryList, pixelsize);
                                else
                                    %Previous segment
                                    oldseg =  CapillaryList{ii}.SegmentList{tipnum-1};
                                    [doanast, gridindex] = CheckAnast(voxelgrid, sameSeg, oldseg, CapillaryList, pixelsize);
                                end

                                if doanast == 1
                                    %Do anastomosis
                                    [voxelgrid, SegmentMatrix, newseg, CapillaryList, CapillaryList{ii}.ActiveList,MatureNodes] = doAnastSM(sameSeg, voxelgrid, SegmentMatrix, gridindex, CapillaryList, CapillaryList{ii}.ActiveList, pixelsize, MatureNodes);
                                    %Update new cell
                                    newseg.updateDirection();
                                    CapillaryList{ii}.SegmentList{tipnum} = newseg;
                                    uniquepos = find(SegmentMatrix(:,7) == ii & SegmentMatrix(:,8) == tipnum);
                                    SegmentMatrix(uniquepos,4:6) =  newseg.Node2;
                                else
                                    %Fill in the grid
                                    [voxelgrid,CapillaryList{ii}.SegmentList{tipnum}] = fillAgentGridTortNodes(voxelgrid, CapillaryList{ii}.SegmentList{tipnum}, pixelsize);
                                end
                            else
                            end
                        else %if a node is on the edge of the grid
                        end
                    end %Capillary
                else
                end

                %% Pushing Stalk
                pstalk = tipnum - 1;
                binaryP = 0; %Set as 0 so that it can be reset if prolif occurs
                if pstalk > 0 %make sure there is a stalk cell
                    if CapillaryList{ii}.SegmentList{pstalk}.CellClock >= CapillaryList{ii}.SegmentList{pstalk}.CycleLength
                        %Proliferate

                        %                  %find first empty space in SegmentList
                        %                  emptylist = cellfun(@(A) any(isempty(A(:))), CapillaryList{i}.SegmentList);
                        %                  Listpos = find (emptylist, 1, 'first');

                        %ADDED V4 - Prolif NOT dependent on VEGF
                        %gradient
                        %vector = CapillaryList{i}.SegmentList{theactive(j)}.findDirectionP(voxelgrid);
                        %New vector is the vector of old tip cell - but 1 in
                        %length
                        vectorlen = pixelsize*CapillaryList{ii}.SegmentList{tipnum}.direction/norm(CapillaryList{ii}.SegmentList{tipnum}.direction);

                        %Test whether the direction is valid -ADDED
                        %V4.2 KERRI
                        tipPosition = CapillaryList{ii}.SegmentList{tipnum}.Node2;
                        tipPosition = round(tipPosition./pixelsize);
                        startPoint = [tipPosition(1)-1, tipPosition(2)-1, tipPosition(3)-1]; %start search at bottom left corner of "cube" around tip

                        %Need to check that the startpoint and startpoint + 2 are not
                        %outside of the grid
                        gridsize = size(voxelgrid.Agent);
                        binary1 = checkBoundaries3(startPoint, gridsize);
                        binary2 = checkBoundaries3(startPoint+2, gridsize);

                        [binarym] = CapillaryList{ii}.SegmentList{pstalk}.Mature; %Make sure it is not mature

                        if binary1 == 0 && binary2 == 0 && binarym == 0%ADDED KERRI V4.2 %The search must have been allowed (node not on edge)

                            %Make the segments longer - V4 DELETED
                            %KERRI
                            %vector = ceil(vector*min_length);

                            %Check whether the new segment will leave the grid
                            newpoint = CapillaryList{ii}.SegmentList{tipnum}.Node2 + vectorlen;
                            sizegrid = size(voxelgrid.Agent);
                            binary = checkBoundaries3(round(newpoint/pixelsize), sizegrid);

                            %Do Prolif if it doesn't leave grid
                            if binary == 0 %Prolif if it doesn't leave the grid

                                %Check whether it (tip is long enough) can prolif %ADDED BY KERRI V3.1
                                [binaryP] = CapillaryList{ii}.SegmentList{tipnum}.canProliferate;
                                %binaryP = 1;


                                if binaryP == 1 %Check it can prolif
                                    %  display('proliferating!') %REMOVED vector V4.2 KERRI
                                    Listpos = tipnum+1;

                                    %Still need to input tip cell here[newSeg] = proliferate(obj, Listposval, prolifval, pixelsize)
                                    [newtip, newstalk, newEbin] = CapillaryList{ii}.SegmentList{tipnum}.proliferate(prolifval, pixelsize); %removed Listpos KAN 3/12/13

                                    %Replace newstalk %KAN 3/12/13
                                    CapillaryList{ii}.SegmentList{tipnum}= newstalk;
                                    uniquepos = find(SegmentMatrix(:,7) == ii & SegmentMatrix(:,8) == tipnum);

                                    %fill in SegmentMatrix
                                    SegmentMatrix(uniquepos,1:3) =  newstalk.Node1;
                                    SegmentMatrix(uniquepos,4:6) =  newstalk.Node2;
                                    SegmentMatrix(uniquepos,7) =  newstalk.CapillaryNum;
                                    SegmentMatrix(uniquepos,8) = newstalk.Listpos;
                                    SegmentMatrix(uniquepos,9) =  newstalk.Mature;

                                    CapillaryList{ii}.SegmentList{pstalk}.CellClock = 0;
                                    %new.updateDirection(); KAN 3/12/13

                                    if newtip.Node1 == newtip.Node2
                                        error('Nodes are the same')
                                    end

                                    %Add newcell to SegmentList
                                    CapillaryList{ii}.SegmentList{Listpos} = newtip;
                                    %fill in SegmentMatrix
                                    epos = find(SegmentMatrix(:,1) == 0, 1 );
                                    if isempty(epos)
                                        epos = length(SegmentMatrix(:,1)) +1;
                                    end
                                    SegmentMatrix(epos,1:3) =  newtip.Node1;
                                    SegmentMatrix(epos,4:6) =  newtip.Node2;
                                    SegmentMatrix(epos,7) =  newtip.CapillaryNum;
                                    SegmentMatrix(epos,8) = newtip.Listpos;
                                    SegmentMatrix(epos,9) =  newtip.Mature;

                                    %
                                    %% Add new cell location to Agent List %ADDED BY KERRI V1.2
                                    %Find new cells nodes and radius
                                    %                                              Node1pt = newcell.Node1;
                                    %                                              Node2pt = newcell.Node2;
                                    %                                              segradius = newcell.Radius;

                                    %% Test Anastomosis %Added V8
                                    %Check Anast fills in the grid
                                    if newtip.Listpos == 1
                                        %This is the first segment
                                        [doanast, gridindex] = CheckAnastBranch(voxelgrid, newtip, CapillaryList, pixelsize);
                                    else
                                        %Previous segment
                                        oldseg =  CapillaryList{ii}.SegmentList{tipnum};
                                        [doanast, gridindex] = CheckAnast(voxelgrid, newtip, oldseg, CapillaryList, pixelsize);
                                    end

                                    if doanast == 1
                                        %Do anastomosis
                                        [voxelgrid, SegmentMatrix, newseg, CapillaryList, CapillaryList{ii}.ActiveList, MatureNodes] = doAnastSM(newtip,   voxelgrid, SegmentMatrix, gridindex, CapillaryList, CapillaryList{ii}.ActiveList, pixelsize, MatureNodes);
                                        %Update new cell
                                        newseg.updateDirection();
                                        CapillaryList{ii}.SegmentList{Listpos} = newseg;
                                        uniquepos = find(SegmentMatrix(:,7) == ii & SegmentMatrix(:,8) == Listpos);

                                        %fill in SegmentMatrix
                                        SegmentMatrix(uniquepos,1:3) =  newseg.Node1;
                                        SegmentMatrix(uniquepos,4:6) =  newseg.Node2;

                                        if newseg.Node1 == newseg.Node2
                                            error('Nodes are the same')
                                        end
                                    else
                                        %Fill in the grid
                                        [voxelgrid,CapillaryList{ii}.SegmentList{Listpos}] = fillAgentGridTortNodes(voxelgrid, CapillaryList{ii}.SegmentList{Listpos}, pixelsize);
                                    end
                                end % if binaryP == 1
                            else
                            end
                        else %if a node is on the edge of the grid
                        end
                    else
                        %Add one to its cell cycle
                        %Add to the cell cycle
                        newcycle = CapillaryList{ii}.SegmentList{pstalk}.CellClock + 1;
                        CapillaryList{ii}.SegmentList{pstalk}.CellClock = newcycle;
                    end
                end %End if ~isempty(theactive)

                %% Branching
                segmentset = tipnum-2; %only phalanx cells
                if segmentset >0
                    for k = 1:segmentset

                        %noempty = find(~cellfun('isempty',CapillaryList{ii}.SegmentList));
                        %test if it can branch: check VEGF and prob
                        [binarybranch, binarycount] = CapillaryList{ii}.SegmentList{k}.canBranch(XYZ, HypCell, pixelsize);
                        %                 display('len:')
                        %                 norm(CapillaryList{i}.SegmentList{ind}.Node2-CapillaryList{i}.SegmentList{ind}.Node1)
                        %Keep count of attempts at branching
                        branchpossible = branchpossible + binarycount;

                        %binarybranch
                        inputseg = CapillaryList{ii}.SegmentList{k};
                        %if yes, branch
                        if binarybranch == 1
                            branchcount = branchcount + 1;
                            Hypbin = HypCell == 1;
                            XYZh = XYZ(Hypbin,:);
                            [CapillaryList, SegmentMatrix,voxelgrid] = branch2(CapillaryList,SegmentMatrix,inputseg,Agentmat,voxelgrid, XYZh, prolifval, pixelsize);
                            CapillaryList{ii}.SegmentList{k}.Branched = 1;
                            %Salt and Pepper Pattern - deactive 2
                            %neighbors for branching V5
                            if k > 1
                                CapillaryList{ii}.SegmentList{k-1}.Branched = 1;
                            end
                            %k will always be less than end
                            %if k < length(noempty) %NEED TO CHANGE IF APOP
                            CapillaryList{ii}.SegmentList{k+1}.Branched = 1;
                        end
                    end %k=length(segmentset)
                end
            end %end for length(caps) V5
        end % for hrsix = 1:6
    end

    %% Tumor cell module
    %randperm
    numberofcells = size(XYZ,1);
    randcell = randperm(numberofcells);
    jr = 0;
    prolif = 0;
    hypcount = 1;
    %Go through all of the current cells
    for jj = 1:numberofcells %proliferate all cells so far but not the new ones


        %Don't update the numberofcells
        %randomly choose the cells
        jr = randcell(jj);

        %Test if cell is hypoxic
        if XYZ(jr,2) <= 10 %200 microns
            HypCell(jr) = 0; %Not hypoxic
            HypCounter(jr) = 0; %Now counter is back to 0
            Agentmat(XYZ(jr,1), XYZ(jr,2),XYZ(jr,3)) = 1;
        else
            %It might be hypoxic
            %Need pos in microns
            truepos = 20*XYZ(jr,:);
            %Find smallest distance to mature vessel
            % thedisth = pdist2(truepos, MatureNodes);
            thedisth = simpleDist(truepos, MatureNodes);
            % assert(isequal(thedisth(:),thedisthtemp(:)))
            dmin = min(thedisth);
            if dmin >= 200
                %not near a mature vessel - hypoxic
                HypCell(jr) = 1; %hypoxic

                Hypbin=HypCell==1;
                XYZh=XYZ(Hypbin,:);

                HypCounter(jr) = HypCounter(jr) + 1;

                % hyploc(hypcount,:) = XYZ(jr,:);
                hypcount = hypcount + 1;
                %disp('hypoxic cell')
                % Agentmat(XYZ(jr,1), XYZ(jr,2),XYZ(jr,3)) = 5;
            else
                HypCell(jr) = 0; %Not hypoxic
                HypCounter(jr) = 0; %Now counter is back to 0
                % Agentmat(XYZ(jr,1), XYZ(jr,2),XYZ(jr,3)) = 1;
            end
        end
        %% Check if there is space
        emptypos = FindEmptySpace(Agentmat,XYZ(jr,:), indmat);
        migr = rand(1);

        if ~isempty(emptypos)

            %% Migration
            %find closest cell
            %Test distance from stromal cells
            pos = XYZ(jr,:);
            if size(XYZm,1) > 1
                thedistm = pdist2(pos, XYZm);
            else
                thedistm = 200;
            end
            if size(XYZf,1) > 1
                thedistf = pdist2(pos, XYZf);
            else
                thedistf = 200;
            end
            %change these cells migration rate
            themig = MigProb(jr);
            %Find close cancer cells
            if min(thedistm) <= 10
                %change migr rate to 2.5X
                if themig < 1
                    MigProb(jr) = migmultm*0.625;
                end
                %change prolif rates
                if CellState(jr) == 1
                    if HypCell(jr) == 0
                        %stem - increase by 1.25
                        DivRate(jr) = sdiv * 1.25;
                    else
                        %hypoxia decrease by 1/2
                        %stem - increase by 1.25
                        DivRate(jr) = sdiv * 1.25 * 0.5;
                    end
                elseif CellState(jr) == 2
                    if HypCell(jr) == 0
                        %progen - increase by 3.5
                        DivRate(jr) = promultm*pdiv*3.5;
                    else
                        %hypoxia decrease by 1/2
                        DivRate(jr) = promultm*0.5*3.5*pdiv;
                    end
                end

            elseif  min(thedistm) > 10 && min(thedistf) > 10
                %migrate should be normal&&
                if themig == 0.625
                    MigProb(jr) = 0.25;
                end
                %change prolif rates
                if CellState(jr) == 1
                    if HypCell(jr) == 0
                        %stem - normal
                        DivRate(jr) = sdiv;
                    else
                        %hypoxia 1/2
                        %stem - normal
                        DivRate(jr) = sdiv * 0.5;
                    end
                elseif CellState(jr) == 2
                    if HypCell(jr) == 0
                        %progen - normal
                        DivRate(jr) = pdiv;
                    else
                        %hypoxia 1/2
                        %progen - normal
                        DivRate(jr) = 0.5*pdiv;
                    end
                end
            elseif min(thedistm) > 10 && min(thedistf) <= 10
                %change migr rate to 2.5X
                if themig < 1
                    MigProb(jr) = migmultm*0.425;
                end
                %change prolif rates
                if CellState(jr) == 1
                    if HypCell(jr) == 0
                        %stem - increase by 1.25
                        DivRate(jr) = sdiv * 1.25;
                    else
                        %hypoxia 1/2
                        %stem - increase by 1.25
                        DivRate(jr) = sdiv * 1.25 * 0.5;
                    end
                elseif CellState(jr) == 2
                    if HypCell(jr) == 0
                        %progen - increase by 3.5
                        DivRate(jr) = promultm*pdiv*3.5;
                    else
                        %hypoxia half
                        %progen - increase by 3.5
                        DivRate(jr) = promultm*0.5*pdiv*3.5;
                    end
                end
            end

            %Hypoxia
            %Test if the cell is hypoxic - if so change Migspeed to 3X and
            %divrate to half
            %cellpos = theCells{j}.Position;
            cellCCR5 = CCR5level(jr);
            if HypCell(jr) == 1
                if cellCCR5 ==1
                    MigSpeed(jr) = 30;
                else
                    MigSpeed(jr) = 3;
                end
            else
                if cellCCR5 ==1
                    MigSpeed(jr) = 10;
                else
                    MigSpeed(jr) = 1;
                end
            end

            %Change these cells proliferation rate

            %Check migration cycles are the same - now migration happens
            %every iteration
            if DoMigration == 1 && migr <= MigProb(jr)%theCells{j}.MigCycle == theCells{j}.MigClock &&

                %Check how many cell spaces to move
                migmoves = MigSpeed(jr) - 1;

                %Record old position
                oldpos = XYZ(jr,:);

                %Update new position
                XYZ(jr,:) = emptypos;

                %for loop for number of moves
                if migmoves > 0
                    for mi = 1:migmoves
                        %find next empty spot
                        emptypos2 = FindEmptySpace(Agentmat,XYZ(jr,:), indmat);
                        if ~isempty(emptypos2)
                            %Update new position
                            XYZ(jr,:) = emptypos2;
                        end
                    end %mi = 1:migmoves
                end %if migmoves >0

                %Remove previous position from Agentmat
                Agentmat(oldpos(1), oldpos(2), oldpos(3)) = 0;

                %Put new position in Agentmat
                newpos = XYZ(jr,:);
                Agentmat(newpos(1), newpos(2), newpos(3)) = 1;

            else
                %It isn't time to migrate - increment migration cycle
                %MigCycle(j) = MigCycle(j) +1;
            end
            % end %THIS END IS FOR TEST ONLY REMOVE !!!!!!!!!!!!!!!!!!!!!
            %% Proliferate Cell
            % [startprolif, newpos] = theCells{j}.DoProlif(Agentmat, indmat); %V4.1
            %% Symmetric Stem Cell Division %V4.1
            if CellState(jr) == 2 %Progen
                if DivCycle(jr) < DivLimit(jr) %a progen cell can only prolif 6 times (Potten, 2003)
                    %Check that there is room to prolif
                    newpos = FindEmptySpace(Agentmat,XYZ(jr,:), indmat);

                    if ~isempty(newpos)
                        % division rate of progens per day 50% (Ashkenazi, 2008)
                        %Active
                        QState(jr) = 1;
                        if rand < DivRate(jr)
                            startprolif = 1; %Do Prolif
                            % disp('Do progen prolif')

                        else
                            % disp('Dont do progen prolif')
                            startprolif = 0; %Do Not Prolif

                        end
                    else
                        startprolif = 0; %Do Not Prolif
                        %Quiesce
                        QState(jr) = 2;
                    end %~isempty(newpos)

                else
                    %disp('Niche reached max prolif - Apoptosis');
                    QState(jr) = 0;%obj.ChangeToDead;
                    startprolif = -1; %Cell is dead
                    newpos = {}; %There is no newpos
                    %Add one to MigCycle --- WHY???
                    %                      cycnum = MigCycle(j)+1;
                    %                      obj.MigCycle = mod(cycnum,4);
                end %obj.Cycle < 60

            else %Stem Cell
                %                 %Check if cell can proliferate
                %                 %Check that there is room to prolif
                newpos = FindEmptySpace(Agentmat, XYZ(jr,:), indmat);

                if  ~isempty(newpos)
                    %             %Stem Cells can proliferate forever
                    rnum = rand(2,1); %rand number
                    %rnum2 = rand(1);
                    %Stem Cells asym divide 90% of the time, 10% sym
                    %Cancer Stem Cell rate of Division per day .02
                    %.02*.9 = .018 prob of asym division per day
                    %NEEDS TO CHANGE - TOO SMALL
                    if rnum(1) < DivRate(jr) %Ref: (Ashkenazi, 2008) & (Enderling, 2007)
                        if rnum(2) < symchange

                            %Do sym division
                            startprolif = 2;
                            QState(jr) = 1;
                            %The progen's Cycle should start at 0
                            %Therefore you need to reset stem cell cycle to 0
                            DivCycle(jr) = 0;

                        else
                            %Do asym division
                            startprolif = 1; %Prolif
                            %Become active
                            QState(jr) = 1; %V4.1
                            %Therefore you need to reset stem cell cycle to 0
                            DivCycle(jr) = 0;
                        end
                        %binaryans = 1; %Do Prolif

                        %                  disp('Cycle Reset')
                        %                  obj.State = 1; %Become active %V4.1

                    else
                        startprolif = 0; %Don't Prolif
                        newpos = {}; %No newpos
                    end
                else
                    % disp('stem becomes quies')
                    startprolif = 0; %Don't Prolif
                    newpos = {}; %No newpos
                    %Quiesce
                    QState(jr) = 2; %V4.1
                end %~isempty(newpos)
            end % if theCells{j}.StemType == 1

            %% Proliferation V4.1

            if startprolif == 1 %CHANGED V4.1 %Produce Progen Cell
                %proliferate
                % display('progen prolif')
                %[newcell, Agentmat] = Prolif(obj,newpos, Agentmat,migmove)
                % if CellState(j) == 2 %% if a progen
                %The cell can create a new cell

                %Define firstcell's properties
                %Find the new position %CHANGED V2.1
                %newpos = obj.FindEmptySpace(grid);
                numberofcells = numberofcells +1;
                prolif = prolif+1;
                empspace = numberofcells;
                XYZ(empspace,:) = newpos;
                CellState(numberofcells) = 2;
                %grid(newpos(1),newpos(2),newpos(3)) = 1; %ADDED V2.1
                DivCycle(empspace) = DivCycle(jr) + 1; %The new cell's cycle is the Mothercell's cycle plus one

                %5% chance of being CCR5+ =
                if rand(1) <= 0.05
                    CCR5level(empspace) = 1;
                    MigSpeed(empspace) = migmove;
                    numberofCCR5 = numberofCCR5 +1;
                    MigProb(empspace) = 1;
                else
                    CCR5level(empspace) = 0;
                    MigSpeed(empspace) = 1;
                    MigProb(empspace) = 0.25;
                end

                %Update the cell cycle +1
                if (DivCycle(jr) < 0) %test Cycle is valid
                    error('Cell entry has bad Cycle (< 0)');
                end
                DivCycle(jr) = DivCycle(jr) + 1;
                QState(empspace) = 1;
                DivLimit(empspace) = dlim;
                DivRate(empspace) = pdiv;%0.5;
                HypCounter(empspace) = 0;
                HypCell(empspace) = 0;
                AntigenHet(empspace) = AntigenHet(jr);
                %CellSymRate = 0;
                %Place cell on Agentmat       %V3.1
                Agentmat(newpos(1), newpos(2), newpos(3)) = 1;

            elseif startprolif == 2  %Produce Stem Cell
                % display('stem prolif')
                %proliferate stem sym
                % [theCells{length(theCells) +1}, Agentmat] = theCells{j}.SymProlif(newpos, Agentmat,migmove, sdiv); %V4.1
                %The cell can create a new cell
                %Define initial metastatic stem cell
                prolif = prolif+1;
                numberofcells = numberofcells+1;

                DivRate(numberofcells) = sdiv;
                %Define firstcell's properties
                XYZ(numberofcells,:) = newpos;
                CellState(numberofcells) = 1;
                DivCycle(numberofcells) = 0; %Stem cell's cycle is always 0
                DivLimit(numberofcells) = 108000;%No limit
                QState(numberofcells) = 1;
                HypCounter(numberofcells) = 0;
                HypCell(numberofcells) = 0;
                AntigenHet(numberofcells) = AntigenHet(jr);
                %5% chance of being CCR5+
                if rand(1) <= 0.05
                    CCR5level(numberofcells) = 1;
                    MigSpeed(numberofcells) = migmove;
                    numberofCCR5 = numberofCCR5 +1;
                    MigProb(numberofcells) = 1;
                else
                    CCR5level(numberofcells) = 0;
                    MigSpeed(numberofcells) = 1;
                    MigProb(numberofcells) = 0.25;
                end

                %Place cell on Agentmat %V4.1
                Agentmat(newpos(1), newpos(2), newpos(3)) = 1;

                %UPDATE SYMDIV
                %Update the symRate
                %CellSymRate(numberofcells) = symchange;
                %Add the number of stems
                numberofstems = numberofstems + 1;

            elseif startprolif == 0
                %display('Cell does not prolif')
            else %should be -1;
                %Keep track of the senescent
                possiblydead = [possiblydead, jr];
            end

        else
            %The Cell is quiesecent - Nothing is updated

        end %~isempty(emptypos)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CAR-T Module
    CARTkills = [];
    if size(XYZcart,1) > 0
        [numRows,numCols] = size(XYZcart);
        numberofcellscart = numRows;
        %randperm
        randcellcart = randperm(numberofcellscart);
        jcart = 0;
        killcount = 0;
        CARTApop = zeros(size(XYZcart));

        %probability of a CAR-T cell becoming exhausted form killing
        %CTexhaustionLevel = .000001; %taken out for speed bc inactivation rate is neglegable

        %Go through all of the current cells
        for jjj = 1:numberofcellscart %proliferate all cells so far but not the new ones
            %Don't update the numberofcells
            %randomly choose the cells
            jcart = randcellcart(jjj);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%             Migration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Version 2
            %Check if there is empty space
            %% Check if there is space


            %Note: jm is the random fibroblast position
            cartmigchance = 1; %was.25;
            cartkillchance = 1;%was.25;
            cartprolifchance = .2; %was .75 %https://www.medrxiv.org/content/10.1101/2020.05.09.20096586v2.full#F4
            cartnumkills = 5;

            %% CAR T Migration
            if rand < cartmigchance
                Cmigrate = 16;
                [emptypos] = FindEmptySpace(Agentmat,XYZcart(jcart,:), indmat);
                if ~isempty(emptypos)
                    %Test distance from tumor cells
                    %cpos =  XYZcart(jcart,:);
                    cpos =  double(XYZcart(jcart,:));
                    oldpos = cpos;

                    thedists = pdist2(cpos, XYZ);
                    [mindist, minidx] = min(thedists);
                    isclose = mindist<10;
                    %minpos = int64(XYZ(minidx,:));
                    minpos = XYZ(minidx,:);
                    v=(minpos-cpos)/sqrt((minpos-cpos)*(minpos-cpos)'); %find the unit vector

                    while Cmigrate > 0
                        if isclose %if the tcell is close to a cancer cell move toward it
                            if ~isempty(FindTarget(Agentmat,cpos, indmat,1))%if next to cancer dont move
                                break;
                            end

                            testpos = cpos + round(v);
                            isfree = Agentmat(testpos(1), testpos(2), testpos(3))==0; %test if it can move where it wants
                            if isfree
                                cpos = testpos;
                                v=(minpos-cpos)/sqrt((minpos-cpos)*(minpos-cpos)');%maybe remove if slow
                            else
                                %move randomly
                                emptypos2 = FindEmptySpace(Agentmat,cpos, indmat);
                                if ~isempty(emptypos2)
                                    cpos = emptypos2;
                                end
                                break;
                            end

                        else %not close to cell; move randomly
                            emptypos2 = FindEmptySpace(Agentmat,cpos, indmat);
                            if ~isempty(emptypos2)
                                cpos = emptypos2;

                                %check if it has moved within distance of a cancer cell
                                thedists = pdist2(cpos, XYZ);
                                [mindist, minidx] = min(thedists);
                                isclose = mindist<10;
                                if isclose %set up vector toward what it is close to
                                    minpos = XYZ(minidx,:);
                                    v=(minpos-cpos)/sqrt((minpos-cpos)*(minpos-cpos)'); %find the unit vector
                                end
                            else
                                break; %no open space. Its going nowhere
                            end
                        end
                        Cmigrate = Cmigrate-1;
                    end %end migration steps

                    %update simulation variables
                    Agentmat(oldpos(1), oldpos(2), oldpos(3)) = 0;
                    Agentmat(cpos(1), cpos(2), cpos(3)) = 4;
                    XYZcart(jcart,:) = cpos;
                end
            end

            %% Killing
            if rand < cartkillchance
                [target1] = FindTarget(Agentmat,XYZcart(jcart,:), indmat,cartnumkills);
                if ~isempty(target1)
                    [~,deadidx] = ismember(target1,XYZ,'rows');
                    deadidx = deadidx(rand(length(deadidx),1) < AntigenHet(deadidx)); %filter out non-antigen presenting cells
                    CARTkills = [CARTkills,deadidx'];
                end

                %car-t has chance of becoming exhausted and is marked
                %to be killed
                %taken out for speed bc inactivation rate is neglegable
                %if rand < CTexhaustionLevel
                %      CARTTimer(jcart) = 0;
                %end
            end

            %% Proliferation
            if rand < cartprolifchance && CARTTimer(jcart)> cartLifespan-cartProlifDuration
                [len_XYZcart, ~] = size(XYZcart); %accounts for 1d array
                empspacecart = len_XYZcart+1;
                newposcart = FindEmptySpace(Agentmat, XYZcart(jcart,:), indmat);
                if isempty(newposcart) == 0
                    cartnum2 = cartnum2 + 1;
                    XYZcart(empspacecart,:) = newposcart;
                    indcart = sub2ind(size(Agentmat),newposcart(:,1), newposcart(:,2), newposcart(:,3));
                    Agentmat(indcart) = 4;
                    CARTTimer(empspacecart) = cartLifespan;
                end
            end
            CARTTimer(jcart) = CARTTimer(jcart)-1;
            %CART random death
            if numberofcellscart < CartBreakpt
                if rand <= cartRandDeath
                    CARTApop(jcart)=1;
                end
            else
                MM = (.225-cartRandDeath)*(numberofcellscart-CartBreakpt)/(7000 + numberofcellscart-CartBreakpt);
                modifiedCartRandDeath = MM+cartRandDeath;
                if rand<modifiedCartRandDeath
                    CARTApop(jcart)=1;
                end
            end
            if CARTTimer(jcart) <= 0
                CARTApop(jcart) = 1;
            end
        end

        for kap = length(CARTApop):-1:1
            if CARTApop(kap) == 1
                cartapos = XYZcart(kap,:);
                cartaind = sub2ind(size(Agentmat),cartapos(:,1), cartapos(:,2), cartapos(:,3));
                Agentmat(cartaind) = 0;
                XYZcart(kap,:) = [];
                cartnum2 = cartnum2-1;
                CARTTimer(kap) = [];
            end
        end

        CARTkills = unique(CARTkills);
        killcount = killcount + length(CARTkills);
        possiblydead2 = [possiblydead2 CARTkills];
    else
        killcount = 0;
    end

    numberofstems = length(find(CellState== 1));
    numberofcells = length(CellState);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Allow macrophage infiltration
    %Find Nodes close to cancer cells
    %find the distance of Nodes to cancer cells
    MNdist = pdist2(MatureNodes,20*XYZ);
    MNmins = min(MNdist,[],2);
    %find close Nodes
    MNclose = MNmins <200;
    closeNodes = MatureNodes(MNclose,:);
    %randperm
    mnrand = randperm(size(closeNodes,1));
    %How many are infiltrated?
    numinf = round(MacroInfProb*size(closeNodes,1));
    if numinf > 1
        %declare death
        MatNod2 = closeNodes(mnrand(1:numinf),:);
        for mm = 1:size(MatNod2,1)
            %if rand(1) < MacroInfProb
            %Find mature nodes Agentmat pos
            mnpos = ceil(MatNod2(mm,:)/20); %NX3
            %Find empty agentmat spaces near mature nodes
            %need to change to inds and back again
            mninds = sub2ind(size(Agentmat),mnpos(:,1),mnpos(:,2),mnpos(:,3));
            theempt = Agentmat(mninds)==0;
            eminds = mninds(theempt);
            if ~isempty(eminds)
                %choose a random empty space
                randind = eminds(randperm(length(eminds),1));
                %make new macro
                [xmn, ymn, zmn] =ind2sub(size(Agentmat),randind);
                addnew = length(Statem)+1;
                XYZm(addnew,:) = [xmn, ymn, zmn];
                Statem(addnew) = 2;

                %add macro to list and grid
                Agentmat(randind) = 2;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TIL Infiltration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Apoptosis Cell
    %Apoptosis must be done separately - TO DO
    %Kill off hypoxic cells
    toohyp = find(HypCounter > 40);

    %Maybe the cell is dead
    if isempty(possiblydead) == 0 %check there are possiblydeads
        %       for j = 1: length(possiblydead) %go through indexes
        %          %Apoptosis
        %          %Check to see if cell is dead
        %          if theCells{possiblydead(j)}.State ~= 0
        %              %if the state is not dead remove it from the list
        %              possiblydead(j) = 0;
        %          end
        %       end
        %They should all be senescent at this point
        possiblydead(possiblydead == 0) = [];
        %randperm
        pdrand = randperm(size(possiblydead,1));
        %How many die?
        numdead = round(sendeath*size(possiblydead,1));
        deathovertime(time) = numdead;
        if numdead > 1
            %declare death
            possiblydead2 = [possiblydead2, possiblydead(pdrand(1:numdead))];
        elseif numdead == 1
            posp = 1;
            %possiblydead2 =[]; %breaks CT killing
            for testd = 1:length(possiblydead)
                rd = rand(1);
                if rd <sendeath
                    possiblydead2(posp) = possiblydead(testd);
                    posp = posp+1;
                end
            end
        else
            %possiblydead2 = []; %breaks CT killing
        end
    end
    %CCR5dead = 0;
    %for pd = 1:length (possiblydead)
    %Append hypoxic death to possibly dead;
    %if isempty(possiblydead) == 1
    if isempty(toohyp) == 0
        possiblydead2 = [possiblydead2, toohyp'];
    else
        % possiblydead2 = [];
    end
    %end
    %get rid of replicates
    possiblydead2 = unique(possiblydead2);
    deathovertime(time) = length(possiblydead2);

    if isempty(possiblydead2) == 0
        %find location
        deadpos = XYZ(possiblydead2,:);
        deadinds = sub2ind(size(Agentmat),deadpos(:,1), deadpos(:,2), deadpos(:,3));
        %remove from grid
        Agentmat(deadinds) = 0;
        %test CCR5 level
        CCR5dead = sum(CCR5level(possiblydead2));

        %end
        %theCells(possiblydead2) = [];
        %theCells = zeros(100,1);
        CellState(possiblydead2) = []; %1 is stem and 2 is progen
        QState(possiblydead2) = []; %Whether the cell is alive/active (1), quiescent (2), or senescent (0)
        %CellSymRate(possiblydead2) = [];
        CCR5level(possiblydead2) = [];
        MigSpeed(possiblydead2) = [];
        MigProb(possiblydead2) = [];
        HypCell(possiblydead2) = [];
        HypCounter(possiblydead2) = [];
        XYZ(possiblydead2,:) = [];
        %MigCycle = zeros(500000,1);
        DivLimit(possiblydead2) = [];
        DivRate(possiblydead2) = [];
        DivCycle(possiblydead2) = [];
        AntigenHet(possiblydead2) = [];
        %subtract dead CCR5 cells
        numberofCCR5 = numberofCCR5 - CCR5dead;
        numberofstems = length(find(CellState== 1));
        %       disp('Apoptosis has occured')
        %       disp(length(theCells))
        numberofcells = length(CellState);
    end
    Hypbin=HypCell==1;
    XYZh=XYZ(Hypbin,:);
    hyptest = sum(Hypbin);

    %Save the numberofcells to an array
    cellsovertime(time) = numberofcells;
    stemsovertime(time) = numberofstems;
    CCR5overtime(time) = numberofCCR5;
    hypoxiaovertime(time)= hyptest;
    macroovertime(time) = numberofcellsm;
    prolifovertime(time) = prolif;
    killsovertime(time)=killcount;
    vascovertime(time) =  sum(sum(sum(voxelgrid.Agent)));

    %Save the numberofcell to txt file
    %SaveDataarray each iteration

    if time == CARTDelay
        sizegrid = numel(Agentmat);
        rorder3 = randperm(sizegrid);
        thecart = rorder3(1:cartnum);
        Agentmat(thecart) = 4;
        [Icart,Jcart,Kcart] = ind2sub(size(Agentmat),thecart);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        theinds = [Icart;Jcart;Kcart];

        XYZcartold = theinds';
        XYZcart = GetCartSpawnLoc(MatureNodes,cartnum,XYZ,indmat);


        Statecart = 4*ones(cartnum,1);%unknown if needed
        CARTTimer = cartLifespan*ones(cartnum*2,1);
        CARTApop = zeros(cartnum,1);
    end

    cartovertime(time) = size(XYZcart,1);

    %         c = 0;
    %         for x = 1:50
    %             for y = 1: 50
    %                 for z = 1:50
    %                     if Agentmat(x,y,z)==4
    %                         c = c+1;
    %                     end
    %                 end
    %             end
    %         end
    %         disp("Agentmat T Cell Count")
    %         disp(num2str(c))
    %         disp('Counter CART Count')
    %         disp(num2str(cartovertime(time)))
    %         disp('XYZcart length')
    %         disp(num2str(length(XYZcart)))

    if numberofcells > 500000
        break
    end

    if numberofcells <= 0
        break
    end

    %Added PlotVascV1.m to code MV/JBU
    if mod(time, print_every) == 0
        tempVG1 = voxelgrid.Agent(2,:,:);
        voxelgrid.Agent(2,:,:) = voxelgrid.Agent(1,:,:)|voxelgrid.Agent(2,:,:);

        tempVG2 = voxelgrid.Agent(:,2,:);
        voxelgrid.Agent(:,2,:) = voxelgrid.Agent(:,1,:)|voxelgrid.Agent(:,1,:);

        tempVG3 = voxelgrid.Agent(:,:,2);
        voxelgrid.Agent(:,:,2) = voxelgrid.Agent(:,:,1)|voxelgrid.Agent(:,:,2);

        tempVG4 = voxelgrid.Agent(end-1,:,:);
        voxelgrid.Agent(end-1,:,:) = voxelgrid.Agent(end,:,:)|voxelgrid.Agent(end-1,:,:);

        tempVG5 = voxelgrid.Agent(:,end-1,:);
        voxelgrid.Agent(:,end-1,:) = voxelgrid.Agent(:,end,:)|voxelgrid.Agent(:,end-1,:);

        tempVG6 = voxelgrid.Agent(:,:,end-1);
        voxelgrid.Agent(:,:,end-1) = voxelgrid.Agent(:,:,end)|voxelgrid.Agent(:,:,end-1);

        figure
        hold on
        xlim([0 gridsize(1)])%Set lims so matlab does not graph over Mv
        ylim([0 gridsize(2)])
        zlim([0 gridsize(3)])
        p=patch(isosurface(voxelgrid.Agent==1,0));
        set(p,'facecolor','red' ,'edgecolor', 'none');
        daspect([1 1 1])
        %rotate3d on;
        view([60 30]);
        camlight
        lighting gouraud
        hold off
        cd Data
        cd(exType)
        cd(experimentfolder)
        cd(subfolder)
        cd Vasculature_Images
        SaveFig = strcat('Figure_', num2str(time), '.fig');
        %print ('-dtiff', '-r300', SaveFig)
        savefig(SaveFig)
        close all
        voxelgrid.Agent(2,:,:) = tempVG1; % Matlab will not graph the first so we logical or with the second row so it plots
        voxelgrid.Agent(:,2,:) = tempVG2; %BSRI 19
        voxelgrid.Agent(:,:,2) = tempVG3;
        voxelgrid.Agent(end-1,:,:) = tempVG4;
        voxelgrid.Agent(:,end-1,:) = tempVG5;
        voxelgrid.Agent(:,:,end-1) = tempVG6;
        cd ..

        %Save Hypoxic cell locations in txt file, note hypcount is the number
        %of hypoxic cells + 1 due to indexing- MV/HF
        cd OverTime
        %  hyploc(hypcount:end, :) = [];
        name6 =strcat('XYZloc_t_',num2str(time),'.txt');
        % dlmwrite(name6, XYZ, 'delimiter', '\t')

        name3b =strcat('XYZh_t_',num2str(time),'.txt');
        % dlmwrite(name3b, XYZh, 'delimiter', '\t')

        namedataH1 =strcat('HypCell_t_',num2str(time),'.txt');
        % dlmwrite(namedataH1, HypCell, 'delimiter', '\t')

        name_cart = strcat('XYZcart_t_',num2str(time),'.txt');
        % if ~isempty(XYZcart>0)
        %     % dlmwrite(name_cart, XYZcart, 'delimiter','\t')
        % else
        %     % dlmwrite(name_cart, [], 'delimiter','\t')
        % end

        namehet = strcat("AntigenExp_t_",num2str(time),'.txt');
        % dlmwrite(namehet,AntigenHet, 'delimiter','\t')

        name4 =strcat('State','_CART',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
        % dlmwrite(name4, CellState, 'delimiter', '\t')

        %Save voxel grid for plotting vasculature later in mat file-MV
        name8 =strcat('voxelgrid_t_',num2str(time));
        save( name8, 'voxelgrid');

        cd ../Run_images
        SaveAddress2 = strcat('CapListp_t_', num2str(time));
        CapMatrix = CapillaryList;
        %     saveCapMatrix(SaveAddress2, CapMatrix);
        save(SaveAddress2, 'CapMatrix');

        cd ../../../../..
    end
end
%origional save data


sdiv2 = sdiv*100;
promult2 = round(promultm*35);
migmult2 = round(migmultm*25);



seedrate2 = seedrate*100;

cd('Data')
cd(exType)
cd(experimentfolder)
cd(subfolder)
if ~exist("fig_data","dir")
    mkdir('fig_data')
end
% if ~exist("end_of_run","dir")
%     mkdir('end_of_run')
% end

cd('fig_data')

dlmwrite(namedataA, cellsovertime, 'delimiter', '\t')

namedataD =strcat('NumberofDeath','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataD, deathovertime, 'delimiter', '\t')

namedataS =strcat('NumberofStems','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataS, stemsovertime, 'delimiter', '\t')

namedataH =strcat('NumberofHypoxic','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataH, hypoxiaovertime, 'delimiter', '\t')

namedataH1 =strcat('HypCell','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataH1, HypCell, 'delimiter', '\t')

namedataM =strcat('NumberofMacro','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataM, macroovertime, 'delimiter', '\t')

namedataC =strcat('NumberofCCR5','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataC, CCR5overtime, 'delimiter', '\t')

namedataK =strcat('NumberofCART','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataK, cartovertime, 'delimiter', '\t')

namedataP =strcat('NumberofProlif','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataP, prolifovertime, 'delimiter', '\t')

namedataKills =strcat('NumberofKills','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataKills, killsovertime, 'delimiter', '\t')

namedataVasc =strcat('NumberofVasc','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(namedataVasc, vascovertime, 'delimiter', '\t')

namehet = strcat("AntigenExp_t_",num2str(time),'.txt');
dlmwrite(namehet,AntigenHet, 'delimiter','\t')


% cd ..
% cd('end_of_run')


name3 =strcat('XYZ','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(name3, XYZ, 'delimiter', '\t')

name3b =strcat('XYZh','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(name3b, XYZh, 'delimiter', '\t')

name4 =strcat('Stateh','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
dlmwrite(name4, CellState, 'delimiter', '\t')

%Also save capillaries
SaveAddress2 = strcat('CapListp', '_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time));
CapMatrix = CapillaryList;
% saveCapMatrix(SaveAddress2, CapMatrix);
save(SaveAddress2, 'CapMatrix');

name2 =strcat('CellMatrix','_STEM','_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time));
% saveAgentmat(name2, Agentmat);
save( name2, 'Agentmat');

% namedata =strcat('rundata',num2str(runnum),'_pdiv',AA_str,'_migmulti',BB_str,'_t',num2str(time),'.txt');
% fid = fopen(namedata,'w');
%     fprintf(fid,'%6s %12.8f\n','SeedRate: ',seedrate);
%     fprintf(fid,'%6s %12.8f\n','SenDeath: ',sendeath);
%     fprintf(fid,'%6s %12.8f\n','MigrationOn: ', DoMigration);
%     fprintf(fid,'%6s %12.8f\n','DivLimit: ',DivLimit(1));
%     fprintf(fid,'%6s %12.8f\n','SymRate: ',0.2);
%     fprintf(fid,'%6s %12.8f\n','Stem Divrate: ',DivRate(1));
%     fprintf(fid,'%6s %12.8f\n','InitPos: ',XYZ(1,:));
%     fprintf(fid,'%6s %12.8f\n','MacroNum: ',macronum);
%     fprintf(fid,'%6s %12.8f\n','FinalIt: ', time);
cd ../../../../..
