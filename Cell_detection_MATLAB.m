%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%1)this script analyses the data acquired using: Cell_detection_Fiji.ijm
%2)it takes the mean of all slices of one cell or the larg est slice of the cell
%3)and calculates the corresponding intensity values for the green/red/blue channels
%4)determines which cells exceed a predefined threshold in the red and blue channels
%5)final data is in FinalData


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% adjust the params below:  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if you want images showing the Green/Red/Blue intensity values choose 1; otherwise 2
IwantImages = 2;

%I want to take the largest plane of the cel=1 ; I want to take the mean=2
IwantZplane = 2;

thresholdRed=1000;
thresholdBlue=10;
thresholdBlueLarge=100;%not included in the analyses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% VARIABELS  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%variables in the FIJI files
numbers=1;
area=2;
meanLum=3;
stdLum=4;
minLum=5;
maxLum=6;
xcoord=7;%7 and 8 are centrum (average of x and y points)
ycoord=8;%11 and 12 are the upper left of the bounding box
circ=15;
medianLum=16;
slice=17;

finalRoisRedVarThres=zeros(300);
finalRoisRedFixThres=zeros(300);

var=3.5;
roinumbers=0;
count=0;

filenames = {'_Blue','_GREEN','_Red'};
fnames = dir('*.txt');
numfids = length(fnames);
vals = cell(1,numfids);

FinalData(1,:) = cellstr({'stack','nr of slices','nr GREEN','nr RED+','nr RED-','nr BLUE+','nr BLUE-','nr RED+/BLUE+',....
    'int GREEN', 'int RED', 'int RED+', 'int RED-','int BLUE', 'int BLUE+', 'int BLUE-',....
    'RED+/BLUE- RED int', 'RED+/BLUE+ RED int', 'RED-/BLUE+ BLUE int','RED+/BLUE+ BLUE int'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% START LOOP  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filesteller = 1;
for ii = 1:3:numfids;
    clear AsortedG
    clear singleRoi
    clear FosPosData
    clear singleRoiZ
    clear FosPosDataZ
    
    filesteller = filesteller+1;
    
    B = importdata((fnames(ii,1).name)); %Blue
    G = importdata((fnames(ii+1,1).name)); %Green
    R = importdata((fnames(ii+2,1).name)); %Red
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%first broad check of x positions lying close to each other%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %AsortedG combines the data from the Green Red and Blue channel
    AsortedG=G.data;
    AsortedG(:,21:24)= R.data(:,[2 3 4 15]);%area mean std circ
    AsortedG(:,25:28)= B.data(:,[2 3 4 15]);%area mean std circ
    
    %clear column 5 and sort file according to the coordinates of the green cells
    AsortedG(1:length(AsortedG),5)=0;%minval -> wordt allemaal nul
    AsortedG=sortrows(AsortedG,xcoord);%xcoord
    AsortedG=sortrows(AsortedG,ycoord);%ycoord
    
    NrOfSlices=max(AsortedG(:,slice));
    
    %fill column 5 with numbers; same number is most likely the same cell
    %this is only a rough check on the basis of xpos, ypos and slice
    cellTeller = 1;
    for t = 1:length(AsortedG);
        xvar = AsortedG(t,xcoord);
        yvar = AsortedG(t,ycoord);
        sliceVar = AsortedG(t,slice);
        
        xvarMin = xvar-var;
        xvarMax = xvar+var;
        yvarMin = yvar-var;
        yvarMax = yvar+var;
        sliceMin = sliceVar-1;
        sliceMax = sliceVar+1;
        
        if AsortedG(t,5)==0;
            %find every entry that lies within the boundaries for xcoord,
            %ycoord and lies one liesabove or below the current
            test=find(AsortedG(:,xcoord)>xvarMin & AsortedG(:,xcoord)<xvarMax & AsortedG(:,ycoord)>yvarMin & AsortedG(:,ycoord)<yvarMax & AsortedG(:,slice)>=sliceMin & AsortedG(:,slice)<=sliceMax);
            
            if length(test)==1; %there is only one roi for this cell thus give a cell number specific for this roi
                cellTeller=cellTeller+1;%to give a specific number increase sliceteller with one
                AsortedG(t,5)=cellTeller;%give this entry that number
                
            else check = find(AsortedG(test,5)>0);%for multiple hits, check whether they already have a number
                if isempty(check);%if not
                    cellTeller=cellTeller+1;%you want to give all these entries a specific number; so increase with one
                    AsortedG(test,5)=cellTeller;%and give this number to the respective entries
                    
                else cellnumber = AsortedG(test(check(1,1),1),5);%in case a number is already present check which nr this is
                    AsortedG(test,5)=cellnumber;%and give the newly identified entries this number
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%calculate mean int and area for the red,green,blue channel  %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %sort again on the basis of cells belonging together (i.e. different
    %slices in one stack) and take the mean to get cell means
    AsortedG=sortrows(AsortedG,2);
    AsortedG=sortrows(AsortedG,5);
    
    if IwantZplane == 1;
        %%here we only include the largest z-plane
        for iii = 1:length(AsortedG);
            check=size(find(AsortedG(:,5)==iii),1);%go through all cell numbers in column 5;find cells with the same code
            if check==0;
            elseif check==1;
                singleRoi(iii,2:size(AsortedG,2))=AsortedG(find(AsortedG(:,5)==iii),2:size(AsortedG,2));
                singleRoi(iii,1)=AsortedG(find(AsortedG(:,5)==iii),1);
                %now take the largest area (or z-plane)
                %area of the NISSL roi is column 2
            else entry=max(find(AsortedG(:,5)==iii));
                singleRoi(iii,2:size(AsortedG,2))=AsortedG(entry,2:size(AsortedG,2));
                %%nanmean(AsortedG(find(AsortedG(:,5)==iii),2:size(AsortedG,2)));% take the mean when multiple cells with the sam code
                singleRoi(iii,1)=min(AsortedG(find(AsortedG(:,5)==iii),1));
            end
        end
    else
        %%here means are created
        for iii = 1:length(AsortedG);
            check=size(find(AsortedG(:,5)==iii),1);%go through all cell numbers in column 5;find cells with the same code
            checkAll(iii,1)=check;%get an idea of how many parts each cell is comprised
            if check==0;
            elseif check==1;
                singleRoi(iii,2:size(AsortedG,2))=AsortedG(find(AsortedG(:,5)==iii),2:size(AsortedG,2));
                singleRoi(iii,1)=AsortedG(find(AsortedG(:,5)==iii),1);
            else singleRoi(iii,2:size(AsortedG,2))=mean(AsortedG(find(AsortedG(:,5)==iii),2:size(AsortedG,2)),"omitmissing");% take the mean when multiple cells with the same code
                singleRoi(iii,1)=min(AsortedG(find(AsortedG(:,5)==iii),1));
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%  fixed threshold and baseline correction             %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sortSingleRoi=singleRoi;
    
    sortSingleRoi(find(sortSingleRoi(:,22)>=thresholdRed),36)=1;%rood
    sortSingleRoi(find(sortSingleRoi(:,26)>=thresholdBlue),37)=1;%blauw
    sortSingleRoi(:,39)=sortSingleRoi(:,36).*sortSingleRoi(:,37);%beide
    
    %sorteren op rood intensiteit en de nummers van de rois die de threshold overschrijden in een final
    %data frame zetten
    sortSingleRoi=sortrows(sortSingleRoi,22);
    finalRoisRedFixThres(1:length(sortSingleRoi(find(sortSingleRoi(:,36)==1),1)),filesteller)=sortSingleRoi(find(sortSingleRoi(:,36)==1),1);
    
    %%%dataset met alleen de FOS+ cellen
    FosPosData = sortSingleRoi(find(sortSingleRoi(:,22)>=thresholdRed),:);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%  create final datasets                   %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currFilename=struct2cell(fnames(ii,1));
    currName = currFilename(1,1);
    
    FinalData(filesteller,1)=currName;
    FinalData(filesteller,2)=num2cell(NrOfSlices);
    FinalData(filesteller,3)= num2cell(max(sortSingleRoi(:,5)));%nr GREEN
    FinalData(filesteller,4)= num2cell(size(find(sortSingleRoi(:,36)==1),1));%nr RED+
    FinalData(filesteller,5)= num2cell(size(find(sortSingleRoi(:,36)==0),1));%nr RED-
    FinalData(filesteller,6)= num2cell(size(find(sortSingleRoi(:,37)==1),1));%nr BLUE+
    FinalData(filesteller,7)= num2cell(size(find(sortSingleRoi(:,37)==0),1));%nr BLUE-
    FinalData(filesteller,8)= num2cell(size(find(sortSingleRoi(:,39)==1),1));%nr RED+/BLUE+
    FinalData(filesteller,9)= num2cell(mean(sortSingleRoi(:,3)));%int GREEN
    FinalData(filesteller,10)= num2cell(mean(sortSingleRoi(:,22)));%int RED
    FinalData(filesteller,11)= num2cell(mean(sortSingleRoi(find(sortSingleRoi(:,36)==1),22)));%int RED+
    FinalData(filesteller,12)= num2cell(mean(sortSingleRoi(find(sortSingleRoi(:,36)==0),22)));%int RED-
    FinalData(filesteller,13)= num2cell(mean(sortSingleRoi(:,26)));%int BLUE
    FinalData(filesteller,14)= num2cell(mean(sortSingleRoi(find(sortSingleRoi(:,37)==1),26)));%int BLUE+
    FinalData(filesteller,15)= num2cell(mean(sortSingleRoi(find(sortSingleRoi(:,37)==0),26)));%int BLUE-
    FinalData(filesteller,16)= num2cell(mean(sortSingleRoi(find(sortSingleRoi(:,36)==1 & sortSingleRoi(:,37)==0),22)));%RED+/BLUE- RED int
    FinalData(filesteller,17)= num2cell(mean(sortSingleRoi(find(sortSingleRoi(:,36)==1 & sortSingleRoi(:,37)==1),22)));%RED+/BLUE+ RED int
    FinalData(filesteller,18)= num2cell(mean(sortSingleRoi(find(sortSingleRoi(:,36)==0 & sortSingleRoi(:,37)==1),26)));%RED-/BLUE+ BLUE int
    FinalData(filesteller,19)= num2cell(mean(sortSingleRoi(find(sortSingleRoi(:,36)==1 & sortSingleRoi(:,37)==1),22)));%RED+/BLUE+ BLUE int'
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%  images                                  %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if  IwantImages==1;
        figure(1)
        sortSingleRoi=sortrows(sortSingleRoi,22);
        scatter([1:size(sortSingleRoi,1)],sortSingleRoi(:,22),'r','filled');
        title('mean intensity red channel');
        hold on
        
        figure(2)
        sortSingleRoi=sortrows(sortSingleRoi,3);
        scatter([1:size(sortSingleRoi,1)],sortSingleRoi(:,3),'g','filled');
        title('mean intensity green channel');
        hold on
        
        figure(3)
        sortSingleRoi=sortrows(sortSingleRoi,26);
        scatter([1:size(sortSingleRoi,1)],sortSingleRoi(:,26),'b','filled');
        title('mean intensity blue channel');
        hold on
        
        figure(4)
        sortSingleRoi=sortrows(sortSingleRoi,30);
        scatter([1:size(sortSingleRoi,1)],sortSingleRoi(:,30),'c','filled');
        title('mean intensity blue large channel');
        hold on
        
        figure(5)
        sortSingleRoi=sortrows(sortSingleRoi,35);
        scatter([1:size(sortSingleRoi,1)],sortSingleRoi(:,35),'k','filled');
        title('mean area puncta');
        hold on
        
        %%%plot maken van de FOS+ cellen tov het aantal puncta
        figure(6)
        scatter(FosPosData(:,22),FosPosData(:,33),'filled');
        title('FosPs vs NrPuncta');
        hold on
        
        %%%plot maken van de FOS+ cellen tov mean area puncta
        figure(7)
        scatter(FosPosData(:,22),FosPosData(:,35),'filled');
        title('FosPos vs AreaPuncta');
        hold on
    end
    
end

if  IwantImages==1;
    saveas(figure(1),strcat('Fixed_FOS.pdf'));
    saveas(figure(2),strcat('Fixed_NISSL.pdf'));
    saveas(figure(3),strcat('Fixed_PV.pdf'));
    saveas(figure(4),strcat('Fixed_PVenlargedRoi.pdf'));
    saveas(figure(5),strcat('Fixed_Puncta.pdf'));
    saveas(figure(6),strcat('Fixed_FOSvsNrPuncta.pdf'));
    saveas(figure(7),strcat('Fixed_FOSvsAreaPuncta.pdf'));
end

