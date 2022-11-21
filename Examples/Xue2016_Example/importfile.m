function [Column1,Column2,Column3,Column4,Column5,Column6,Column7,Column8,Column9,Column10,Column11,Column12,Column13,Column14,Column15,Column16,Column17,Column18,Column19,Column20,Column21,Column22,Column23,Column24,Column25,Column26,Column27,Column28,Column29,Column30,Column31,Column32,Column33,Column34,Column35,Column36,Column37] = importfile(filename, startRow, endRow, FormatCode)

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:

% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

% import the max number of columns needed. See original data files to determine max columns needed.

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
switch FormatCode
    
  case 'S' % Well S
        DateTime = dataArray{:, 1};
        WL = dataArray{:, 2};
        Empty1 = dataArray{:, 3};
        Empty2 = dataArray{:, 4};
        Empty3 = dataArray{:, 5};
        Empty4 = dataArray{:, 6};
        Empty5 = dataArray{:, 7};
        Empty6 = dataArray{:, 8};
        Empty7 = dataArray{:, 9};
        Empty8 = dataArray{:, 10};
        Empty9 = dataArray{:, 11};
        Empty10 = dataArray{:, 12};
        Empty11 = dataArray{:, 13};
        Empty12 = dataArray{:, 14};
        Empty13 = dataArray{:, 15};
        Empty14 = dataArray{:, 16};
        Empty15 = dataArray{:, 17};
        Empty16 = dataArray{:, 18};
        Empty17 = dataArray{:, 19};
        Empty18 = dataArray{:, 20};
        Empty19 = dataArray{:, 21};
        Empty20 = dataArray{:, 22};
        Empty21 = dataArray{:, 23};
        Empty22 = dataArray{:, 24};
        Empty23 = dataArray{:, 25};
        Empty24 = dataArray{:, 26};
        Empty25 = dataArray{:, 27};
        Empty26 = dataArray{:, 28};
        Empty27 = dataArray{:, 29};
        Empty28 = dataArray{:, 30};
        Empty29 = dataArray{:, 31};
        Empty30 = dataArray{:, 32};
        Empty31 = dataArray{:, 33};
        Empty32 = dataArray{:, 34};
        Empty33 = dataArray{:, 35};
        Empty34 = dataArray{:, 36};
        Empty35 = dataArray{:, 37};
        
        Column1 = DateTime;
        Column2 = WL;
         Column3 = Empty1;
        Column4 = Empty2;
        Column5 = Empty3;
        Column6 = Empty4;
        Column7 = Empty5;
        Column8 = Empty6;
        Column9 = Empty7;
        Column10 = Empty8;
        Column11 = Empty9;
        Column12 = Empty10;
        Column13 = Empty11;
        Column14 = Empty12;
        Column15 = Empty13;
        Column16 = Empty14;
        Column17 = Empty15;
        Column18 = Empty16;
        Column19 = Empty17;
        Column20 = Empty18;
        Column21 = Empty19;
        Column22 = Empty20;
        Column23 = Empty21;
        Column24 = Empty22;
        Column25 = Empty23;
        Column26 = Empty24;
        Column27 = Empty25;
        Column28 = Empty26;
        Column29 = Empty27;
        Column30 = Empty28;
        Column31 = Empty29;
        Column32 = Empty30;
        Column33 = Empty31;
        Column34 = Empty32;
        Column35 = Empty33;
        Column36 = Empty34;
        Column37 = Empty35;
        
case 'A' % Well A
        DateTime = dataArray{:, 1};
        ETmin = dataArray{:, 2};
        EThours = dataArray{:, 3};
        ETdays = dataArray{:, 4};
        dtw = dataArray{:, 5};
        drawdown = dataArray{:, 6};
        GWElev = dataArray{:, 7};
        Empty1 = dataArray{:, 8};
        Empty2 = dataArray{:, 9};
        Empty3 = dataArray{:, 10};
        Empty4 = dataArray{:, 11};
        Empty5 = dataArray{:, 12};
        Empty6 = dataArray{:, 13};
        Empty7 = dataArray{:, 14};
        Empty8 = dataArray{:, 15};
        Empty9 = dataArray{:, 16};
        Empty10 = dataArray{:, 17};
        Empty11 = dataArray{:, 18};
        Empty12 = dataArray{:, 19};
        Empty13 = dataArray{:, 20};
        Empty14 = dataArray{:, 21};
        Empty15 = dataArray{:, 22};
        Empty16 = dataArray{:, 23};
        Empty17 = dataArray{:, 24};
        Empty18 = dataArray{:, 25};
        Empty19 = dataArray{:, 26};
        Empty20 = dataArray{:, 27};
        Empty21 = dataArray{:, 28};
        Empty22 = dataArray{:, 29};
        Empty23 = dataArray{:, 30};
        Empty24 = dataArray{:, 31};
        Empty25 = dataArray{:, 32};
        Empty26 = dataArray{:, 33};
        Empty27 = dataArray{:, 34};
        Empty28 = dataArray{:, 35};
        Empty29 = dataArray{:, 36};
        Empty30 = dataArray{:, 37};
        
        Column1 = DateTime;
        Column2 = ETmin;
        Column3 = EThours;
        Column4 = ETdays;
        Column5 = dtw;
        Column6 = drawdown;
        Column7 = GWElev;
        Column8 = Empty1;
        Column9 = Empty2;
        Column10 = Empty3;
        Column11 = Empty4;
        Column12 = Empty5;
        Column13 = Empty6;
        Column14 = Empty7;
        Column15 = Empty8;
        Column16 = Empty9;
        Column17 = Empty10;
        Column18 = Empty11;
        Column19 = Empty12;
        Column20 = Empty13;
        Column21 = Empty14;
        Column22 = Empty15;
        Column23 = Empty16;
        Column24 = Empty17;
        Column25 = Empty18;
        Column26 = Empty19;
        Column27 = Empty20;
        Column28 = Empty21;
        Column29 = Empty22;
        Column30 = Empty23;
        Column31 = Empty24;
        Column32 = Empty25;
        Column33 = Empty26;
        Column34 = Empty27;
        Column35 = Empty28;
        Column36 = Empty29;
        Column37 = Empty30;
        
        
     case 'B' %Well B
        DateTime = dataArray{:, 1};
        ETmin = dataArray{:, 2};
        EThours = dataArray{:, 3};
        ETdays = dataArray{:, 4};
        FeetH2OPorta = dataArray{:, 5};
        FeetH2OPortb = dataArray{:, 6};
        FeetH2OPortc = dataArray{:, 7};
        FeetH2OPortd = dataArray{:, 8};
        FeetH2OPorte = dataArray{:, 9};
        FeetH2OPortf = dataArray{:, 10};
        FeetH2OPortg = dataArray{:, 11};
        FeetH2OPorth = dataArray{:, 12};
        InchHg = dataArray{:, 13};
        dtwPorta = dataArray{:, 14};
        dtwPortb = dataArray{:, 15};
        dtwPortc = dataArray{:, 16};
        dtwPortd = dataArray{:, 17};
        dtwPorte = dataArray{:, 18};
        dtwPortf = dataArray{:, 19};
        dtwPortg = dataArray{:, 20};
        dtwPorth = dataArray{:, 21};
        drawdownPorta = dataArray{:, 22};
        drawdownPortb = dataArray{:, 23};
        drawdownPortc = dataArray{:, 24};
        drawdownPortd = dataArray{:, 25};
        drawdownPorte = dataArray{:, 26};
        drawdownPortf = dataArray{:, 27};
        drawdownPortg = dataArray{:, 28};
        drawdownPorth = dataArray{:, 29};
        GWElevPorta = dataArray{:, 30};
        GWElevPortb = dataArray{:, 31};
        GWElevPortc = dataArray{:, 32};
        GWElevPortd = dataArray{:, 33};
        GWElevPorte = dataArray{:, 34};
        GWElevPortf = dataArray{:, 35};
        GWElevPortg = dataArray{:, 36};
        GWElevPorth = dataArray{:, 37};
        
        Column1 = DateTime;
        Column2 = ETmin;
        Column3 = EThours;
        Column4 = ETdays;
        Column5 = FeetH2OPorta;
        Column6 = FeetH2OPortb;
        Column7 = FeetH2OPortc;
        Column8 = FeetH2OPortd;
        Column9 = FeetH2OPorte;
        Column10 = FeetH2OPortf;
        Column11 = FeetH2OPortg;
        Column12 = FeetH2OPorth;
        Column13 = InchHg;
        Column14 = dtwPorta;
        Column15 = dtwPortb;
        Column16 = dtwPortc;
        Column17 = dtwPortd;
        Column18 = dtwPorte;
        Column19 = dtwPortf;
        Column20 = dtwPortg;
        Column21 = dtwPorth;
        Column22 = drawdownPorta;
        Column23 = drawdownPortb;
        Column24 = drawdownPortc;
        Column25 = drawdownPortd;
        Column26 = drawdownPorte;
        Column27 = drawdownPortf;
        Column28 = drawdownPortg;
        Column29 = drawdownPorth;
        Column30 = GWElevPorta;
        Column31 = GWElevPortb;
        Column32 = GWElevPortc;
        Column33 = GWElevPortd;
        Column34 = GWElevPorte;
        Column35 = GWElevPortf;
        Column36 = GWElevPortg;
        Column37 = GWElevPorth; 
end

