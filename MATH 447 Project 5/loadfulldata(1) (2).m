%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/tyrusberry/Downloads/full_data.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2020/03/27 20:27:38
clear;clc;
loadLocations;

%% Initialize variables.
filename = 'full_data.csv';
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: datetimes (%{yyyy-MM-dd}D)
%	column2: categorical (%C)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%{yyyy-MM-dd}D%C%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', 0.0, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
fulldata = table(dataArray{1:end-1}, 'VariableNames', {'date','location','new_cases','new_deaths','total_cases','total_deaths'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% fulldata.date=datenum(fulldata.date);

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

Data = zeros(size(fulldata));
[a,b,Data(:,2)] = unique(fulldata(:,2));

for i=1:size(fulldata,1)
    Data(i,1) = datenum(fulldata{i,1});
end

Data(:,1) = Data(:,1)-min(Data(:,1))+1;
T = max(Data(:,1));

Data(:,3:6) = fulldata{:,3:6};
DataMatrix = Data;

N = size(a,1);

P = zeros(N,1);         %%% Initial susceptible population for each country
Continent = zeros(N,1); %%% contenent label for each country
AllCases = zeros(N,T);  %%% Cases for each country
AllDeaths = zeros(N,T); %%% Deaths for each country


validLocations = 0;
[c,d,continents] = unique(locations(:,3));



for i = 1:N
    for j = 1:size(locations,1)
        if (a{i,1} == locations{j,2})&&(~isnan(locations{j,5}))
            validLocations = validLocations+1;
            Countries(validLocations) = a{i,1};
            P(validLocations) = locations{j,5};
            Continent(validLocations) = continents(j);
            Dates = Data(Data(:,2)==i,1);
            NewCases = zeros(T,1);
            NewCases(Dates) = Data(Data(:,2)==i,3);
            NewDeaths = zeros(T,1);
            NewDeaths(Dates) = Data(Data(:,2)==i,4);
            AllCases(validLocations,:) = cumsum(NewCases);
            AllDeaths(validLocations,:) = cumsum(NewDeaths); 
        end
    end
end

N=validLocations;
P = P(1:N);
Continent = Continent(1:N);
AllCases = AllCases(1:N,:);
AllDeaths = AllDeaths(1:N,:);


clear DataMatrix NewDeaths b fulldata validLocations Dates c i continents j Data NewCases a d locations

