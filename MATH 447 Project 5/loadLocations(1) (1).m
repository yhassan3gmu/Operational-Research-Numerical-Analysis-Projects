%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/tyrusberry/Dropbox/Teaching/447/Projects/COVID19PDE/locations.csv
%
% Auto-generated by MATLAB on 22-Apr-2020 21:40:05

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["countriesAndTerritories", "location", "continent", "population_year", "population"];
opts.VariableTypes = ["string", "string", "categorical", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["countriesAndTerritories", "location"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["countriesAndTerritories", "location", "continent"], "EmptyFieldRule", "auto");

% Import the data
locations = readtable("locations.csv", opts);


%% Clear temporary variables
clear opts

