%{
    function that read a star catalogue downloaded form heasarc.gsfc.nasa.gov
    and turns it into an array with 3 columns [right ascension,
    declination, magnitude] to be used in the star sensor model

    % Authors: Giuseppe Brentino
    
%}

clearvars; close all; clc;

temp_cat = readtable('star_catalogue.csv');
temp_cat = temp_cat(2:end,:);
%temp_cat = temp_cat(1:3:height(temp_cat),:); % reduce the catalogue to 1/3 of its original length
catalogue = zeros(height(temp_cat),width(temp_cat));

% convert from from string to double

for i = 1:height(temp_cat)

    ra = table2cell(temp_cat(i,1));
    dec = table2cell(temp_cat(i,2));
    mag = table2cell(temp_cat(i,3));

    ra = replace(ra,',','.');
    dec = replace(dec,',','.');

    Ra = str2double(ra{1});
    Dec = str2double(dec{1});
    
    % star position versor in ECI
    catalogue(i,1) = cos(Dec)*cos(Ra);
    catalogue(i,2) = cos(Dec)*sin(Ra);
    catalogue(i,3) = sin(Dec);
    % star magnitude
    catalogue(i,4) = mag{1};
    
end

save .\data\star_catalogue catalogue

