%{
    function that read a star catalogue downloaded form heasarc.gsfc.nasa.gov
    and turns it into an array with 3 columns [right ascension,
    declination, magnitude] to be used in the star sensor model

    % Authors: Giuseppe Brentino
    
%}

clearvars; close all; clc;

temp_cat = readtable('star_catalogue.csv');
temp_cat = temp_cat(2:end,:);
catalogue = zeros(height(temp_cat),width(temp_cat));

% convert from from string to double

for i = 1:height(temp_cat)

    ra = table2cell(temp_cat(i,1));
    dec = table2cell(temp_cat(i,2));
    mag = table2cell(temp_cat(i,3));

    ra = replace(ra,',','.');
    dec = replace(dec,',','.');

    catalogue(i,1) = str2double(ra{1});
    catalogue(i,2) = str2double(dec{1});
    catalogue(i,3) = mag{1};
end

save star_catalogue catalogue

%% compute stars in the sensor

clearvars -except catalogue; close all; clc;

real_ra = 0;
real_deg = 0;
fov = 10; % half the field ov view
catalogue_temp = catalogue;
visible_stars = [];
for i = 1:length(catalogue(:,1))

        if real_ra-fov < 0 && catalogue(i,1) > 360-fov
            catalogue_temp(i,1) =  catalogue(i,1) - 360;
        elseif real_ra+fov > 360 && catalogue(i,1) < fov
            catalogue_temp(i,1) =  catalogue(i,1) + 360;
        end

        if real_deg-fov < - 90 
            catalogue_temp(i,2) = -90 + ( -90 - catalogue(i,2) );
        elseif real_deg+fov > 90
            catalogue_temp(i,2) = 90 - (catalogue(i,2) - 90);
        end
 

    if ( catalogue_temp(i,1) >= real_ra-fov &&  catalogue_temp(i,1) <= real_ra+fov )... %% check ra in square
            && (catalogue_temp(i,2) >= real_deg-fov &&  catalogue_temp(i,2) <= real_deg+fov)... %% check dec in square
            && norm(catalogue_temp(i,2)-real_deg,catalogue_temp(i,1)-real_ra)<=2*fov %% reduce to circle

        visible_stars = [visible_stars; catalogue(i,:) i];
    end

end

%% search star
% hp: relative ra and deg correspond to the position of the star in the sensor
% where v = -rel_deg, u = ra

i = 1; % star to look for in the catalogue