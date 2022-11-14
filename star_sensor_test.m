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

    catalogue(i,1) = str2double(ra{1});
    catalogue(i,2) = str2double(dec{1});
    catalogue(i,3) = mag{1};
end

save star_catalogue catalogue

%% compute stars in the sensor

clearvars -except catalogue; close all; clc;

real_ra = 0;
real_deg = 90;
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

mag_toll = 5e-2;
deg_toll = 2/3600;
ra_toll = 2/3600;

ref_star = visible_stars;          % reference star in the sensor

%CHECK

for i = 1:size(ref_star,1)
    if ref_star(i) > 360 - fov
        ref_star(i) = ref_star(i)-360;
    end
end

u_star = ref_star(:,1) - real_ra;          % ref star posizion in sensor x axis
v_star = real_deg - ref_star(:,2);         % ref star posizion in sensor y axis


flag.find_mag = 1;          % 1 if i want to compare the reference magnitude with the catalogue stars
flag.find_pos = 0;          %  % 1 if i want to check the geometry
flag.found = 0;
k=1;                        % counter for geometry check
count_geom = 2;

while k < size(ref_star,1)
    for i = 1:size(catalogue,1)

        % compare the magnitude of the reference star with the ones in the
        % catalogue
        if catalogue(i,3) >= ref_star(1,3)-mag_toll && ...
                catalogue(i,3) <= ref_star(1,3)+mag_toll && flag.find_mag == 1

            flag.find_mag = 0;
            flag.find_pos = 1;
            possible_match{k,1} = [catalogue(i,:) i];
        end

        %check relative geometry
        if flag.find_pos == 1

            rel_pos = [u_star(k+1)-u_star(k);v_star(k+1)-v_star(k)];

            % look for a star in the catalogue with the same relative position
            % wrt the star with index i
            for j = 1:size(catalogue,1)

                if catalogue(j,2) >= catalogue(i,2)+rel_pos(2)-deg_toll && ...
                        catalogue(j,2) <= catalogue(i,2)+rel_pos(2)+deg_toll && ...
                        catalogue(j,3) >= ref_star(k+1,3)-mag_toll && ...
                        catalogue(j,3) <= ref_star(k+1,3)+mag_toll
                
                    % check ra
                    diff_ra = catalogue(i,1)+rel_pos(1);

                    if  diff_ra > 360
                        diff_ra = catalogue(i,1)+rel_pos(1) -360;
                    elseif diff_ra < - 360
                          diff_ra = catalogue(i,1)+rel_pos(1) + 360;
                    end

                    if catalogue(j,1) >= diff_ra-ra_toll && ...
                    catalogue(j,1) <= diff_ra+ra_toll

                    flag.found = 1;

                    possible_match{k,count_geom} = [catalogue(j,:) i];
                    end

                end

            end

        end

    end
            k = k+1;
end
%% use later to resize the catalogue at each time step to improve computational speed

figure()
hold on
grid on
c = visible_stars(:,3);
scatter(u_star,v_star,[],c,'filled');
colorbar('Direction','reverse')

% % % x = 1:size(catalogue,1);
% % % p = polyfit(x,catalogue(:,2),1);
% % % val = polyval(p,x);
% % % 
% % % plot(x,val)

