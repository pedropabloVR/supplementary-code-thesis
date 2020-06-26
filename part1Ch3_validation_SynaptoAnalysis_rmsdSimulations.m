% Supplementary code accompanying Chapter 3 of PVR thesis to generate the
% plots in Figs. 3.15 and 3.16.  

% Simulations to model the sensitivites of the root mean square distance
% from centroid (RMSD) as a measure of cluster size. 

% Pedro Vallejo Ramirez
% Created: 20/11/2019
% Updated: 24/06/2020

close all 
clear all

% Simulation parameters

R           = 100; % circle radius in nm
R_array     = 100:100:1000; % nm
locs        = 50000; % ball park number of localisations for a dense blob (change from 1000 to 50000 for results in thesis)
center      = R; % place the center at the same distance of the radius 
jitter_max  = 100; % nm
jitter_increments = 0:10:jitter_max;
jitter_arr  = randi([0 jitter_max],1,locs); % simulate an array of uncertainty errors  
show_circ   = 0; 
show_square = 0;
save        = 0;
vertices    = 5;
show_masks  = 0;

% Variables for making nice figures in matlab

width               = 10;  % width in centimeters
height              = 10;  % height in centimeters
font_size           = 24;
lw                  = 1.0; %linewidth
msz                 = 4;  %markersize

Path = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd simulations';

%% Simulations with real masks from synaptosome experiments
% Code to generate simulation shown in Fig. 3.16 (b) in thesis 

% Import masks 20191120

PathMask    = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/4Cb_phys_a1_mask_red.tif';
im          = imread(PathMask);
mask        = (im>0); %convert into binary mask
stats       = regionprops(mask,'Area','Centroid','BoundingBox');
imlabel     = bwlabel(mask);
vislabels(imlabel); % pick your favorite blob

RMSD_blob = zeros(1,size(stats,1));
RMSD_jitter = zeros(1,size(stats,1));
error_blob = zeros(1,size(stats,1));


for i = 1:size(stats,1)
    % extract blob 20 with a rectangle of ~100 pixels around it
    blob_stats  = stats(i);
    %expand bounding box by twice its largest dimension
    if blob_stats.BoundingBox(4) > blob_stats.BoundingBox(3)
        width       = blob_stats.BoundingBox(4)*2;
    else 
        width       = blob_stats.BoundingBox(3)*2;   
    end 
    crop_rect   = [blob_stats.Centroid(1)-width/2 blob_stats.Centroid(2)-width/2 width width];

    mask_crop   = imcrop(mask,crop_rect); % obtain cropped blob
    stat_crop   = regionprops(mask_crop);

    % now fill it up with localisations 
    numPoints   = locs;
    center      = width/2;
    xRandom     = center+center* rand(1, numPoints) - width / 4;
    yRandom     = center+center* rand(1, numPoints) - width / 4;
    % all these points lie close to the bounding box of the blob
% 
%     figure
%     imshow(mask_crop)
%     hold on;
%     scatter(yRandom,xRandom);

    % discard points which do not fall inside the mask
    in_points   = mask_crop(ceil(xRandom(:)),ceil(yRandom(:)));
    diagonal    = diag(in_points); % the diagonal of the matrix contains the indices of all points inside the blobs
    xRandom_in  = xRandom(diagonal);
    yRandom_in  = yRandom(diagonal);
    locs_placed = sum(diagonal); % number of localisations actually placed inside the mask

    % jitter the localisations inside the mask
    jitter_sign = getSign(locs_placed); % generate random 1 and -1 array
    t           = 2*pi*rand(1,locs_placed); % create angle array
    jitter_array  = jitter_arr(diagonal); %select values corresponding to the x,y points inside mask
    x_jitter    = xRandom_in+jitter_array.*jitter_sign.*cos(t); 
    y_jitter    = yRandom_in+jitter_array.*jitter_sign.*sin(t);
    
    if show_masks
        figure
        imshow(mask_crop,'InitialMagnification', 1500)
        title(['Blob id: ' num2str(i) ' Area: ' num2str(blob_stats.Area) ' Locs: ' num2str(locs_placed) ' Max Jitter: ' num2str(jitter_max)]);
        hold on
        scatter(yRandom_in,xRandom_in);
        scatter(y_jitter,x_jitter,2,'+','r');
         set(gca,...
        'FontSize',font_size,...
        'Units','Normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontName','Arial');
    end 
    h = gcf;
    % now compute the RMSD for the unjittered and the jittered spots
    center_blob =size(mask_crop,2)/2;
    RMSD_blob(i) = sqrt(mean( (xRandom_in - center_blob).^2 + (yRandom_in - center_blob).^2));
    RMSD_jitter(i) = sqrt(mean( (x_jitter - center_blob).^2 + (y_jitter - center_blob).^2));
    error_blob(i) = (abs(RMSD_blob(i) - RMSD_jitter(i))/RMSD_blob(i)) *100;
    area_blob(i) =  blob_stats.Area;  
    %make_animation( h,i,[Path filesep '_jitterMax_10_10000locs.gif'])
    %pause(0.2) %you can enter the time in pause to change the loop
    
end 
% figure;scatter(area_blob,error_blob);title('Area vs Error in RMSD after jitter');
% xlabel('Area (pixels)');
% ylabel('RMSD error');


figure;scatter(area_blob.*(11.7^2),error_blob,'+','r' );title('Area vs Error in RMSD after jitter');
xlabel('Area (nm^2)');
ylabel('RMSD error');
         set(gca,...
        'FontSize',font_size,...
        'Units','Normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontName','Arial');
fig1 = gcf;

if save  
    print(fig1,[Path filesep '_AreaVsRMSDError_afterJitter'],'-dpng','-r300') % print at 300 dpi resolution
end

%% Record a single blob being filled with localisations and compute RMSD vs amount of increasing jitter
% Code to generate simulation shown in Fig. 3.16 (a) in thesis 

RMSD_blob = zeros(1,size(jitter_increments,2));
RMSD_jitter = zeros(1,size(jitter_increments,2));
error_blob = zeros(1,size(jitter_increments,2));

% iterate over different amounts of jitter

for i = 1:size(jitter_increments,2)
    % extract blob 20 with a rectangle of ~100 pixels around it
    blob_stats  = stats(20);
    %expand bounding box by twice its largest dimension
    if blob_stats.BoundingBox(4) > blob_stats.BoundingBox(3)
        width       = blob_stats.BoundingBox(4)*2;
    else 
        width       = blob_stats.BoundingBox(3)*2;   
    end 
    crop_rect   = [blob_stats.Centroid(1)-width/2 blob_stats.Centroid(2)-width/2 width width];

    mask_crop   = imcrop(mask,crop_rect); % obtain cropped blob
    stat_crop   = regionprops(mask_crop);


    % now fill it up with localisations 
    numPoints   = locs;
    center      = width/2;
    xRandom     = center+center* rand(1, numPoints) - width / 4;
    yRandom     = center+center* rand(1, numPoints) - width / 4;
    % all these points lie close to the bounding box of the blob
% 
%     figure
%     imshow(mask_crop)
%     hold on;
%     scatter(yRandom,xRandom);

    % discard points which do not fall inside the mask
    in_points   = mask_crop(ceil(xRandom(:)),ceil(yRandom(:)));
    diagonal    = diag(in_points); % the diagonal of the matrix contains the indices of all points inside the blobs
    xRandom_in  = xRandom(diagonal);
    yRandom_in  = yRandom(diagonal);
    locs_placed = sum(diagonal); % number of localisations actually placed inside the mask

    % jitter the localisations inside the mask
    jitter_inc  = jitter_increments(i);
    jitter_arr  = randi([0 jitter_inc],1,locs_placed); % simulate an array of uncertainty errors  
    jitter_sign = getSign(locs_placed); % generate random 1 and -1 array
    t           = 2*pi*rand(1,locs_placed); % create angle array
    x_jitter    = xRandom_in+jitter_arr.*jitter_sign.*cos(t); 
    y_jitter    = yRandom_in+jitter_arr.*jitter_sign.*sin(t);
    

    if show_masks
        figure
        imshow(mask_crop) %'InitialMagnification', 800)
        title(['Blob id: ' num2str(20) ' Area: ' num2str(blob_stats.Area) ' Locs: ' num2str(locs_placed) ' Max Jitter: ' num2str(jitter_inc)]);
        hold on
        scatter(yRandom_in,xRandom_in);
        scatter(y_jitter,x_jitter,2,'+','r');
         set(gca,...
        'FontSize',font_size,...
        'Units','Normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontName','Arial');
     truesize([800 800])
 
        %xlabel(['Jitter: ' num2str(jitter_inc) ' Area: ' num2str(blob_stats.Area)])
    end 
    h = gcf;
    % now compute the RMSD for the unjittered and the jittered spots
    center_blob =size(mask_crop,2)/2;
    RMSD_blob(i) = sqrt(mean( (xRandom_in - center_blob).^2 + (yRandom_in - center_blob).^2));
    RMSD_jitter(i) = sqrt(mean( (x_jitter - center_blob).^2 + (y_jitter - center_blob).^2));
    error_blob(i) = (abs(RMSD_blob(i) - RMSD_jitter(i))/RMSD_blob(i)) *100;
%     make_animation( h,i,[Path filesep '_JiitterIncrements_20Max_10000locs.gif'])
%     pause(0.2) %you can enter the time in pause to change the loop
%     
end 

figure;scatter(jitter_increments.*11.7,error_blob);title('Jitter vs Error in RMSD');
xlabel('Jitter (nm)');
ylabel('RMSD error');
         set(gca,...
        'FontSize',font_size,...
        'Units','Normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontName','Arial');
fig1 = gcf;

if save  
    print(fig1,[Path filesep '_JitterVsRMSDError'],'-dpng','-r300') % print at 300 dpi resolution
end

%% Modeling a compact vs disperse distribution of points within the circle
% Code for generating Figure 3.15 in thesis 

% testing the sensitivity of the RMSD measure to the relative spread of
% points inside a ground truth circle

RMSD_jitter_circ = zeros(1,size(jitter_increments,2));
error_jitter_circ = zeros(1,size(jitter_increments,2));

for i = 1:size(jitter_increments,2)
    
    % randomly place localisations (x,y, intensity) inside the boundaries of
    % the circle, jittered by an amount sigma equal to the uncertainty of the 
    % localisation precision
    
    R           = R;
    center      = R;
    t           = 2*pi*rand(1,locs); % create angle array
    r           = R*sqrt(rand(1,locs)); % array of possible radii inside the circle  
    jitter_inc  = jitter_increments(i);
    jitter_arr  = randi([0 jitter_inc],1,locs); % simulate an array of uncertainty errors  
    jitter_sign = getSign(locs);
    x           = center + r.*cos(t)+jitter_arr.*jitter_sign.*cos(t); % add points randomly inside circle, jittered around the circle
    y           = center + r.*sin(t)+jitter_arr.*jitter_sign.*sin(t); 
    
    if show_circ
        figure;
        scatter(x,y,'+');
        hold on;
        [h,x_circ,y_circ] = makeCircle(R,R,R); % draw theoretical circle with localisations inside
        xlim([-100,300])
        ylim([-100,300])
        axis equal
         set(gca,...
                'FontSize',font_size,...
                'Units','Normalized',...
                'FontUnits','points',...
                'FontWeight','normal',...
                'FontName','Arial');
        title(['Circle of radius ' num2str(R) '. Max jitter: ' num2str(jitter_inc) ]);
    end     
     h = gcf;
    % calculate the RMSD of the spatial point pattern from this and compare to
    % the ground truth size of the circle. 

    % Get mean center of mass
    Xc = mean(x);
    Yc = mean(y);
    R_mean = 2*R/3;
    % Calculate root mean squared distance
    RMSD_jitter_circ(i) = sqrt(mean((x - Xc).^2 + (y - Yc).^2));
    error_jitter_circ(i) = (abs(R_mean - RMSD_jitter_circ(i))/R_mean) *100; % error between the average distance from the center of the circle and the calculated RMSD from the spread of points
    fprintf(" radius: %0.2f nm, 2R/3:%0.2f nm, Num. points %0.2f, RMSD: %2.2f nm, Percentage error: %2.2f%% \n", R,R_mean,locs, RMSD_jitter_circ(i),error_jitter_circ(i));
    if show_circ
        make_animation( h,i,[Path filesep 'radius_100_jitterMax_100_10000locs.gif'])
        pause(0.2) %you can enter the time in pause to change the loop
    end
end

%% Plot the error as a function of the theoretical radius (Fig. 3.15 in thesis)
figure;
plot(jitter_increments,error_jitter_circ);
xlabel('Maximum jitter (nm)');
ylabel('RMSD vs theory % error')
ylim([0 50])
 set(gca,...
                'FontSize',font_size,...
                'Units','Normalized',...
                'FontUnits','points',...
                'FontWeight','normal',...
                'FontName','Arial');
title(['Number of points: ' num2str(locs) ]);
fig1 = gcf;

if save  
    print(fig1,[Path filesep '_percentageError_Jitteredcircle'],'-dpng','-r300') % print at 300 dpi resolution
end


%% Obsolete code 2019 11 21

%% iterate over different radii values to calculate the mean difference
% betweent the theoretical radius and the RMSD calculation
% RMSD = zeros(1,size(R_array,2));
% error = zeros(1,size(R_array,2));
% RMSD_square = zeros(1,size(R_array,2));
% error_square = zeros(1,size(R_array,2));

% for i = 1:size(R_array,2)
%     
%     % randomly place localisations (x,y, intensity) inside the boundaries of
%     % the circle, jittered by an amount sigma equal to the uncertainty of the 
%     % localisation precision
%     
%     R           = R_array(i);
%     center      = R;
%     t           = 2*pi*rand(1,locs); % create angle array
%     r           = R*sqrt(rand(1,locs)); % array of possible radii inside the circle  
%     jitter_sign = getSign(locs);
%     x           = center + r.*cos(t)+jitter_arr.*jitter_sign.*cos(t); % add points randomly inside circle, jittered around the circle
%     y           = center + r.*sin(t)+jitter_arr.*jitter_sign.*sin(t); 
%     
%     if show_circ
%         figure;
%         scatter(x,y,'+');
%         hold on;
%         [h,x_circ,y_circ] = makeCircle(R,R,R); % draw theoretical circle with localisations inside
%         axis equal
%     end     
%     
%     % calculate the RMSD of the spatial point pattern from this and compare to
%     % the ground truth size of the circle. 
% 
%     % Get mean center of mass
%     Xc = mean(x);
%     Yc = mean(y);
%     R_mean = 2*R/3;
%     % Calculate root mean squared distance
%     RMSD(i) = sqrt(mean((x - Xc).^2 + (y - Yc).^2));
%     error(i) = (abs(R_mean - RMSD(i))/R_mean) *100; % error between the average distance from the center of the circle and the calculated RMSD from the spread of points
%     fprintf(" radius: %0.2f nm, 2R/3:%0.2f nm, RMSD: %2.2f nm, Percentage error: %2.2f%% \n", R,R_mean, RMSD(i),error(i));
% 
%     %% now try for a square
%     % place in a loop to try for different side lengths of squares and
%     % introduce jitter. 
% 
%     width = 2*R; % same width as diameter of circle
%     x = 0;
%     y = 0; % centered at the origin
%     % fill the bounding rectangle with localisations
%     x_square_rand = center+width * rand(1, locs) - width / 2+jitter_arr.*jitter_sign.*cos(t);
%     y_square_rand = center+width * rand(1, locs) - width / 2+jitter_arr.*jitter_sign.*sin(t);
%     if show_square
%         figure
%         rectangle('Position',[x y width width]) % draw rectangle around bounding box
%         grid on;
%         hold on;
%         plot(x_square_rand, y_square_rand, '+', 'MarkerSize', 10);
%         xlabel('X');
%         ylabel('Y');
%         axis equal;
%     end 
% 
%     x_square_mean = mean(x_square_rand);
%     y_square_mean = mean(y_square_rand);
% 
%     % formula for the average distance of any point inside a square to the center of the square
%     avg_dist_inside_square = (width/6)*(sqrt(2)+log(1+sqrt(2)));  
% 
%     % get RMSD from perfectly filled square
%     RMSD_square(i) = sqrt(mean( (x_square_rand - center).^2 + (y_square_rand - center).^2));
%     error_square(i) = (abs(avg_dist_inside_square - RMSD_square(i))/avg_dist_inside_square) *100;
% 
% end
% 
% % Plot the error as a function of the theoretical radius
% figure;
% plot(R_array,error);
% xlabel('Radius of theoretical circle (nm))');
% ylabel('Percentage error of the RMSD vs theory')
%  set(gca,...
%                 'FontSize',font_size,...
%                 'Units','Normalized',...
%                 'FontUnits','points',...
%                 'FontWeight','normal',...
%                 'FontName','Arial');
% axis equal
% fig1 = gcf;
% if save  
%     Path = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/synaptosomes/rmsd simulations';
%     print(fig1,[Path filesep '_percentageError_circle'],'-dpng','-r300') % print at 300 dpi resolution
% end
% 
% figure;
% plot(R_array*2,error_square);
% xlabel('Side length of theoretical square (nm))');
% ylabel('Percentage error of the RMSD vs theory')
%  set(gca,...
%                 'FontSize',font_size,...
%                 'Units','Normalized',...
%                 'FontUnits','points',...
%                 'FontWeight','normal',...
%                 'FontName','Arial');


%% Deprecated code 2019 11 01

% Generalizing the simulation for a polygon with N sides

% % Specify polygon variables
% numVert = 5; % vertices
% radius = R; % radius (same as circle)
% radVar = 1; % variance in the spikiness of vertices
% angVar = 1; % variance in spacing of vertices around the unit circle
% 
% % create the matrix of angles for equal separation of points
% angleMatrix = 0: 2*pi/numVert: 2*pi;
% 
% coeff  = 0.5; %use a coefficient instead of a random variable for generating the polygon
% % use a vector instead of a for loop to generate the polygon bounds
% vert_array = 1:numVert;
% x_poly = center+(radius + radius*rand(1)*radVar) * cos(angleMatrix + (2*pi/numVert)*rand(1)*angVar);
% y_poly = center+(radius + radius*rand(1)*radVar) * sin(angleMatrix + (2*pi/numVert)*rand(1)*angVar);
% 
% % Graph the polygon and connect the final point to the first point
% figure
% plot([x_poly'; x_poly(1)],[y_poly'; y_poly(1)],'ro-')
% hold on
% pgon = polyshape(x_poly,y_poly); % fit a matlab polygon to the randomly generated vertices
% plot(pgon); 
% [xlim,ylim] = boundingbox(pgon);
% plot(xlim,ylim,'r*',xlim,fliplr(ylim),'g*')
% size_x = round(xlim(2)-xlim(1));
% size_y = round(ylim(2)-ylim(1));
% 
% % center of the bounding box
% 
% % Use the bounding box to generate a mask (matrix) with 1s inside the
% % polygon and zeros outside
% n = round(xlim(2)-xlim(1));
% m = round(ylim(2)-ylim(1));
% mask_pgon = poly2mask(x_poly,y_poly,m+center,n+center);
% figure
% imshow(flipud(mask_pgon))
% hold on
% plot(x_poly,y_poly,'b','LineWidth',2)
% plot([x_poly'; x_poly(1)],[y_poly'; y_poly(1)],'ro-')
% 
% hold off
% 
% % alternatively, use a matlab in-built function to generate a regular polygon
% pgon2 = nsidedpoly(numVert,'Center',[center center],'Radius',R);
% 
% % Now fill the polygon with jittered points
% % need to generate a set of random radius and angle values within the boundaries of the
% % polygon
% 
% t           = 2*pi*rand(1,locs); % create angle array
% r           = (radius + radius*coeff*radVar)*sqrt(rand(1,locs)); % array of possible radii inside the circle  
% jitter_sign = getSign(locs);
% 
% %x_polyIn = center+(radius + radius*rand(1)*radVar) * cos(angleMatrix(k) + (2*pi/numVert)*rand(1)*angVar)+jitter_arr.*jitter_sign.*cos(t);
% %y_polyIn = center+(radius + radius*rand(1)*radVar) * sin(angleMatrix(k) + (2*pi/numVert)*rand(1)*angVar)+jitter_arr.*jitter_sign.*cos(t);
% theta = (t + (2*pi/numVert)*rand(1)*angVar);
% x_polyIn = center+(r + r.*rand(1)*radVar).*cos(theta); %+jitter_arr.*jitter_sign.*cos(theta);
% y_polyIn = center+(r + r.*rand(1)*radVar).*sin(theta); %+jitter_arr.*jitter_sign.*cos(theta);
% scatter(x_polyIn,y_polyIn)
% 
% 
% x           = center + r.*cos(t)+jitter_arr.*jitter_sign.*cos(t); % add points randomly inside circle, jittered around the circle
% y           = center + r.*sin(t)+jitter_arr.*jitter_sign.*sin(t); 


%% Solution with a single square
%  
% % place in a loop to try for different side lengths of squares and
% % introduce jitter. 
% 
% width = 2*R; % same width as diameter of circle
% x = 0;
% y = 0; % centered at the origin
% % fill the bounding rectangle with localisations
% figure
% rectangle('Position',[x y width width]) % draw rectangle around bounding box
% grid on;
% x_square_rand = center+width * rand(1, locs) - width / 2;
% y_square_rand = center+width * rand(1, locs) - width / 2;
% hold on;
% plot(x_square_rand, y_square_rand, '+', 'MarkerSize', 10);
% xlabel('X');
% ylabel('Y');
% axis equal;
% 
% x_square_mean = mean(x_square_rand);
% y_square_mean = mean(y_square_rand);
% 
% % formula for the average distance of any point inside a square to the center of the square
% avg_dist_inside_square = (width/6)*(sqrt(2)+log(1+sqrt(2)));  
% 
% % get RMSD from perfectly filled square
% RMSD_square = sqrt(mean( (x_square_rand - center).^2 + (y_square_rand - center).^2));
% error = ((avg_dist_inside_square - RMSD)/avg_dist_inside_square) *100;
% %fprintf(" side length: %0.2f nm, RMSD: %2.2f nm, Percentage error: %2.2f%% \n", R, RMSD(i),error(i));
% 
% 
% %% Get average distance from any point in the square perimeter to the center of the square
% % define vector along half the vertices of the square
% x_vert = x:10:width/2;
% vert_x = [x_vert;zeros(1,size(x_vert,2))];
% y_vert = y:10:width/2;
% vert_y = [y_vert;zeros(1,size(y_vert,2))];
% 
% avg_vert_dist = sqrt((vert_x(1,:)-center).^2 +(vert_x(2,:)-center).^2);
% % this is the average distance from any point along the perimeter of the square to the center of the square 
% rms_vert_dist = mean(avg_vert_dist);


function make_animation( h,index,filename )
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if index == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end

end 