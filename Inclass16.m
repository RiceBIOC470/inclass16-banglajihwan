% Inclass16

%The folder in this repository contains code implementing a Tracking
%algorithm to match cells (or anything else) between successive frames. 
% It is an implemenation of the algorithm described in this paper: 
%
% Sbalzarini IF, Koumoutsakos P (2005) Feature point tracking and trajectory analysis 
% for video imaging in cell biology. J Struct Biol 151:182?195.
%
%The main function for the code is called MatchFrames.m and it takes three
%arguments: 
% 1. A cell array of data called peaks. Each entry of peaks is data for a
% different time point. Each row in this data should be a different object
% (i.e. a cell) and the columns should be x-coordinate, y-coordinate,
% object area, tracking index, fluorescence intensities (could be multiple
% columns). The tracking index can be initialized to -1 in every row. It will
% be filled in by MatchFrames so that its value gives the row where the
% data on the same cell can be found in the next frame. 
%2. a frame number (frame). The function will fill in the 4th column of the
% array in peaks{frame-1} with the row number of the corresponding cell in
% peaks{frame} as described above.
%3. A single parameter for the matching (L). In the current implementation of the algorithm, 
% the meaning of this parameter is that objects further than L pixels apart will never be matched. 

% Continue working with the nfkb movie you worked with in hw4. 

% Part 1. Use the first 2 frames of the movie. Segment them any way you
% like and fill the peaks cell array as described above so that each of the two cells 
% has 6 column matrix with x,y,area,-1,chan1 intensity, chan 2 intensity
file1 = 'nfkb_movie1.tif' 
reader = bfGetReader(file1); 
chan_1 = 1;
chan_2 = 2;
time_1 = 1;
time_2 = 2; 
zplane = 6; 

iplane_1 = reader.getIndex(zplane-1, chan_1-1,time_1-1) +1;
iplane_2 = reader.getIndex(zplane-1, chan_1-1,time_2-1) +1;
imgT1 = bfGetPlane(reader, iplane_1);
imgT2 = bfGetPlane (reader, iplane_2);

iplane_3 = reader.getIndex(zplane-1, chan_2-1,time_1-1) +1;
iplane_4 = reader.getIndex(zplane-1, chan_2-1,time_2-1) +1;
imgC2T1 = bfGetPlane(reader, iplane_3);
imgC2T2 = bfGetPlane (reader, iplane_4);

masks = readIlastikFile('nfkb_movie1_Simple Segmentation.h5'); % I got the readIlastikFile function from Github and used Ilastik to segment.
mask_t1 = masks(:,:,1);
mask_t2 = masks(:,:,2);

stats = regionprops(mask_t1, 'Area'); 
figure; hist([stats.Area], 40);
xlabel ('Cell Area', 'FontSize', 24);
ylabel ('Frequency', 'FontSize', 24);
% min area is around 500
minarea= 500;
mask_t1 = imfill(mask_t1, 'holes');
mask_t1 = bwareaopen(mask_t1, minarea);

mask_t2 = imfill(mask_t2, 'holes');
mask_t2 = bwareaopen(mask_t2, minarea);

%making peaks
stats_c1t1 = regionprops(mask_t1, imgT1, 'Centroid', 'Area', 'MeanIntensity');
stats_c1t2 = regionprops(mask_t2, imgT2, 'Centroid', 'Area', 'MeanIntensity');

stats_c2t1 = regionprops(mask_t1, imgC2T1, 'Centroid', 'Area', 'MeanIntensity');
stats_c2t2 = regionprops(mask_t2, imgC2T2, 'Centroid', 'Area', 'MeanIntensity');

% for time point 1
c1t1_xy = cat(1, stats_c1t1.Centroid);
%c1t1_x = c1t1_xy1(:,1);
%c1t1_y = c1t1_xy1(:,2);
c1t1_a = cat(1,stats_c1t1.Area);
c1t1_mi = cat(1, stats_c1t1.MeanIntensity);
c2t1_mi = cat(1, stats_c2t1.MeanIntensity);
tmp_1 = -1*ones(size(c1t1_a)); 

peaks{1} = [c1t1_xy, c1t1_a, tmp_1, c1t1_mi, c2t1_mi]


% for time point 2
c1t2_xy = cat(1, stats_c1t2.Centroid);
%c1t2_x = c1t2_xy(:,1);
%c1t2_y = c1t2_xy(:,2);
c1t2_a = cat(1,stats_c1t2.Area);
c1t2_mi = cat(1, stats_c1t2.MeanIntensity);
c2t2_mi = cat(1, stats_c2t2.MeanIntensity);
tmp_2 = -1*ones(size(c1t2_a)); 

peaks{2} = [c1t2_xy, c1t2_a, tmp_2, c1t2_mi, c2t2_mi]

% Part 2. Run match frames on this peaks array. ensure that it has filled
% the entries in peaks as described above. 
peaks_matched = MatchFrames(peaks, 2, 50);
peaks_matched{1}


% Part 3. Display the image from the second frame. For each cell that was
% matched, plot its position in frame 2 with a blue square, its position in
% frame 1 with a red star, and connect these two with a green line. 

figure; imshow(imgT2, [200,1000]); hold on;
for ii = 1:size(peaks{1})
    plot(peaks{1}(ii,1), peaks{1}(ii,2),'r*','MarkerSize' , 24)  
    
   
   plot(peaks{2}(ii,1), peaks{2}(ii,2),'cs','MarkerSize' , 24)
   nextind = peaks_matched{1}(ii,4);
   if nextind > 0 
   plot([peaks{2}(nextind,1) peaks{1}(ii,1)], [peaks{2}(nextind,2) peaks{1}(ii,2)], 'g');
   end 
end 


