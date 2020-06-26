%----Supplementary code accompanying Part2 Chapter1 of PVR thesis----%

% Filtered back-projection simulations with phantom image

% Pedro Vallejo Ramirez
% Created: 18/09/2017

% Updated on 2020/01/14

% This script is used to load a phantom image and test filtered
% back-projection on it. 


% Load the phantom and define the amount of number of projections (256)
proj = 128;
p = phantom(256);
L = length(p);
delta_theta = 360/proj;
theta = 0:delta_theta:359;


%% Radon transform
%Take the radon transform of the image
R = radon(p,theta);
output_size = max(size(p));
iRad = zeros(output_size,output_size,L);
iRad_nofilter = zeros(output_size,output_size,L);


%iRad = iradon(R(:,100),delta_theta,output_size);

% Show Radon transform
%figure;imshow(R,[])
% plot first projection in the radon transform
%figure; plot(R(:,5))

% Truncate the sinogram and reconstruct it one by one
for i = 1:proj
    %truncate sinogram
    Rt = R(:,1:i); 
    iRad(:,:,i) = iradon(Rt,delta_theta,output_size);
    iRad_nofilter(:,:,i) = iradon(Rt,delta_theta,output_size,'linear','none');
    
%     imshow(iRad(:,:,i),[]);
%     pause
end

%% Display a few steps in the FBP process
% for 16 projections
figure
imshow(iRad(:,:,1),[]);
figure
imshow(iRad(:,:,5),[]);
figure
imshow(iRad(:,:,4),[]);
figure
imshow(iRad(:,:,10),[]);
figure
imshow(iRad(:,:,16),[]);

% for 256 projections
figure
imshow(iRad(:,:,1),[]);
figure
imshow(iRad(:,:,20),[]);
figure
imshow(iRad(:,:,128),[]);

figure
imshow(iRad_nofilter(:,:,128),[]);
%%
% Can I show it's equivalent to reconstruct an image from a single
% projection == to taking the FT of the image and multiplying it by a line
% filter at the given angle theta? 

% get fourier transform of original iamge
ft = myfft2(p);
figure;imshow(ft,[])

% make a mask with the same size and only 1s in a line through the
% origin
mask =  zeros(256);
for j = 1:size(mask,2) 
    mask(129:129,j) = 1;
end 
mask90 = rot90(mask);
% multiply the mask by the FT
filtered_phantom = ft.*mask90;

% this works
ft_try2 = fftshift(fft2(p));
ft_filt = ft_try2 .* mask90;
figure;imshow((ifft2(ft_filt)),[]);

% multiply the mask by the FT
filtered_phantom = ft.*mask;
% this works
ft_try2 = fftshift(fft2(p));
ft_filt = ft_try2 .* mask;
figure;imshow((ifft2(ft_filt)));

% The inversion yields a distorted reconstruction of the sample, as
% expected, but it's not necesarily equivalent to the FBP smearing of only
% a single projection. 


% iFT the filtered version to recover the image of the object
object_filtered = ifft2(filtered_phantom);
figure;imshow(abs(object_filtered),[]);

%% Saving the images
%Rad(Rad<0) = 0;
%iRad = uint32(iRad); % Can't save a single with imwrite. 

%[path, name, type] = fileparts(filepath);
%saveTiff(iRad, 'FBPmovie');


