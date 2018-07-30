function grayscaleimage = getTIFF2gray(gpuAvailable)

ImageName = uigetfile('.tiff','Pick reference image',[]);  %choosing the fixed image for the registration
disp(ImageName);
info = imfinfo(ImageName);
num_images = numel(info);


for i = 1:num_images; % converting the tiff stack into an grayscale image for elastix to work
    A = imread(ImageName,i,'Info', info);
    B = (A(:,:,3));
%     if gpuAvailable == 1;
%         grayImage(:,:,i) = gpuArray(B);
%     else grayImage(:,:,i) = B;
        grayImage(:,:,i) = B;
end
grayscaleimage{1} = grayImage;


for i = 1:num_images; % converting the tiff stack into an grayscale image for elastix to work
    A = imread(ImageName,i,'Info', info);
    B = (A(:,:,2));
%     if gpuAvailable == 1;
%         grayImage(:,:,i) = gpuArray(B);
%     else grayImage(:,:,i) = B;
        grayImage(:,:,i) = B;
end
grayscaleimage{2} = grayImage;


for i = 1:num_images; % converting the tiff stack into an grayscale image for elastix to work
    A = imread(ImageName,i,'Info', info);
    B = (A(:,:,1));
%     if gpuAvailable == 1;
%         grayImage(:,:,i) = gpuArray(B);
%     else grayImage(:,:,i) = B;
        grayImage(:,:,i) = B;
end


grayscaleimage{3} = grayImage;

end
