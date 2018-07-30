function grayscaleimages = getTIFFs2gray(gpuAvailable)

ImageName = uigetfile('.tiff','Pick images',[],'MultiSelect', 'on');  %choosing the images for the registration

numSamples = size(ImageName,2);


if ischar(ImageName) == 1
        disp(ImageName);
        info = imfinfo(ImageName);
        num_images = numel(info);
        
        for i = 1:num_images; % converting the tiff stack into an grayscale image for elastix to work
            A = imread(ImageName,i,'Info', info);
            B = (A(:,:,3));
            C = (A(:,:,2));
            D = (A(:,:,1));
                blueImage (:,:,i) = B;
                greenImage(:,:,i) = C;
                redImage  (:,:,i) = D;
        end
            grayscaleimages{1} = blueImage;
            grayscaleimages{2} = greenImage;
            grayscaleimages{3} = redImage;
            clear blueImage greenImage redImage
        
        
else
    for j = 1:numSamples
        disp(ImageName{j});
        info = imfinfo(ImageName{j});
        num_images = numel(info);
        
        
        for i = 1:num_images % converting the tiff stack into an grayscale image for elastix to work
            A = imread(ImageName{1,j},i,'Info', info);
            B = (A(:,:,3));
            C = (A(:,:,2));
            D = (A(:,:,1));
%             if gpuAvailable == 1
%                 grayImage(:,:,i) = gpuArray(B);
%             else grayImage(:,:,i) = B;
                blueImage (:,:,i) = B;
                greenImage(:,:,i) = C;
                redImage  (:,:,i) = D;
                
        end
            
            grayscaleimages{1,j} = blueImage;
            grayscaleimages{2,j} = greenImage;
            grayscaleimages{3,j} = redImage;
            clear blueImage greenImage redImage
%         end
    end
end

end

