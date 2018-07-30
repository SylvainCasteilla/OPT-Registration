  clear Images vector V ifolderOutput resultImages resultImage;
% 
try
    d = gpuDevice;
    gpuAvailable = d.SupportsDouble;
    disp('GPU available!');
catch
    gpuAvailable = false;
    disp('No Nvidia GPU detected.');
end

%% choosing tiff stacks for the registration,
%%% opening a file chooser for the reference image, and then one for the other ones
%%% tiff stacks' channels are separated and converted to grayscale 3D matrices
%% 
      fixImage = GetTIFF2Gray(gpuAvailable); %choosing the fixed image
      movImage = GetTIFFs2Gray(gpuAvailable); %choosing the moving images in a cell array
% % 
tic
%%% Images will be the cell containing the channel on which the registration is done
%%% Images{i,1} is the blue channel used for alignment and registrations,
%%% Images{i,2} is the green channel and Images{i,3} is the red channel
%% 
if size(movImage,1) > 1
    numSamples = size(movImage,2);
    images = cell(numSamples+1,3); 
    for j = 1:3
        images{1,j} = fixImage{j};
    end
    for i = 1:numSamples
        for j = 1:3
            images{i+1,j} = movImage{j,i};
        end
    end

elseif size(movImage,1) == 1
    numSamples = 1;
    images = cell(1,2);
    for j = 1:3
        images{1,j} = fixImage{j};
    end
    for j = 1:3
        images{2,j} = movImage{j};
    end
end
clear fixImage movImage ;
numImages = size(images,1);

%% selecting the input folder, which contain the parameters files for the 2 steps of the regiatration
folderInput  = uigetdir ([] ,'Select the input folder containing parameters files');
disp(folderInput);

%% choosing the parameters file for elastix
paramFile = uigetfile('.txt','select the first registration parameter file',folderInput);
paramFile = {strcat(folderInput,'\',paramFile)};

%% choosing the parameters file for the second elastix registration
paramFile_2 = uigetfile('.txt','select the second registration parameter file',folderInput);
paramFile_2 = {strcat(folderInput,'\',paramFile_2)};


%% selecting the output folder for computation
%%% the different output folders for images registrations will be written
%%% in this folder
folderOutput = uigetdir ([] ,'Select the output folder for your registration'); 
disp(folderOutput);




%% finding maximum length and size of all the images
imagesHeight = zeros(size(numImages,1));
imagesLength = zeros(size(numImages,1));
for i = 1 : numImages
    imagesHeight(i) = size(images{i,1},1);
    imagesLength(i) = size(images{i,1},2);
end

resultImage = zeros(numImages,2);

Mx = max(imagesHeight);
My = max(imagesLength);
for i = 1 : numImages
    resultImage(i,:) = [Mx - imagesHeight(i),My - imagesLength(i)] ; 
end

%%  pretreatment of the images
h = waitbar(0,'cleaning images...');
direction = zeros(size(images));
for i = 1 : numImages
    %%% the images are thresold to eliminate the pixels of low intensity,
    %%% mostly in the background of the images
    for j = 1:3
        images{i,j}(images{i,j}<9) = 0;
    end
    
    %%% the corners of the images is also cut, to cut some residuals
    %%% capillaries pixels
    triangles = zeros(size(images{i,1}));

        for j = 0:180
            for k = 0:80
                 if (4*k+j) <= 180
                    triangles(k+1,j+1,:) = 1;
                    triangles(imagesHeight(i)-k,j+1,:) = 1;
                    triangles(k+1,imagesLength(i)-j,:) = 1;
                    triangles(imagesHeight(i)-k,imagesLength(i)-j,:) = 1;
                 end
                 if (k+2*j) <= 80 
                    triangles(k+1,j+1,:) = 1;
                    triangles(imagesHeight(i)-k,j+1,:) = 1;
                    triangles(k+1,imagesLength(i)-j,:) = 1;
                    triangles(imagesHeight(i)-k,imagesLength(i)-j,:) = 1;
                 end
            end
        end
        for j = 1:3
            images{i,j}(triangles == 1) = 0;
        end
        waitbar(i/numImages)
end
clear triangles
close(h)
%%% images are all resized to the larger dimensions of the images
%%% the direction (up or down) they're facing is determined
h = waitbar(0,'resizing images...');

for i = 1:numImages    
    for j = 1:3
        images{i,j} = padarray(images{i,j},[resultImage(i,1),resultImage(i,2),0],0,'post');
        direction(i) = head_direction(images{i,1});
    


%%% images are tranlated to the center of the expanded images
%%% if the direction of an image is not down, it is rotated to be in the
%%% same direction as the others
        images{i,j} = imtranslate(images{i,j}, [resultImage(i,2)/2,resultImage(i,1)/2,0]);
        
        if direction(i) == 1
            images{i,j} = imrotate3(images{i,j},180,[1,0,0],'linear','crop','FillValues',0);
            images{i,j} = imrotate3(images{i,j},180,[0,0,1],'linear','crop','FillValues',0);
            images{i,j} = imtranslate(images{i,j}, [0,0,-2]);
        end
    end
    waitbar(i/numImages)

end
close(h)
%% saving the pretreated images to a folder called Images
h = waitbar(0,'saving resized images...');

imagesFolder = fullfile(folderOutput,'Images',filesep);
mkdir(imagesFolder);
for i = 1:numImages     
    ImageName = sprintf('image_%d',i);
    ImageName = strcat(imagesFolder,ImageName,'.tif');
    imwrite(images{i,1}(:,:,1),ImageName);
     for j = 2:1024 
        imwrite(images{i,1}(:,:,j),ImageName,'WriteMode','append'); 
     end
    waitbar(i/numImages)

end
close(h)

toc

%%  repositioning the moving images to match the fixed image

%%% initialization of the differents angles and vector used
angle1 = [];
angle2 = [];
angle3 = [];
vector = zeros(numSamples,3);


h = waitbar(0,'aligning images...');

for i = 1:numSamples
    tempImage = images{i+1,1};
    
%%% InitialAlignment determine the center of gravity of the two compared
%%% images and align them together, giving vector, and then rotate the
%%% images along the 3 axis to find the angles of rotation
        [Vx,Vy,Vz,angleX,angleY,initialAngle] = InitialAlignment(images{1},tempImage,gpuAvailable);

    angle1(i) = angleX;
    angle2(i) = angleY;
    angle3(i) = initialAngle;
    vector(i,:) = [Vy, Vx, Vz];
    for j = 1:3
        images{i+1,j} = imtranslate(images{i+1,j}, vector(i,:));
        images{i+1,j} = imrotate3(images{i+1,j},angleX,[1,0,0],'linear','crop','FillValues',0);
        images{i+1,j} = imrotate3(images{i+1,j},angleY,[0,1,0],'linear','crop','FillValues',0);
        images{i+1,j} = imrotate3(images{i+1,j},initialAngle,[0,0,1],'linear','crop','FillValues',0);
    end
    waitbar(i/numSamples)

end
  clear tempImage;
close(h)

%% saving the pretreated images to a folder called aligned_Images
h = waitbar(0,'saving aligned images...');

alignedImagesFolder = fullfile(folderOutput,'Aligned_Images',filesep);
mkdir(alignedImagesFolder);
for i = 1:numImages
    ImageName = sprintf('image_TR_%d',i);
    ImageName = strcat(alignedImagesFolder,ImageName,'.tif');
     imwrite(images{i,1}(:,:,1),ImageName);
     for j = 2:1024 
        imwrite(images{i,1}(:,:,j),ImageName,'WriteMode','append'); 
     end
         waitbar(i/numImages)

end
close(h)

%% computing elastix on rough parameters to obtain new aligned image and reference
h = waitbar(0,'computing elastix...');

for i = 1:numSamples
    
    %%% we create an output folder for each image (avoid elastix to overwrite the result image)
    ifolderOutput_1{i} = fullfile(folderOutput,'Registration_1', sprintf('image%d',i));
    
    %%% computing elastix on a rigid transform
    elastix(images{i+1,1},images{1,1},ifolderOutput_1{i},paramFile);
    waitbar(i/numSamples)
end
close(h)
toc

%% obtaining the result images for the second registration
reg1 = fullfile(folderOutput,'rigid_registration',filesep);
mkdir(reg1);
resultImage = uint8(zeros(size(images{1,1})));
for i = 1:numSamples
    
    resultName = strcat(ifolderOutput_1{1,i},'\result.0.mhd');
    resultInfo = mhd_read_header(resultName);
    disp(resultName);
    num_images_result = resultInfo.Dimensions(3);

    
    %%% converting the result mhd header into a grayscale image of the
    %%% right format so elastix works
    resultImage = mhd_read(resultInfo);
    resultImage = uint8(resultImage);
    
    %%% converting the result tiff stack into a grayscale image of the
    %%% right format so elastix works
    
    
    resultImages{i} = resultImage;
    clear resultImage;
    
end

%% saving the result images in a folder called reg1
h = waitbar(0,'saving first registration...');

for i = 1:numSamples
    ImageName = sprintf('image_result_%d',i);
    ImageName = strcat(reg1,ImageName,'.tif');
     imwrite(resultImages{i}(:,:,1),ImageName);
     for j = 2:1024 
        imwrite(resultImages{i}(:,:,j),ImageName,'WriteMode','append'); 
     end
     waitbar(i/numSamples)
end
close(h)
%% creating an average of the first aligned images 
averageAlignedImage = uint8(zeros(size(resultImages{i})));
for i =  1:numSamples
    averageAlignedImage = resultImages{i}+averageAlignedImage ;
end


%averageAlignedImage = (averageAlignedImage/numSamples);
averageImage = fullfile(folderOutput,'averageImage.tiff');
imwrite(averageAlignedImage (:,:,1),averageImage);
for j = 2:num_images_result
    imwrite(averageAlignedImage (:,:,j),averageImage,'WriteMode','append');
end


%% second elastix registration
h = waitbar(0,'computing elastix...');

for i = 1:numSamples
    
    %%% we create an output folder for each image (avoid elastix to overwrite the result image)
    ifolderOutput_2{i} = fullfile(folderOutput,'Registration_2', sprintf('image%d',i));
    
    %%% computing elastix
    elastix(resultImages{i},images{1,1},ifolderOutput_2{i},paramFile_2);
    
    waitbar(i/numSamples)
end
close(h)

averageRegisteredImage = uint8(zeros(size(resultImages{i})));

%% saving the result images in a folder called reg2
reg2 = fullfile(folderOutput,'affine_registration',filesep);
mkdir(reg2);
h = waitbar(0,'saving second registration...');

for i = 1:numSamples
    
    resultName = strcat(ifolderOutput_2{1,i},'\result.0.mhd');
    resultInfo = mhd_read_header(resultName);
    disp(resultName);
    num_images_result = resultInfo.Dimensions(3);

    
    %%% converting the result mhd header into a grayscale image of the
    %%% right format so elastix works
    resultImage = mhd_read(resultInfo);
    resultImage = uint8(resultImage);
    ImageName = sprintf('image_result_%d_reg2',i);
    ImageName = strcat(reg2,ImageName,'.tif');
    imwrite(resultImage(:,:,1),ImageName);
    for j = 2:num_images_result 
        imwrite(resultImage(:,:,j),ImageName,'WriteMode','append');
    end 
    
    resultImages{i} = resultImage;
    
    
    averageRegisteredImage = resultImage+averageRegisteredImage;
    clear resultImage;
    
    waitbar(i/numSamples)
end
close(h)


%% saving an average of the registered images
regImage = fullfile(folderOutput,'FinalImage.tiff');

imwrite(averageRegisteredImage (:,:,1),regImage);
for j = 2:num_images_result 
    imwrite(averageRegisteredImage (:,:,j),regImage,'WriteMode','append');
end


%% using transformix to obtain the final images on the three channels
h = waitbar(0,'applying transformations to all channels...');

for i = 1:numSamples
    for j = 1:3
        %%% transformParametersName is the folder where the
        %%% transformparameters file is located (elastix output)
        transformParametersName = strcat(ifolderOutput_1{1,i});
        
        %%% transformix apply the transformation to the image given and
        %%% produce an image 
        images{i+1,j} = uint8(transformix(images{i+1,j},transformParametersName));
        
        %%% transformParametersName is the folder where the
        %%% transformparameters file is located (elastix output)
        transformParametersName = strcat(ifolderOutput_2{1,i});
        
        %%% transformix apply the transformation to the image given and
        %%% produce an image 
        images{i+1,j} = uint8(transformix(images{i+1,j},transformParametersName));
    end
    waitbar(i/numSamples)
end
close(h)


%% saving the colored images in the folder Tif_Images
tifImages = fullfile(folderOutput,'Tif_Images',filesep);
colorTiff = zeros(Mx, My, 3, 1024, 'uint8');
mkdir(tifImages);
h = waitbar(0,'saving color result images...');


for i = 1:numImages
    ImageName = sprintf('image_result_%d',i);
    ImageName = strcat(tifImages,ImageName,'.tiff');
    for j = 1:3
        colorTiff(:,:,j,:) = images{i,4-j}(:,:,:);
    end
    options.overwrite = true;
        options.color = true;
        options.message = false;
        saveastiff(colorTiff, ImageName, options);
    waitbar(i/numImages)
end
close(h)
toc

% %     ImageName = sprintf('image_TR_2bis',i);
% %     ImageName = strcat(ImageName,'.tif');
% %      imwrite(Images{2,2}(:,:,1),ImageName);
% %      for j = 2:1024 
% %         imwrite(Images{2,2}(:,:,j),ImageName,'WriteMode','append'); 
% %      end