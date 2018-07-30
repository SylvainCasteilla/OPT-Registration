function [Vx,Vy,Vz,angle1,angle2,angle3] = InitialPosition(img1,img2,gpuAvailable)
% InitialPosition will compute the first alignment parameters for the
% different images
%%% the different translations are mostly found with the following
%%% algorithm : the coordinates of the center of gravity are determined for
%%% both images then they are translated to superpose the centers
%%% 
%%% for the angles on the different axis, the algorithm compares the number
%%% of non-zeros pixels the two images have in common. maximazing the
%%% number of pixels in common gives the angle of rotation

clear Image1 Image2 
    tic
    image_height_1 = size(img1,1);
    image_width_1 = size(img1,2);
    numframes = size(img1,3);
    image_height_2 = size(img2,1);
    image_width_2 = size(img2,2);
        Image1 = img1;
        Image2 = img2;
        

    Mx = max(image_height_1,image_height_2);
    My = max(image_width_1,image_width_2);
  
  
%% translating the fishes to have the same number of considered frames 
        for i=1:numframes
    
            I_1 = Image1(:,:,i);
            Intensity_1z(i) = sum(sum(sum(I_1)))/(Mx*My);
            P_Intensity_1z(i) = Intensity_1z(i)*i;
    
            I_2 = Image2(:,:,i);
            Intensity_2z(i) = sum(sum(sum(I_2)))/(Mx*My);
            P_Intensity_2z(i) = Intensity_2z(i)*i;
        end
    
    %%% head is the position where the intensity "gradient" is the most
    %%% important, giving the last frame containing fish pixels
    head_1 = find(round(abs(diff(Intensity_1z))*1000) == max(round(abs(diff(Intensity_1z))*1000)),1);
    head_2 = find(round(abs(diff(Intensity_2z))*1000) == max(round(abs(diff(Intensity_2z))*1000)),1);

    Vz = head_1-head_2;
    Image2 = imtranslate(Image2,[0,0,Vz]);
   
    
%% rotation on x and y axis on a 60 degree window   
    for i = -30:30
        nbpixelsInCommon = 0;

        for j = 100:10:My-100
            rotIm = imrotate(squeeze(Image2(:,j,:)),i,'bilinear','crop');
            nbpixelsInCommon = size(find(squeeze(Image1(:,j,:))& rotIm(:,:)),1)+nbpixelsInCommon;
        end
        optimizer1(i+31) = nbpixelsInCommon;
    end
    
    M = max(optimizer1);
    angle1 = find(optimizer1 == M,1)-31;
    

    
    for i = -30:30
        nbpixelsInCommon = 0;

        for j = 100:10:Mx-100
            rotIm = imrotate(squeeze(Image2(j,:,:)),i,'bilinear','crop');
            nbpixelsInCommon = size(find(squeeze(Image1(j,:,:))& rotIm(:,:)),1)+nbpixelsInCommon;
        end
        optimizer2(i+31) = nbpixelsInCommon;
    end
    
    M = max(optimizer2);
    angle2 = find(optimizer2 == M,1)-31;
    
    %%% if the angle is smaller than 5 degrees, there is no need to rotate
    
    if abs(angle1) > 5 
        Image2 = imrotate3(Image2,angle1,[1,0,0],'linear','crop','FillValues',0);
        Image2 = imrotate3(Image2,angle2,[0,1,0],'linear','crop','FillValues',0);
    else
        angle1 = 0;
        angle2 = 0;
    end

    
    
    
 %% findind the coordinates of the center of gravity of the two images, 
 %% and translating the second image to have the same center as the first one 
        for i=1:numframes
    
            I_1 = Image1(:,:,i);
            Intensity_1z(i) = sum(sum(sum(I_1)))/(Mx*My);
            P_Intensity_1z(i) = Intensity_1z(i)*i;
    
            I_2 = Image2(:,:,i);
            Intensity_2z(i) = sum(sum(sum(I_2)))/(Mx*My);
            P_Intensity_2z(i) = Intensity_2z(i)*i;
        end
    
        center_1_z = round(sum(P_Intensity_1z)/sum(Intensity_1z));
        center_2_z = round(sum(P_Intensity_2z)/sum(Intensity_2z));

        for i=1:Mx
    
            I_1 = Image1(i,:,:);
            Intensity_1x(i) = sum(sum(sum(I_1)))/(numframes*My);
            P_Intensity_1x(i) = Intensity_1x(i)*i;
    
            I_2 = Image2(i,:,:);
            Intensity_2x(i) = sum(sum(sum(I_2)))/(Mx*My);
            P_Intensity_2x(i) = Intensity_2x(i)*i;
        end
        center_1_x = round(sum(P_Intensity_1x)/sum(Intensity_1x));
        center_2_x = round(sum(P_Intensity_2x)/sum(Intensity_2x));
        
        for i=1:My
    
            I_1 = Image1(:,i,:);
            Intensity_1y(i) = sum(sum(sum(I_1)))/(Mx*numframes);
            P_Intensity_1y(i) = Intensity_1y(i)*i;
    
            I_2 = Image2(:,i,:);
            Intensity_2y(i) = sum(sum(sum(I_2)))/(Mx*numframes);
            P_Intensity_2y(i) = Intensity_2y(i)*i;
        end
        center_1_y = round(sum(P_Intensity_1y)/sum(Intensity_1y));
        center_2_y = round(sum(P_Intensity_2y)/sum(Intensity_2y));
    
       
        Vx = center_1_x - center_2_x;
        Vy = center_1_y - center_2_y;
        Vz1 = center_1_z - center_2_z;
        
        
     Image2 = imtranslate(Image2,[Vy,Vx,0]);
     Vz = Vz;

%% rotating on the z axis
    for i = 0:359
        nbpixelsInCommon = 0;

        for j = 1:16:head_1
            rotIm = imrotate(squeeze(Image2(:,:,j)),i,'bilinear','crop');
            nbpixelsInCommon = size(find(squeeze(Image1(:,:,j))& rotIm(:,:)),1)+nbpixelsInCommon;
        end


        optimizer(i+1) = nbpixelsInCommon;
    end
    
    
    M = max(optimizer);
    angle3 = find(optimizer == M,1)-1;
    Image2 = imrotate3(Image2,angle3,[0,0,1],'linear','crop','FillValues',0);


%% second translation on x and y axis to be more precise                
        for i=1:Mx   
            I_2 = Image2(i,:,:);
            Intensity_2x(i) = sum(sum(sum(I_2)))/(Mx*My);
            P_Intensity_2x(i) = Intensity_2x(i)*i;
        end
        
        center_2_x = round(sum(P_Intensity_2x)/sum(Intensity_2x));
        
        for i=1:My    
            I_2 = Image2(:,i,:);
            Intensity_2y(i) = sum(sum(sum(I_2)))/(Mx*numframes);
            P_Intensity_2y(i) = Intensity_2y(i)*i;
        end
        
        center_2_y = round(sum(P_Intensity_2y)/sum(Intensity_2y));
    
        Vx1 = center_1_x - center_2_x;
        Vy1 = center_1_y - center_2_y;
        
     Image2 = imtranslate(Image2,[Vy1,Vx1,0]);
     Vx = Vx + Vx1;
     Vy = Vy + Vy1;
     
     
%% second rotation following the z axis     
   
    for i = 0:359
        nbpixelsInCommon = 0;

        for j = 1:16:head_1
            rotIm = imrotate(squeeze(Image2(:,:,j)),i,'bilinear','crop');
            nbpixelsInCommon = size(find(squeeze(Image1(:,:,j))& rotIm(:,:)),1)+nbpixelsInCommon;
        end


        optimizer(i+1) = nbpixelsInCommon;
    end
    
    
    M = max(optimizer);
    angle3 = mod(find(optimizer == M,1)-1 +angle3,360);

    toc
end