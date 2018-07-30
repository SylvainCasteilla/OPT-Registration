function direction = head_direction (Image)

image_height = size(Image,1);
image_width = size(Image,2);
numframes = size(Image,3);

Intensity = zeros(1,numframes);
P_Intensity = zeros(1,numframes);

if isa(Image,'uint8')
    temp = Image;
else
    temp = gather(Image);
end


for i=1:numframes
    
    I = temp(:,:,i);
    Intensity(i) = sum(sum(sum(I)))/(image_height*image_width);
    P_Intensity(i) = Intensity(i)*i;
    
end

center = ceil(sum(P_Intensity)/sum(Intensity));
head = find(round(abs(diff(Intensity))*1000) == max(round(abs(diff(Intensity))*1000)));

if head - center < 0
    direction = 1;
    %disp('fish is head up');
else
    direction = -1;
    %disp('fish is head down');
    
end


 
     
    
    