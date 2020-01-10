function [centroids_nuclei]=track_nuclei(bwnuclei,zplanes)
%joins nuclei that overlap in z into one object


fprintf('Connecting Objects in Z\n')


IdentityIndex=[];
    for i = 1:zplanes-1 %for each slice

        %   S2 is the labeled object of each slice
        CurrSlice = bwnuclei(:,:,i);
        NxtSlice = bwnuclei(:,:,i+1);
        CurrSliceObj = regionprops(CurrSlice, 'PixelList', 'Centroid');      %objects in current slice
        NxtSliceObj = regionprops(NxtSlice, 'PixelList', 'Centroid');        %objects in next slice

        for j=1:length(CurrSliceObj)  %for each object in CurrSlice


            a = CurrSliceObj(j).PixelList;      %look at pixel list of the object
            N_old = 0;                          %restarts number of overlap = 0
            IdentityIndex(j,i)=j;               %index for each object   
            for k=1:length(NxtSliceObj)         %compare of each object in NxtSlice

                b = NxtSliceObj(k).PixelList;   %look at the pixel list
                c = intersect(a,b,'rows');
                N = length(c(:,1));         %number of overlapping pixels
                if N > N_old
                    N_old = N;
                    Identity(j,i) = k; %remember the greatest overlappping kth object in Z+1 corresponding to the jth object in current Z
                else
                    if N_old==0
                    Identity(j,i)=0;
                    end
                end

            end

        end

    end
    
        
    IdentityTransformed = [];
    IdentityTransformed(:,1)=IdentityIndex(:,1);

    for j = 2:length(Identity(1,:)) %each column of Identity (34-1 columns)

       for i= 1:length(IdentityTransformed(:,1))   %each row
           if IdentityTransformed(i,j-1)~=0

            IdentityTransformed(i,j)=Identity(IdentityTransformed(i,j-1),j-1);
           end
       end
        A=IdentityIndex(:,j);
        B=Identity(:,j-1);
        C=setdiff(A,B); %find objects that doesn't have a overlap partner
        lastindex=find(IdentityTransformed(:,j)~=0,1,'last');
        lengthC=length(C);
        IdentityTransformed(lastindex+1:lastindex+lengthC,j)=C(:);

    end

    ObjNum = length(IdentityTransformed(:,1));
    Zslices = length(IdentityTransformed(1,:));
    color_code = rand(3,ObjNum);

    % Initialize Empty Matrix

%     NewIdentity_color = zeros(size(bwnuclei,1),size(bwnuclei,2),size(bwnuclei,3));
    NewIdentity_Index = zeros(size(bwnuclei,1),size(bwnuclei,2),size(bwnuclei,3));

    for j = 1:Zslices
        S2_current_slice = bwnuclei(:,:,j);
        S2_region_props = regionprops(S2_current_slice, 'PixelList');

        for i = 1:ObjNum


            if IdentityTransformed(i,j) ~=0

            CurrentCell = IdentityTransformed(i,j);


            pixel_coord_X = S2_region_props(CurrentCell).PixelList(:,1);
            pixel_coord_Y = S2_region_props(CurrentCell).PixelList(:,2);

                for k = 1:length(pixel_coord_X)

                    %   Assign Color to each Index
%                     NewIdentity_color(pixel_coord_X(k),pixel_coord_Y(k),1,j)= color_code(1,i);
%                     NewIdentity_color(pixel_coord_X(k),pixel_coord_Y(k),2,j)= color_code(2,i);
%                     NewIdentity_color(pixel_coord_X(k),pixel_coord_Y(k),3,j)= color_code(3,i);

                    %   Assigning new Index
                    NewIdentity_Index(pixel_coord_X(k),pixel_coord_Y(k),j)= i;
                end

            end

        end
    end

    %Regionprop NewIdentity_Index Matrix to get object info
    NewPixel_Props= regionprops(NewIdentity_Index, 'PixelList', 'Centroid');
    
    %filter out super small objects
for i=1:length(NewPixel_Props)
    good(i)=length(NewPixel_Props(i).PixelList)>50;
end
NewPixel_Props=NewPixel_Props(good);
    
    for i=1:length(NewPixel_Props)
        centroids_nuclei(i,:)=NewPixel_Props(i).Centroid;
    end
[displayimg]=make_colored_display(NewPixel_Props,zplanes);    
    
function [displayimg]=make_colored_display(object_properties,zplanes)

    displayimg=zeros(1024,1024,3,zplanes);
    for i=1:length(object_properties)
        color=rand(1,3);
        pix=object_properties(i).PixelList;
        for j=1:length(pix)
            loc=pix(j,:);
            displayimg(loc(1),loc(2),:,loc(3))=color;
        end
    end
    implay(displayimg)