%Nucleus segmentation
%based on convex and concave analysis
function [segmented SplitImage bNoNeedSplit]=NewSegmentBasedConvexConcave(BinaryImage,SplitImage,threshold)

%DistanceThr=27;%threshold for CMU dna images, tunable for different images
% DistanceThr=50;%for images from licolas
DistanceThr=threshold;%for the largest component
%BW=imfill(BinaryImage,'holes'); % filling holes may cause filling
%background.
BW=BinaryImage;
BW=padarray(BW,[1 1],0);  %pad image with 0 in order to avoid out of indices
SplitImage=padarray(SplitImage,[1 1],0);

r=size(BW,1);

if ~islogical(BW)
  BW=im2bw(BW,graythresh(BW));
end
L=bwlabel(BW);

[Boundary bwLabel NumberObject Neighborhood]=bwboundaries(BW);
ConvexHull=regionprops(L,'ConvexHull');%note that the coordinates of the convex hull points are fractional number.

%for debug
% figure;imshow(BW);
% f=JQi_Set_Property_Without_Border(size(BW));imshow(BW,'border','tight');

%Convert the point coordinates from fractional number to integer number
s=size(BW);
TempImage=zeros(s+2);%avoid out of range subscripts
for i=1:length(ConvexHull)
    CHi=ConvexHull(i).ConvexHull;
%     hold on;plot(CHi(:,1),CHi(:,2),'-g','linewidth',2);
    
    column=CHi(:,1)+1;
    row=CHi(:,2)+1;
    fc=floor(column);fr=floor(row);cc=ceil(column);cr=ceil(row);
    TempImage(sub2ind(s+2,fr,fc))=1;TempImage(sub2ind(s+2,fr,cc))=1;
    TempImage(sub2ind(s+2,cr,fc))=1;TempImage(sub2ind(s+2,cr,cc))=1;
end
TempImage=TempImage(2:(s(1)+1),2:(s(2)+1));
ConvexHullIntegerPoints=BW&TempImage;% the integer coordinates of convex hull points
%for dubug
% [x y]=find(ConvexHullIntegerPoints>0);
% plot(y,x,'ro','linewidth',2);
%bNoComponent=false;% indicate no component 
bNoNeedSplit=true;;%bool indicating no need to split
for i=1:length(Boundary)
    if sum(Neighborhood(i,:))% excluding holes
        continue;
    end
    bi=Boundary{i};
    rbi=size(bi,1); %number of row
    count=0;
    maxDistance=0;
    %record indices for each concavest point
    Concavest_Point_Ind=[];
    for j=1:rbi
        if (ConvexHullIntegerPoints(bi(j,1),bi(j,2))==1)&(count==0);
            firstind=j;
            count=count+1;
%             hold on;plot(bi(j,2),bi(j,1),'r+');
            continue;
        end
        if (ConvexHullIntegerPoints(bi(j,1),bi(j,2))==1)&(count==1);
            secondind=j;
            count=count+1;
%             hold on;plot(bi(j,2),bi(j,1),'g+');
        end
        if count==2
            pointA=[bi(firstind,1) bi(firstind,2)];
            pointB=[bi(secondind,1) bi(secondind,2)];
            points=bi(firstind+1:secondind-1,:); %boundary points between Point A and B
            points=padarray(points,[0 1],'post');
            pointA=[pointA 0];pointB=[pointB 0];
            distance=sqrt(sum((cross(points-repmat(pointB,size(points,1),1),repmat(pointA-pointB,size(points,1),1))).^2,2))./norm(pointA-pointB);
            [maxValue maxInd]=max(distance);
            Concavest_Point_Ind=[Concavest_Point_Ind maxInd+firstind];%recording indices for each concavest point
            if maxValue>maxDistance
                maxDistance=maxValue;maxDistanceIndex=maxInd+firstind;
                maxFirstInd=firstind;maxSecondInd=secondind;
%                 hold on;plot(bi(maxDistanceIndex,2),bi(maxDistanceIndex,1),'b+');
            end
            count=1;
            firstind=secondind;
        end
            
    end
    % bi(maxFirstInd:maxSecondInd,:) are boundary points which correspond
    % to the same convex hull straight line as the deepest point
    %bi(1:maxFirstInd,:) and bi(maxSecondInd:end,:) are other distant
    %points which correspond to other convex hull straight lines
    if maxDistance>DistanceThr
        maxDistancePoint=bi(maxDistanceIndex,:);
%         points=bi(1:maxFirstInd,:);
        
%         distance=sqrt(sum((points-repmat(maxDistancePoint,size(points,1),1)).^2,2));
%         [minDistance1 minInd1]=min(distance);
%         
%         points=bi(maxSecondInd:end,:);
%         distance=sqrt(sum((points-repmat(maxDistancePoint,size(points,1),1)).^2,2));
%         [minDistance2 minInd2]=min(distance);
%         
%         if(minDistance1>=minDistance2)
%             minDistanceIndex=minInd2+maxSecondInd;
%         else
%             minDistanceIndex=minInd1;
%         end
        HoleInd=find(Neighborhood(:,i));%find out holes enclosed by boundary i
        points=cell2mat(Boundary(HoleInd));
        bi_Without_Boundary_Curve_SCP_Stand=bi; %compute splitting path between steepest point and boundary points
        bi_Without_Boundary_Curve_SCP_Stand(maxFirstInd:maxSecondInd,:)=[];
        points=[points;bi_Without_Boundary_Curve_SCP_Stand];
%         Concavest_Point_Ind(Concavest_Point_Ind==maxDistanceIndex)=[];%compute splitting path between steppest point and concavest points
%         points=[points;bi(Concavest_Point_Ind,:)];
        
%         points=cell2mat(Boundary);
%         points(maxFirstInd:maxSecondInd,:)=[];
        distance=sqrt(sum((points-repmat(maxDistancePoint,size(points,1),1)).^2,2));
        [minDistance2 minInd2]=min(distance);
        minDistanceIndex=minInd2;
%         if(maxDistanceIndex>=1&&
%         maxDistanceIndex<=rbi&&minDistanceIndex>=1&& minDistanceIndex<=rbi)
        if(maxDistanceIndex>=1&& maxDistanceIndex<=rbi&&minDistanceIndex>=1&& minDistanceIndex<=size(points,1))
%             [x y]=bresenham(bi(maxDistanceIndex,2),bi(maxDistanceIndex,1),bi(minDistanceIndex,2),bi(minDistanceIndex,1));
            [x y]=bresenham(bi(maxDistanceIndex,2),bi(maxDistanceIndex,1),points(minDistanceIndex,2),points(minDistanceIndex,1));
            TempBW=zeros(size(BW));
            TempBW(sub2ind(size(BW),y,x))=1;
            TempBW=imdilate(TempBW,[0 1 0;1 1 1;0 1 0]);
            [x y]=find(TempBW==1);
            %%
%         
%             hold on;plot(y,x,'r-');
%             hold on;plot([y(1) y(end)],[x(1) x(end)],'r-','linewidth',3);
          %for debugging
          %%for debugging
            Neighbor_Offset=[-r-1 -1 r-1 -r 0 r -r+1 1 r+1];%3*3 neighbors including center.
            P1=[bi(maxDistanceIndex,1) bi(maxDistanceIndex,2)];P2=[points(minDistanceIndex,1) points(minDistanceIndex,2)];
            Ind1=sub2ind(size(BW),P1(1),P1(2));Ind2=sub2ind(size(BW),P2(1),P2(2));
            NeighborInd1=Ind1+Neighbor_Offset;NeighborValue1=BW(NeighborInd1)==0;% check background
            NeighborInd2=Ind2+Neighbor_Offset;NeighborValue2=BW(NeighborInd2)==0;
            BlackNeighborInd1=NeighborInd1(NeighborValue1);
            BlackNeighborInd2=NeighborInd2(NeighborValue2);
            if ~isempty(BlackNeighborInd1)
                [P1_r P1_c]=ind2sub(size(BW),BlackNeighborInd1(1));
            end
            if ~isempty(BlackNeighborInd2)
                [P2_r P2_c]=ind2sub(size(BW),BlackNeighborInd2(1));
            end
%             hold on;plot([P1_c P2_c],[P1_r P2_r],'-r','linewidth',3);
            %%
%             hold on;plot([bi(maxDistanceIndex,2) bi(minDistanceIndex,2)],[bi(maxDistanceIndex,1) bi(minDistanceIndex,1)],'r-','linewidth',3);
           % hold on;plot([y(1)+2 y(end)],[x(1)+1 x(end)+2],'r-','linewidth',3);
            BW(sub2ind(size(BW),x,y))=0;
            SplitImage(sub2ind(size(SplitImage),x,y))=0;
            bNoNeedSplit=false;%indicate next iteration to split
       
        end
    else % convex component,no need to split
        BW(bwLabel==i)=0;
    end
end
segmented=BW(2:end-1,2:end-1);%pad operation performed at the entry of this function.
SplitImage=SplitImage(2:end-1,2:end-1);
bNoNeedSplit=bNoNeedSplit;% return value;