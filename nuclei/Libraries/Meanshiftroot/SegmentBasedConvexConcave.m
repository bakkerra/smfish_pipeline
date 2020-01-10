%Nucleus segmentation
%based on convex and concave analysis
function segmented=SegmentBasedConvexConcave(BinaryImage,threshold )

DistanceThr=threshold;%threshold default 3

BW=imfill(BinaryImage,'holes');
L=bwlabel(BW);

Boundary=bwboundaries(BW);
ConvexHull=regionprops(L,'ConvexHull');%note that the coordinates of the convex hull points are fractional number.

%for debug
% figure;imshow(BW);

%Convert the point coordinates from fractional number to integer number
s=size(BW);
TempImage=zeros(s+2);%avoid out of range subscripts
for i=1:length(ConvexHull)
    CHi=ConvexHull(i).ConvexHull;
%     hold on;plot(CHi(:,1),CHi(:,2));
    
    column=CHi(:,1)+1;
    row=CHi(:,2)+1;
    fc=floor(column);fr=floor(row);cc=ceil(column);cr=ceil(row);
    TempImage(sub2ind(s+2,fr,fc))=1;TempImage(sub2ind(s+2,fr,cc))=1;
    TempImage(sub2ind(s+2,cr,fc))=1;TempImage(sub2ind(s+2,cr,cc))=1;
end
TempImage=TempImage(2:(s(1)+1),2:(s(2)+1));
ConvexHullIntegerPoints=BW&TempImage;% the integer coordinates of convex hull points
%for dubug
%[x y]=find(ConvexHullIntegerPoints>0);
%plot(y,x,'r*');

for i=1:length(Boundary)
    bi=Boundary{i};
    rbi=size(bi,1); %number of row
    count=0;
    maxDistance=0;
    
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
            points=bi(firstind+1:secondind-1,:);
            points=padarray(points,[0 1],'post');
            pointA=[pointA 0];pointB=[pointB 0];
            distance=sqrt(sum((cross(points-repmat(pointB,size(points,1),1),repmat(pointA-pointB,size(points,1),1))).^2,2))./norm(pointA-pointB);
            [maxValue maxInd]=max(distance);
            if maxValue>maxDistance
                maxDistance=maxValue;maxDistanceIndex=maxInd+firstind;
                maxFirstInd=firstind;maxSecondInd=secondind;
%                 hold on;plot(bi(maxDistanceIndex,2),bi(maxDistanceIndex,1),'b+');
            end
            count=1;
            firstind=secondind;
        end
            
    end
    if maxDistance>DistanceThr
        points=bi(1:maxFirstInd,:);
        maxDistancePoint=bi(maxDistanceIndex,:);
        distance=sqrt(sum((points-repmat(maxDistancePoint,size(points,1),1)).^2,2));
        [minDistance1 minInd1]=min(distance);
        
        points=bi(maxSecondInd:end,:);
        distance=sqrt(sum((points-repmat(maxDistancePoint,size(points,1),1)).^2,2));
        [minDistance2 minInd2]=min(distance);
        
        if(minDistance1>=minDistance2)
            minDistanceIndex=minInd2+maxSecondInd;
        else
            minDistanceIndex=minInd1;
        end
        if(maxDistanceIndex>=1&& maxDistanceIndex<=rbi&&minDistanceIndex>=1&& minDistanceIndex<=rbi)
            [x y]=bresenham(bi(maxDistanceIndex,2),bi(maxDistanceIndex,1),bi(minDistanceIndex,2),bi(minDistanceIndex,1));
            TempBW=zeros(size(BW));
            TempBW(sub2ind(size(BW),y,x))=1;
            TempBW=imdilate(TempBW,[0 1 0;1 1 1;0 1 0]);
            [x y]=find(TempBW==1);
        
%             hold on;plot(y,x,'r-');
            BW(sub2ind(size(BW),x,y))=0;
        end
    end
end
segmented=BW;