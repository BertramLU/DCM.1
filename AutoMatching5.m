% auto matching algorithm

% set the parameters
DOMDataPath='ChutouDOM.tif';
TrajectoryDataPath='trajectory.pcd';
ChannelCenter=[353468.268, 3470283.409];
ChannelGreyColor=210;
ChannelGreyColorMix=20;
ScalingInterval=[0.8, 2.0];
ScalingStep=0.01;
DOMFourierSmoothIterationCount=300;

% read data
% read geotiff-type data
[A0b,R0] = readgeoraster(DOMDataPath);
A0=rgb2gray(A0b);
% read trajectory point cloud
B0=pcread(TrajectoryDataPath);

% trajectory data pre-process
B1=B0.Location(:,1:2);
TraLength=sqrt((B0.XLimits(1,1)-B0.XLimits(1,2)).^2+(B0.YLimits(1,1)-B0.YLimits(1,2)).^2);
AvrTra(1,1)=sum(B1(:,1))./B0.Count;
AvrTra(1,2)=sum(B1(:,2))./B0.Count;
B1(:,3)=B1(:,1)-AvrTra(1,1);
B1(:,4)=B1(:,2)-AvrTra(1,2);
RefineXCount=floor(TraLength.*0.03./R0.CellExtentInWorldX).*2+1;
RefineYCount=floor(TraLength.*0.03./R0.CellExtentInWorldY).*2+1;

% DOM data pre-process
ChannelCenterPixel=abs(floor((ChannelCenter-[R0.XWorldLimits(1,1), R0.YWorldLimits(1,2)])./[R0.CellExtentInWorldX, R0.CellExtentInWorldY]));
ChannelCenterPixelXBegin=ChannelCenterPixel(1,1)-floor(TraLength./R0.CellExtentInWorldX);
ChannelCenterPixelXEnd=ChannelCenterPixel(1,1)+floor(TraLength./R0.CellExtentInWorldX);
ChannelCenterPixelYBegin=ChannelCenterPixel(1,2)-floor(TraLength./R0.CellExtentInWorldY);
ChannelCenterPixelYEnd=ChannelCenterPixel(1,2)+floor(TraLength./R0.CellExtentInWorldY);
A1=A0(ChannelCenterPixelYBegin:ChannelCenterPixelYEnd, ChannelCenterPixelXBegin:ChannelCenterPixelXEnd);
% convert geotiff data to point location
xmin=R0.XWorldLimits(1,1);
ymax=R0.YWorldLimits(1,2);
dx=R0.CellExtentInWorldX;
dy=R0.CellExtentInWorldY;
pointCountX=ChannelCenterPixelXEnd-ChannelCenterPixelXBegin+1;
pointCountY=ChannelCenterPixelYEnd-ChannelCenterPixelYBegin+1;
pointCount=0;
pointLocationY=0;
pointIndex=0;
% divide the channel first
A2(pointCountY,pointCountX)=0;
ChannelGreyColorMin=ChannelGreyColor-ChannelGreyColorMix;
ChannelGreyColorMax=ChannelGreyColor+ChannelGreyColorMix;
for i=1:pointCountY
    for j=1:pointCountX
        if A1(i,j)>=ChannelGreyColorMin && A1(i,j)<=ChannelGreyColorMax
            A2(i,j) =255;
            pointCount=pointCount+1;
        else
            A2(i,j)=0;
        end
    end
end

% adapted channel extraction
A1T = adaptthresh(A1, 0.05,"Statistic", "gaussian", "NeighborhoodSize",2*floor(size(A1)/8)+1);
A1BW = imbinarize(A1,A1T);
% figure
% imshowpair(A1, A1BW, 'montage');
[A1L,n]=bwlabel(A1BW, 8);
A1S=regionprops(A1L, 'all');
A1max_area=max([A1S.Area]);
A1max_idx=find([A1S.Area]==A1max_area);
A1bb=A1S(A1max_idx).BoundingBox;
A1LongestComponent=imcrop(A1L,A1bb);
A1LongestComponent=A1LongestComponent~=A1max_idx;
A1LongestComponent=double(A1LongestComponent);
A1LNew=A1L~=A1max_idx;

% convert to point location
DOMPointLocation(pointCount,2)=0;
for i=1:pointCountY
    pointLocationY= ymax-dy*(i+ChannelCenterPixelYBegin-1);
    for j=1:pointCountX
        if A1LNew(i,j)==0
            pointIndex = pointIndex + 1;
            DOMPointLocation(pointIndex,1)=xmin+dx*(j+ChannelCenterPixelXBegin-1);
            DOMPointLocation(pointIndex,2)=pointLocationY;           
        end
    end
end
% downsample the DOM point
DOMDownSampleRate=floor(pointCount./(B0.Count.*5));
if DOMDownSampleRate==0
    DOMDownSampleRate=1;
end
DOMPointDSCount=size(DOMPointLocation,1)./DOMDownSampleRate;
DOMPointDSCount=floor(DOMPointDSCount);
DOMPointDS(DOMPointDSCount,2)=0;
for i=1:DOMPointDSCount
    DOMPointDS(i,1:2)=DOMPointLocation(i*DOMDownSampleRate,1:2);
end
% decentralize the DOM point
AvrDOM(1,1)=sum(DOMPointDS(:,1))./DOMPointDSCount;
AvrDOM(1,2)=sum(DOMPointDS(:,2))./DOMPointDSCount;
DOMPointDS(:,3)=DOMPointDS(:,1)-AvrDOM(1,1);
DOMPointDS(:,4)=DOMPointDS(:,2)-AvrDOM(1,2);

% calculate rough direction of the DOM and trajectory point, use DOMPointDS and B1
SumDirectionDOM(1,2)=0;
for i=1:(DOMPointDSCount-1)
    for j=(i+1):DOMPointDSCount
        DirectionVector=DOMPointDS(i,3:4)-DOMPointDS(j,3:4);
        SumDirectionDOM=SumDirectionDOM+DirectionVector;
    end
end
AvrDirectionDOM=SumDirectionDOM/(DOMPointDSCount.*(DOMPointDSCount-1)./2);
SumDirectionTra(1,2)=0;
for i=1:(B0.Count-1)
    for j=(i+1):B0.Count
        DirectionVector=B1(i,3:4)-B1(j,3:4);
        SumDirectionTra=SumDirectionTra+DirectionVector;
    end
end
AvrDirectionTra=SumDirectionTra/(B0.Count.*(B0.Count-1)./2);

% transform all point by parametric equation, DOMPointDS(:,3:4) and B1(:,3:4)
C1=DOMPointDS(:,3:4);
C1AvrDirection=AvrDirectionDOM(1,1:2);
C1k=C1AvrDirection(1,2)./C1AvrDirection(1,1);
C1Alpha=C1AvrDirection(1,1)./abs(C1AvrDirection(1,1));
C1Beta=C1AvrDirection(1,2)./abs(C1AvrDirection(1,2));
C1Gama=C1k./abs(C1k);
C1ParaT1=C1Alpha./sqrt (C1k.^2+1);
C1ParaT2=C1Beta.*C1Gama./sqrt (C1k.^2+1);
C1T(:,1)=C1ParaT1.*(C1(:,1)+C1k.*C1(:,2));
C1T(:,2)=C1ParaT2.*(C1k.*C1(:,1)-C1(:,2));
C1TSort=sortrows(C1T,1,"ascend");

C2=B1(:,3:4);
C2AvrDirection=AvrDirectionTra(1,1:2);
C2k=C2AvrDirection(1,2)./C2AvrDirection(1,1);
C2Alpha=C2AvrDirection(1,1)./abs(C2AvrDirection(1,1));
C2Beta=C2AvrDirection(1,2)./abs(C2AvrDirection(1,2));
C2Gama=C2k./abs(C2k);
C2ParaT1=C2Alpha./sqrt (C2k.^2+1);
C2ParaT2=C2Beta.*C2Gama./sqrt (C2k.^2+1);
C2T(:,1)=C2ParaT1.*(C2(:,1)+C2k.*C2(:,2));
C2T(:,2)=C2ParaT2.*(C2k.*C2(:,1)-C2(:,2));
C2TSort=sortrows(C2T,1,"ascend");

% iterate DOM point cloud by Fourier transform
C1Tdft_y0=C1TSort(:,1:2);
for i=1:DOMFourierSmoothIterationCount
    C1Tdft_yi= fft(C1Tdft_y0(:,2));
    C1Tdft_yi (15:end-15) = 0;
    C1Tsmoothed_yi = ifft(C1Tdft_yi);
    C1Tsmoothed_yi=real(C1Tsmoothed_yi);
    C1Tdft_x0=C1Tdft_y0(:,1);
    C1Tdft_y0(:,3)=(C1Tdft_y0(:,2)-C1Tsmoothed_yi).^2;
    C1Tdft_y0=sortrows(C1Tdft_y0,3,"ascend");
    C1Tdft_y0=C1Tdft_y0(1:end-5,1:2);
    C1Tdft_y0=sortrows(C1Tdft_y0,1,"ascend");
end

% normalize the DOM point
C3=[];
C3(:,1)=C1Tdft_x0;
C3(:,2)=C1Tsmoothed_yi;
C3Downsample=downsample(C3,floor(size(C3,1)./size(C2TSort,1)));
C3Downsample=sortrows(C3Downsample,1,"ascend");
C3DownsampleCount=size(C3Downsample,1);
for i=2:size(C3Downsample,1)
    C3DownsampleDx=C3Downsample(i,1)-C3Downsample(i-1,1);
    if C3DownsampleDx>1
        C3DownsampleCountPlus=floor(C3DownsampleDx);
        C3DownsampleDy=(C3Downsample(i,2)-C3Downsample(i-1,2))./C3DownsampleDx;
        for j=1:C3DownsampleCountPlus
            C3DownsampleCount=C3DownsampleCount+1;
            C3Downsample(C3DownsampleCount,1)=C3Downsample(i-1,1)+j;
            C3Downsample(C3DownsampleCount,2)=C3Downsample(i-1,2)+j.*C3DownsampleDy;
        end
    end
end
C3Downsample=sortrows(C3Downsample,1,"ascend");
C3DOMFinal=floor(C3Downsample);
C3DOMFinalXReference=C3DOMFinal(1,1);
for i=2:size(C3DOMFinal,1)
    C3DOMFinalX=C3DOMFinal(i,1);
    if C3DOMFinalX==C3DOMFinalXReference
        C3DOMFinal(i,1)=0;
        C3DOMFinal(i,2)=0;
    else
        C3DOMFinalXReference=C3DOMFinalX;        
    end
end
C3DOMFinal(all(C3DOMFinal==0,2),:) = [];
C1Tdft_yi= fft(C3DOMFinal(:,2));
C1Tdft_yi (15:end-15) = 0;
C1Tsmoothed_yi = ifft(C1Tdft_yi);
C3DOMFinal(:,2)=real(C1Tsmoothed_yi);

% divide accurate channel point
tree0 = KDTreeSearcher(C3DOMFinal);
C2DOMChannelCount=size(C1Tdft_y0,1);
C2DOMChannel=[];
C2DOMChannel(size(C1TSort,1),3)=0;
K = 1;
for i=1:size(C1TSort,1)
    queryPoint = [C1TSort(i,1), C1TSort(i,2)];
    [K, distances] = knnsearch(tree0, queryPoint);
    C2DOMChannel(i,1:2)=C1TSort(i,1:2);
    C2DOMChannel(i,3)=distances;
end
C2DOMChannel=sortrows(C2DOMChannel,3,"ascend");
C2DOMChannel(C2DOMChannelCount+1:end,:)=[];
C2DOMChannel1Count=2000;
C2DOMChannel1=downsample(C2DOMChannel,floor(size(C2DOMChannel,1)./2000));
C2DOMChannel1Count=size(C2DOMChannel1,1);

% smooth DOM point cloud by Fourier transform
C2DOMChannel=sortrows(C2DOMChannel,1,"ascend");
C1Tdft_yi= fft(C2DOMChannel(:,2));
C1Tdft_yi (15:end-15) = 0;
C1Tsmoothed_yi = ifft(C1Tdft_yi);
C1Tsmoothed_yi=real(C1Tsmoothed_yi);

% normalize the DOM point
C3=[];
C3(:,1)=C2DOMChannel(:,1);
C3(:,2)=C1Tsmoothed_yi;
C3Downsample=downsample(C3,floor(size(C3,1)./size(C2TSort,1)));
C3Downsample=sortrows(C3Downsample,1,"ascend");
C3DownsampleCount=size(C3Downsample,1);
for i=2:size(C3Downsample,1)
    C3DownsampleDx=C3Downsample(i,1)-C3Downsample(i-1,1);
    if C3DownsampleDx>1
        C3DownsampleCountPlus=floor(C3DownsampleDx);
        C3DownsampleDy=(C3Downsample(i,2)-C3Downsample(i-1,2))./C3DownsampleDx;
        for j=1:C3DownsampleCountPlus
            C3DownsampleCount=C3DownsampleCount+1;
            C3Downsample(C3DownsampleCount,1)=C3Downsample(i-1,1)+j;
            C3Downsample(C3DownsampleCount,2)=C3Downsample(i-1,2)+j.*C3DownsampleDy;
        end
    end
end
C3Downsample=sortrows(C3Downsample,1,"ascend");
C3DOMFinal=floor(C3Downsample);
C3DOMFinalXReference=C3DOMFinal(1,1);
for i=2:size(C3DOMFinal,1)
    C3DOMFinalX=C3DOMFinal(i,1);
    if C3DOMFinalX==C3DOMFinalXReference
        C3DOMFinal(i,1)=0;
        C3DOMFinal(i,2)=0;
    else
        C3DOMFinalXReference=C3DOMFinalX;        
    end
end
C3DOMFinal(all(C3DOMFinal==0,2),:) = [];
C1Tdft_yi= fft(C3DOMFinal(:,2));
C1Tdft_yi (15:end-15) = 0;
C1Tsmoothed_yi = ifft(C1Tdft_yi);
C3DOMFinal(:,2)=real(C1Tsmoothed_yi);

C3DOMFinalCov=C3DOMFinal(:,2);
C3DOMFinalCov=downsample(C3DOMFinalCov,5);
KNNTreeData=[];
KNNTreeData(1:(size(C3DOMFinal,1)+C2DOMChannel1Count),2)=0;
KNNTreeData(1:size(C3DOMFinal,1),1:2)=C3DOMFinal(:,1:2);
KNNTreeData((size(C3DOMFinal,1)+1):(size(C3DOMFinal,1)+C2DOMChannel1Count),1:2)=C2DOMChannel1(:,1:2);
KNNTreeData=sortrows(KNNTreeData,1,"ascend");
tree = KDTreeSearcher(KNNTreeData);
tree1 = KDTreeSearcher(C3DOMFinal);

% calculate the rough rotation parameter of trajectory
C3Tra1=sortrows(C2TSort,1,"ascend");
C3TraCount=size(C3Tra1,1);
C3TraRoughAngle=[0 10 20 30 140 160 170 180 190 200 210 330 340 350];
C3CovDirectionC=[];
for i=1:14
    C3TraRoughTheta=C3TraRoughAngle(i);
    C3TraRoughRot = [cosd(C3TraRoughTheta) -sind(C3TraRoughTheta); sind(C3TraRoughTheta) cosd(C3TraRoughTheta)];
    C3Tra2=(C3TraRoughRot*C3Tra1')';
    % standerdize the trajectory
    for j=2:size(C3Tra2,1)
        C3TraDx=C3Tra2(j,1)-C3Tra2(j-1,1);
        if C3TraDx>1
            C3TraCountPlus=floor(C3TraDx);
            C3TraDy=(C3Tra2(j,2)-C3Tra2(j-1,2))./C3TraDx;
            for k=1:C3TraCountPlus
                C3TraCount=C3TraCount+1;
                C3Tra2(C3TraCount,1)=C3Tra2(j-1,1)+k;
                C3Tra2(C3TraCount,2)=C3Tra2(j-1,2)+k.*C3TraDy;
            end
        end
    end
    C3Tra2=sortrows(C3Tra2,1,"ascend");
    C3TraFinal=floor(C3Tra2);
    C3TraFinalXReference=C3TraFinal(1,1);
    for j=2:size(C3TraFinal,1)
        C3TraFinalX=C3TraFinal(j,1);
        if C3TraFinalX==C3TraFinalXReference
            C3TraFinal(j,1)=0;
            C3TraFinal(j,2)=0;
        else
            C3TraFinalXReference=C3TraFinalX;        
        end
    end
    C3TraFinal(all(C3TraFinal==0,2),:) = [];
    C3TraFinal=sortrows(C3TraFinal,1,"descend");
    C3TraFinalCov=C3TraFinal(:,2);

    C3TraFinalCov=downsample(C3TraFinalCov,5);
    C3DOMTraDifSize=size(C3DOMFinalCov,1)-size(C3TraFinalCov,1)+1;
    C3DOMTraDif=[];
    C3DOMTraDif(C3DOMTraDifSize,1)=0;
    for j=1:C3DOMTraDifSize
        C3DOMTraDif(j,1)=norm(C3DOMFinalCov(j:(j+size(C3TraFinalCov,1)-1),1));
        C3DOMTraDif(j,2)=(j-1).*5+1;
    end
    C3CovDirectionA=conv(C3DOMFinalCov,C3TraFinalCov);
    C3CovDirectionB=(C3CovDirectionA(size(C3TraFinalCov,1):size(C3DOMFinalCov,1),1)./C3DOMTraDif(:,1))./norm(C3TraFinalCov);
    C3CovDirectionB(:,2)=C3DOMTraDif(:,2);
    C3CovDirectionB(:,3)=C3TraRoughAngle(i);
    C3CovDirectionC=[C3CovDirectionC; C3CovDirectionB];
end
C3CovDirectionC(:,4)=abs(C3CovDirectionC(:,1));
C3CovDirectionC=sortrows(C3CovDirectionC,4,"descend");
C3TraRoughTheta=C3CovDirectionC(1,3);
C3TraRoughRot = [cosd(C3TraRoughTheta) -sind(C3TraRoughTheta); sind(C3TraRoughTheta) cosd(C3TraRoughTheta)];

C2TSortNew=(C3TraRoughRot*C2TSort')';

% calculate the optimal transform parameters
ScalingCount=floor((ScalingInterval(1,2)-ScalingInterval(1,1))./ScalingStep)+1;
C3CovMax(ScalingCount,4)=0;
for i0=1:ScalingCount
    C3CovMax(i0,1)=ScalingInterval(1,1)+(i0-1).*ScalingStep; % scaling factor
    C3Tra=sortrows(C2TSortNew,1,"ascend");
    C3Tra=C3Tra(:,1:2).*C3CovMax(i0,1);
    C3TraCount=size(C3Tra,1);
    for i=2:size(C3Tra,1)
        C3TraDx=C3Tra(i,1)-C3Tra(i-1,1);
        if C3TraDx>1
            C3TraCountPlus=floor(C3TraDx);
            C3TraDy=(C3Tra(i,2)-C3Tra(i-1,2))./C3TraDx;
            for j=1:C3TraCountPlus
                C3TraCount=C3TraCount+1;
                C3Tra(C3TraCount,1)=C3Tra(i-1,1)+j;
                C3Tra(C3TraCount,2)=C3Tra(i-1,2)+j.*C3TraDy;
            end
        end
    end
    C3Tra=sortrows(C3Tra,1,"ascend");
    C3TraFinal=floor(C3Tra);
    C3TraFinalXReference=C3TraFinal(1,1);
    for i=2:size(C3TraFinal,1)
        C3TraFinalX=C3TraFinal(i,1);
        if C3TraFinalX==C3TraFinalXReference
            C3TraFinal(i,1)=0;
            C3TraFinal(i,2)=0;
        else
            C3TraFinalXReference=C3TraFinalX;        
        end
    end
    C3TraFinal(all(C3TraFinal==0,2),:) = [];
    C3TraFinal=sortrows(C3TraFinal,1,"descend");
    C3TraFinalCov=[];
    C3TraFinalCov=C3TraFinal(:,2);
    
    % find corresponding DOM points for trajectory
    C3TraFinalCov=downsample(C3TraFinalCov,5);
    C3DOMTraDifSize=size(C3DOMFinalCov,1)-size(C3TraFinalCov,1)+1;
    C3DOMTraDif=[];
    C3DOMTraDif(C3DOMTraDifSize,1)=0;
    for j=1:C3DOMTraDifSize
        C3DOMTraDif(j,1)=norm(C3DOMFinalCov(j:(j+size(C3TraFinalCov,1)-1),1));
        C3DOMTraDif(j,2)=(j-1).*5+1;
    end
    C3CovA=conv(C3DOMFinalCov,C3TraFinalCov);
    C3CovB=(C3CovA(size(C3TraFinalCov,1):size(C3DOMFinalCov,1),1)./C3DOMTraDif(:,1))./norm(C3TraFinalCov);
    C3CovB(:,2)=C3DOMTraDif(:,2);
    C3CovB(:,3)=abs(C3CovB(:,1));
    C3CovB=sortrows(C3CovB,3,"descend");
    C3CovMax(i0,2)=C3CovB(1,2); % corresponding DOM point start location
    C3CovMax(i0,3)=C3CovB(1,3); % correlation coefficient
    
    C3SRTTra=sortrows(C3TraFinal,1,"ascend");
    C3SRTTraBegin=1;
    C3SRTTraEnd=size(C3SRTTra,1);
    C3SRTDOMBegin=C3CovB(1,2);
    C3SRTDOMEnd=C3SRTDOMBegin+size(C3SRTTra,1)-1;
    if C3SRTDOMBegin<1
        C3SRTTraBegin=abs(C3SRTDOMBegin)+2;
        C3SRTDOMBegin=1;
    end
    if C3SRTDOMEnd>size(C3DOMFinal,1)
        C3SRTTraEnd=size(C3SRTTra,1)-(C3SRTDOMEnd-size(C3DOMFinal,1));
        C3SRTDOMEnd=size(C3DOMFinal,1);
    end
    C3SRTDOM=C3DOMFinal(C3SRTDOMBegin:C3SRTDOMEnd,:);
    C3SRTTra=C3SRTTra(C3SRTTraBegin:C3SRTTraEnd,:);
    C3SRTTraAvrX=sum(C3SRTTra(:,1))./size(C3SRTTra,1);
    C3SRTTraAvrY=sum(C3SRTTra(:,2))./size(C3SRTTra,1);
    C3SRTDOMAvrX=sum(C3SRTDOM(:,1))./size(C3SRTDOM,1);
    C3SRTDOMAvrY=sum(C3SRTDOM(:,2))./size(C3SRTDOM,1);
    C3SRTTraDC=[];
    C3SRTTraDC(:,1)=C3SRTTra(:,1)-C3SRTTraAvrX;
    C3SRTTraDC(:,2)=C3SRTTra(:,2)-C3SRTTraAvrY;
    C3SRTDOMDC=[];
    C3SRTDOMDC(:,1)=C3SRTDOM(:,1)-C3SRTDOMAvrX;
    C3SRTDOMDC(:,2)=C3SRTDOM(:,2)-C3SRTDOMAvrY;


    C3SRTS=(C3SRTTraDC')*C3SRTDOMDC;
    [C3SRTU, C3SRTE, C3SRTV] = svd(C3SRTS);
    C3SRTDetUV=det(C3SRTV*(C3SRTU'));
    C3SRTEE=eye(2);
    C3SRTEE(2,2)=C3SRTDetUV;
    C3SRTRotM=C3SRTV*C3SRTEE*(C3SRTU');

    C3SRTT=([C3SRTDOMAvrX C3SRTDOMAvrY]'-(C3SRTRotM*[C3SRTTraAvrX C3SRTTraAvrY]'))';
    C3SRTTraNew=(C3SRTRotM*(C3SRTTra(:,1:2)'))';
    C3SRTTraNew(:,1)=C3SRTTraNew(:,1)+C3SRTT(1,1);
    C3SRTTraNew(:,2)=C3SRTTraNew(:,2)+C3SRTT(1,2);
    C3SRTTraNewX=sum(C3SRTTraNew(:,1))./size(C3SRTTraNew,1);
    C3SRTTraNewY=sum(C3SRTTraNew(:,2))./size(C3SRTTraNew,1);
    
    C3SRTDOM1=[];
    C3SRTDOM1(1:size(C3SRTTraNew,1),2)=0;
    K = 1;
    for i=1:size(C3SRTTraNew,1)
        queryPoint = [C3SRTTraNew(i,1), C3SRTTraNew(i,2)];
        [K, distances] = knnsearch(tree, queryPoint);
        C3SRTDOM1(i,1:2)=KNNTreeData(K,1:2);
        %C3SRTDOM1(i,1:2)=C3DOMFinal(K,1:2);
    end
    C3SRTDOM1AvrX=sum(C3SRTDOM1(:,1))./size(C3SRTDOM1,1);
    C3SRTDOM1AvrY=sum(C3SRTDOM1(:,2))./size(C3SRTDOM1,1);

    C3SRTTraDC=[];
    C3SRTTraDC(:,1)=C3SRTTraNew(:,1)-C3SRTTraNewX;
    C3SRTTraDC(:,2)=C3SRTTraNew(:,2)-C3SRTTraNewY;
    C3SRTDOMDC=[];
    C3SRTDOMDC(:,1)=C3SRTDOM1(:,1)-C3SRTDOM1AvrX;
    C3SRTDOMDC(:,2)=C3SRTDOM1(:,2)-C3SRTDOM1AvrY;

    C3SRTS1=(C3SRTTraDC')*C3SRTDOMDC;
    [C3SRTU, C3SRTE, C3SRTV] = svd(C3SRTS1);
    C3SRTDetUV=det(C3SRTV*(C3SRTU'));
    C3SRTEE=eye(2);
    C3SRTEE(2,2)=C3SRTDetUV;
    C3SRTRotM1=C3SRTV*C3SRTEE*(C3SRTU');
    C3SRTT1=([C3SRTDOM1AvrX C3SRTDOM1AvrY]'-(C3SRTRotM1*[C3SRTTraNewX C3SRTTraNewY]'))';

    C3SRTTraNew1=(C3SRTRotM*(C3Tra(:,1:2)'))';
    C3SRTTraNew1(:,1)=C3SRTTraNew1(:,1)+C3SRTT(1,1);
    C3SRTTraNew1(:,2)=C3SRTTraNew1(:,2)+C3SRTT(1,2);
    C3SRTTraNew2=(C3SRTRotM1*(C3SRTTraNew1(:,1:2)'))';
    C3SRTTraNew2(:,1)=C3SRTTraNew2(:,1)+C3SRTT1(1,1);
    C3SRTTraNew2(:,2)=C3SRTTraNew2(:,2)+C3SRTT1(1,2);
    C3SRTTraNew2X=sum(C3SRTTraNew2(:,1))./size(C3SRTTraNew2,1);
    C3SRTTraNew2Y=sum(C3SRTTraNew2(:,2))./size(C3SRTTraNew2,1);
    
    C3SRTDOM2=[];
    C3SRTDOM2(1:size(C3SRTTraNew2,1),2)=0;
    K = 1;
    for i=1:size(C3SRTTraNew2,1)
        queryPoint = [C3SRTTraNew2(i,1), C3SRTTraNew2(i,2)];
        [K, distances] = knnsearch(tree, queryPoint);
        C3SRTDOM2(i,1:2)=KNNTreeData(K,1:2);
    end
    C3SRTDOM2AvrX=sum(C3SRTDOM2(:,1))./size(C3SRTDOM2,1);
    C3SRTDOM2AvrY=sum(C3SRTDOM2(:,2))./size(C3SRTDOM2,1);

    C3SRTTraNew2(:,1)=C3SRTTraNew2(:,1)+(C3SRTDOM2AvrX-C3SRTTraNew2X);
    C3SRTTraNew2(:,2)=C3SRTTraNew2(:,2)+(C3SRTDOM2AvrY-C3SRTTraNew2Y);

    C3SRTTraNew3Count=200;
    C3SRTTraNew3a=downsample(C3SRTTraNew2,floor(size(C3SRTTraNew2,1)./200));
    C3SRTTraNew3=C3SRTTraNew3a;
    C3SRTTraNew3Count=size(C3SRTTraNew3,1);
    C3SRTTraDOMDiffMin=[];
    C3SRTTraDOMDiffMin(RefineXCount.*RefineYCount,3)=0;
    RefineXMin=-1.*R0.CellExtentInWorldX.*(RefineXCount-1)./2;
    RefineYMin=-1.*R0.CellExtentInWorldY.*(RefineYCount-1)./2;
    for i=1:RefineXCount
        C3SRTTraNew3(:,1)=C3SRTTraNew3a(:,1)+RefineXMin+R0.CellExtentInWorldX.*i;
        for j=1:RefineYCount
            C3SRTTraNew3(:,2)=C3SRTTraNew3a(:,2)+RefineYMin+R0.CellExtentInWorldY.*j;
            C3SRTTraDOMDiff=0;
            for k=1:size(C3SRTTraNew3,1)
                queryPoint = [C3SRTTraNew3(k,1), C3SRTTraNew3(k,2)];
                [K, distances] = knnsearch(tree, queryPoint);
                C3SRTTraNew3DOM(k,1:2)=KNNTreeData(k,1:2);
                C3SRTTraDOMDiff=C3SRTTraDOMDiff+distances;
            end
            C3SRTTraDOMDiffMin((RefineXCount.*(i-1)+j),1)=i;
            C3SRTTraDOMDiffMin((RefineXCount.*(i-1)+j),2)=j;
            C3SRTTraDOMDiffMin((RefineXCount.*(i-1)+j),3)=C3SRTTraDOMDiff./size(C3SRTTraNew3,1);
        end
    end
    C3SRTTraDOMDiffMin=sortrows(C3SRTTraDOMDiffMin,3,"ascend");
    
    C3CovMax(i0,4)=C3SRTTraDOMDiffMin(1,3);
    C3CovMax(i0,5)=C3SRTT(1,1);
    C3CovMax(i0,6)=C3SRTT(1,2);
    C3CovMax(i0,7)=C3SRTRotM(1,1);
    C3CovMax(i0,8)=C3SRTRotM(1,2);
    C3CovMax(i0,9)=C3SRTRotM(2,1);
    C3CovMax(i0,10)=C3SRTRotM(2,2);
    C3CovMax(i0,11)=C3SRTT1(1,1)+(C3SRTDOM2AvrX-C3SRTTraNew2X);
    C3CovMax(i0,12)=C3SRTT1(1,2)+(C3SRTDOM2AvrY-C3SRTTraNew2Y);
    C3CovMax(i0,13)=C3SRTRotM1(1,1);
    C3CovMax(i0,14)=C3SRTRotM1(1,2);
    C3CovMax(i0,15)=C3SRTRotM1(2,1);
    C3CovMax(i0,16)=C3SRTRotM1(2,2);
    C3CovMax(i0,17)=RefineXMin+R0.CellExtentInWorldX.*C3SRTTraDOMDiffMin(1,1);
    C3CovMax(i0,18)=RefineYMin+R0.CellExtentInWorldY.*C3SRTTraDOMDiffMin(1,2);
end
C3CovMax=sortrows(C3CovMax,4,"ascend");

C3CovMaxIndex=1;
C3SRTTOpt(1,1)=C3CovMax(C3CovMaxIndex,5);
C3SRTTOpt(1,2)=C3CovMax(C3CovMaxIndex,6);
C3SRTRotMOpt=[C3CovMax(C3CovMaxIndex,7) C3CovMax(C3CovMaxIndex,8); C3CovMax(C3CovMaxIndex,9) C3CovMax(C3CovMaxIndex,10)];
C3SRTTraNewOpt=(C3SRTRotMOpt*(C3CovMax(C3CovMaxIndex,1).*C2TSortNew(:,1:2)'))';
C3SRTTraNewOpt(:,1)=C3SRTTraNewOpt(:,1)+C3SRTTOpt(1,1);
C3SRTTraNewOpt(:,2)=C3SRTTraNewOpt(:,2)+C3SRTTOpt(1,2);
C3SRTTOpt1(1,1)=C3CovMax(C3CovMaxIndex,11)+C3CovMax(C3CovMaxIndex,17);
C3SRTTOpt1(1,2)=C3CovMax(C3CovMaxIndex,12)+C3CovMax(C3CovMaxIndex,18);
C3SRTRotMOpt1=[C3CovMax(C3CovMaxIndex,13) C3CovMax(C3CovMaxIndex,14); C3CovMax(C3CovMaxIndex,15) C3CovMax(C3CovMaxIndex,16)];
C3SRTTraNewOpt1=(C3SRTRotMOpt1*C3SRTTraNewOpt')';
C3SRTTraNewOpt1(:,1)=C3SRTTraNewOpt1(:,1)+C3SRTTOpt1(1,1);
C3SRTTraNewOpt1(:,2)=C3SRTTraNewOpt1(:,2)+C3SRTTOpt1(1,2);

% Reproject on DOM
TraRotM1=AvrDirectionTra(1,1)./norm(AvrDirectionTra);
TraRotM2=AvrDirectionTra(1,2)./norm(AvrDirectionTra);
TraRotM=[TraRotM1 TraRotM2; -TraRotM2 TraRotM1];
B1(:,5:6)=(TraRotM*B1(:,3:4)')';
DOMRotM1=AvrDirectionDOM(1,1)./norm(AvrDirectionDOM);
DOMRotM2=AvrDirectionDOM(1,2)./norm(AvrDirectionDOM);
DOMRotM=[DOMRotM1 DOMRotM2; -DOMRotM2 DOMRotM1];
DOMPointDS(:,5:6)=(DOMRotM*DOMPointDS(:,3:4)')';
C3SRTTraRevised(:,1)=(C1ParaT2.*C3SRTTraNewOpt1(:,1)+C1ParaT1.*C1k.*C3SRTTraNewOpt1(:,2))./(C1ParaT1.*C1ParaT2.*(C1k.^2+1));
C3SRTTraRevised(:,2)=(C1ParaT2.*C1k.*C3SRTTraNewOpt1(:,1)-C1ParaT1.*C3SRTTraNewOpt1(:,2))./(C1ParaT1.*C1ParaT2.*(C1k.^2+1));
C3SRTTraRevised(:,3)=C3SRTTraRevised(:,1)+AvrDOM(1,1);
C3SRTTraRevised(:,4)=C3SRTTraRevised(:,2)+AvrDOM(1,2);
C3SRTTraRevised(:,5)=(C3SRTTraRevised(:,3)-xmin)./R0.CellExtentInWorldX;
C3SRTTraRevised(:,6)=(ymax-C3SRTTraRevised(:,4))./R0.CellExtentInWorldY;
C3SRTTraRevised(:,5)=floor(C3SRTTraRevised(:,5));
C3SRTTraRevised(:,6)=floor(C3SRTTraRevised(:,6));
C3SRTTraRevised(:,7)=floor(C3SRTTraRevised(:,1)./R0.CellExtentInWorldX)+ChannelCenterPixelXBegin;
C3SRTTraRevised(:,8)=floor(C3SRTTraRevised(:,2)./R0.CellExtentInWorldY)+ChannelCenterPixelYBegin;
A0c=[];
A0c=A0b;
for i=1:size(C3SRTTraRevised,1)
    A0c((C3SRTTraRevised(i,6)-3):(C3SRTTraRevised(i,6)+3),(C3SRTTraRevised(i,5)-3):(C3SRTTraRevised(i,5)-3))=255;
end

scatter(C3SRTTraNewOpt1(:,1),C3SRTTraNewOpt1(:,2));

scatter(C1TSort(:,1),C1TSort(:,2));
hold on
scatter(KNNTreeData(:,1),KNNTreeData(:,2));
scatter(C2TSort(:,1),C2TSort(:,2));
scatter(C2TSortNew(:,1),C2TSortNew(:,2));
scatter(C3SRTTraNewOpt(:,1),C3SRTTraNewOpt(:,2));
scatter(C3SRTTraNewOpt1(:,1),C3SRTTraNewOpt1(:,2));
scatter(C3SRTTraNew(:,1),C3SRTTraNew(:,2));

figure
imshow(A0c,[]);

[RMat, TVec, Tra1New] = TraSVDTransform('trajectory.pcd', C3SRTTraRevised(:,3:4), C3CovMax(C3CovMaxIndex,1));
figure
scatter3(Tra1New(:,1),Tra1New(:,2),Tra1New(:,3));
hold on
scatter(C3SRTTraRevised(:,3),C3SRTTraRevised(:,4));
GlobalMapZoomSmooth=pcread('GlobalMapzoomsmooth.pcd');
GlobalMapZoomSmoothNewLocation=GlobalMapZoomSmooth.Location;
GlobalMapZoomSmoothNewLocation(:,1:2)=(RMat*(GlobalMapZoomSmoothNewLocation(:,1:2)'))';
GlobalMapZoomSmoothNewLocation(:,1)=GlobalMapZoomSmoothNewLocation(:,1)+TVec(1,1);
GlobalMapZoomSmoothNewLocation(:,2)=GlobalMapZoomSmoothNewLocation(:,2)+TVec(1,2);
% GlobalMapZoomSmoothNew=pointCloud(single(GlobalMapZoomSmoothNewLocation),"Intensity",GlobalMapZoomSmooth.Intensity);
figure
pcshow(GlobalMapZoomSmoothNew);
pcwrite(GlobalMapZoomSmoothNew,'GlobalMapZoomSmoothNew0.pcd');

