function [mody]=PGBEMD_v5(Im1,sigma,denoise)
%INPUT PARAMETERS:
%Im1 - input image
%sigma - parameter used by BM3D denoising, default = 60
%denoise - set as 'bm3d' to perform BM3D denoising (recommended)
%-------------------CONVERTING TO GRAYSCALE--------------------------
[rows, columns, numberOfColorChannels] = size(Im1);
if numberOfColorChannels > 1
  Im1=rgb2gray(Im1);
end
%-----------------SCALING IMAGE TO 0-1 RANGE-------------------------
Im1=double(Im1);
Im1=Im1-min(min(Im1));
Im1=Im1*(1/(max(max(Im1))));
%-----------------------BM3D DENOISING-------------------------------
if strcmp(denoise,'bm3d')==1
    profile='np';
    %PSNR=0;
    %Im2=Im1;
    [PSNR,Im2]=BM3D(Im1,Im1,sigma,profile,0);
end
%----------------COARSE BACKGROUND FILTRATION------------------------
[M, N]=size(Im1);
Im4=SmoothMeanpgFast(Im2,0.25*M);
Im3=Im2-SmoothMeanpg(Im4,5);
%------------------DENSITY MAP CALCULATION---------------------------
density=Densitypg(Im3);
density=density*2+2;
figure, imagesc(density);
%--------------------PGBEMD DECOMPOSITION----------------------------
[mody]=PGBEMD(Im3,density);
mody{1}=mody{1}*(1/max(max(mody{1})));
end

%---------------------INTERNAL FUNCTIONS-----------------------------
function [efemd] = PGBEMD(Im,gestosc)
%main decomposition
[M,N] = size(Im);
MaxDil = imdilatepg(Im,gestosc); %MAX filtration
MinDil = -imdilatepg(-Im,gestosc); %MIN filtration
MaxMap = ~(Im - MaxDil); %binary map of maxima
MinMap = ~(Im - MinDil); %binary map of minima

NumbExtrema = round(0.5*(sum(sum(MaxMap)) + sum(sum(MinMap)))); %mean number of extrema
AvDist = round(sqrt(N*M/NumbExtrema)); %mean distance between extrema
maskDil=AvDist/2;%mask for dilation operation
SmoothMean=SmoothMeanpgFast(0.5*(MaxDil+MinDil),maskDil);%envelope smoothing operation

IMF = Im - SmoothMean;%extracting first mode
residue = SmoothMean;%saving the residue
efemd{1} = IMF;
efemd{2}=residue;
end

function[Im2]=imdilatepg(Im,gestosc)
    %grayscale dilation
    Im2=Im;
    [x0,y0] = size(Im);
    for i=1:x0
        for j=1:y0
            
            ex1=i-round((gestosc(i,j)-1)/2);
            ex2=i+round((gestosc(i,j)-1)/2);
            ey1=j-round((gestosc(i,j)-1)/2);
            ey2=j+round((gestosc(i,j)-1)/2);
            if ex1<1
                ex1=1;
            end
            if ex2>x0
                ex2=x0;
            end
            if ey1<1
                ey1=1;
            end
            if ey2>y0
                ey2=y0;
            end
            neighborhood=Im(ex1:ex2,ey1:ey2);
            Im2(i,j)=max(max(neighborhood));
        end
    end
end

function [Im2] = SmoothMeanpg(Im,masksize)
    %image smoothing
    Im2=Im;
    [x0,y0] = size(Im);
    for i=1:x0
        for j=1:y0
            ex1=i-round((masksize-1)/2);
            ex2=i+round((masksize-1)/2);
            ey1=j-round((masksize-1)/2);
            ey2=j+round((masksize-1)/2);
            if ex1<1
                ex1=1;
            end
            if ex2>x0
                ex2=x0;
            end
            if ey1<1
                ey1=1;
            end
            if ey2>y0
                ey2=y0;
            end
            neighborhood=Im(ex1:ex2,ey1:ey2);
            Im2(i,j)=sum(sum(neighborhood))/numel(neighborhood);
        end
    end
end

function [Im2] = SmoothMeanpgFast(Im,masksize)
    %fast image smoothing
    Im2=Im;
    [x0,y0] = size(Im);
    s1=5;%sampling rate
    s2=(s1+1)/2;%initial coordinate
    x1=s2;
    y1=s2;
    while x1<x0
        x1=x1+s1;%x1 and y1 are the final coordinates
    end
    while y1<y0
        y1=y1+s1;
    end

    for i=s2:s1:x1
        for j=s2:s1:y1
            ex1=i-round((masksize-1)/2);
            ex2=i+round((masksize-1)/2);
            ey1=j-round((masksize-1)/2);
            ey2=j+round((masksize-1)/2);
            if ex1<1
                ex1=1;
            end
            if ex2>x0
                ex2=x0;
            end
            if ey1<1
                ey1=1;
            end
            if ey2>y0
                ey2=y0;
            end
            neighborhood=Im(ex1:ex2,ey1:ey2);
            Im2(i,j)=sum(sum(neighborhood))/numel(neighborhood);
        end
    end
        %fulfilling the rest of the image
    for i=1:x0
        for j=1:y0
            k=i/s1;
            k=ceil(k);
            k=(k-1)*s1+s2;
            l=j/s1;
            l=ceil(l);
            l=(l-1)*s1+s2;
            Im2(i,j)=Im2(k,l);
        end
    end
    Im2=Im2(1:x0,1:y0);
    %density=density(1:x0,1:y0);
end

function [density] = Densitypg(Im)
    %calculation of the density map
    [x0,y0] = size(Im);
    density=zeros(x0,y0);
    maxIm=max(max(Im));
    minIm=min(min(Im));
    %maxDeviation=0.1*(maxIm-minIm);
    maxMaskSize=0.1*(x0+y0);
    
    mean1=sum(sum(Im))/numel(Im);
    meanmatrix=zeros(x0,y0);
    meanmatrix=meanmatrix+mean1;
    deviation=abs(Im-meanmatrix);
    maxDeviation=0.5*sum(sum(deviation))/numel(deviation);
    
    s1=5;%sampling rate
    s2=(s1+1)/2;%initial coordinates
    
    x1=s2;
    y1=s2;
    while x1<x0
        x1=x1+s1;
    end
    while y1<y0
        y1=y1+s1;
    end
    masksize=1;
    
    for i=s2:s1:x1
        for j=s2:s1:y1
            AvDeviation=AvDeviationpg(masksize,Im,i,j);
            %bigger and bigger neighborhood investigation
            if (AvDeviation<maxDeviation)
                if masksize>=maxMaskSize
                    masksize=masksize-2;
                end
                while(AvDeviation<maxDeviation)&&(masksize<maxMaskSize)
                masksize=masksize+2;
                AvDeviation=AvDeviationpg(masksize,Im,i,j);
                end
            %smaller and smaller neighborhood investigation 
            elseif (AvDeviation>=maxDeviation)
                masksize=masksize+2;
                while(AvDeviation>=maxDeviation)&&(masksize>1)
                masksize=masksize-2;
                AvDeviation=AvDeviationpg(masksize,Im,i,j);
                end
                masksize=masksize+2;
            end
            density(i,j)=masksize;
        end
    end
    %fulfilling the rest of the density map
    for i=1:x0
        for j=1:y0
            k=i/s1;
            k=ceil(k);
            k=(k-1)*s1+s2;
            l=j/s1;
            l=ceil(l);
            l=(l-1)*s1+s2;
            density(i,j)=density(k,l);
        end
    end
    density=density(1:x0,1:y0);
end

function [AvDeviation] = AvDeviationpg(masksize,Im,i,j)
%calculation of average deviation parameter
[x0,y0] = size(Im);
ex1=i-round((masksize-1)/2);
ex2=i+round((masksize-1)/2);
ey1=j-round((masksize-1)/2);
ey2=j+round((masksize-1)/2);
if ex1<1
    ex1=1;
end
if ex2>x0
    ex2=x0;
end
if ey1<1
    ey1=1;
end
if ey2>y0
    ey2=y0;
end
neighborhood=Im(ex1:ex2,ey1:ey2);
[m0,n0] = size(neighborhood);
mean=sum(sum(neighborhood))/numel(neighborhood);
meanmatrix=zeros(m0,n0);
meanmatrix=meanmatrix+mean;
neighborhood=abs(neighborhood-meanmatrix);
AvDeviation=sum(sum(neighborhood))/numel(neighborhood);
end