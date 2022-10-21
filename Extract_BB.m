% Fei Ma
function BB = Extract_BB(image_rgb)
%extract the binary image of tongue from the image in rgb type
A = image_rgb;
%mean shift filter%
%     mA = edison_wrapper(A, @RGB2Luv, 'steps', 1, 'SpatialBandWidth', 8, 'RangeBandWidth', 10);
%     A = Luv2RGB(mA);
grayA=rgb2gray(A);
grayA = adapthisteq(grayA);

%use wavelet filter to extract the shadow part
q=6;
[t, phi, p.csi, w, PHI_1, PSI_1, h0, h1, g] = Rational_wavelet_FB(q);
b = phi;
h = ftrans2(b);
I = imfilter(grayA,h);
level_i = graythresh(I);
I1 = im2bw(I,level_i);
I1=~I1;

se = strel('disk',8);
I1 = imdilate(I1,se);
I1=  imfill(I1,'holes');
I1=~I1;

BB=grayA;
idx=find(I1);
BB(idx)=255;

BB=rgb2hsv(A);
H=BB(:,:,1);
S=BB(:,:,2);
V=BB(:,:,3);

% image threshold using Otsu's method to convert to binary image
level_h = graythresh(H);
H1 = im2bw(H,level_h);

level_s = graythresh(S);
S1 = im2bw(S,level_s);

level_v = graythresh(V);
V1 = im2bw(V,level_v);

%and of H and S
AA = (H1) & (S1);

BB = mathprocess(AA);
BB = BB & I1;
end


function B = mathprocess(AA)

%dilate first
se = strel('disk',8);
AA = imdilate(AA,se);
%fill the holes
AA = imfill(AA,'holes');
%erode back to original shape
AA = imerode(AA,se);

%remove small pieces
AA = imclearborder(AA);
%find the biggest piece
CC = bwconncomp(AA);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
B=zeros(size(AA));
B(CC.PixelIdxList{idx}) = 1;
end
