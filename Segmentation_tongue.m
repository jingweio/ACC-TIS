function [xs,ys,B] = Segmentation_tongue
%% load image
[FileName,Images_path] = uigetfile({'*.jpg';'*.png';'*.JPG';'*.bmp';'*.pgm'},'Select the wound image');
image_rgb = imread(strcat(Images_path,FileName));

%% extract binary image
disp('Binary image extraction..');
BB = Extract_BB(image_rgb);
BB = erode_dilate(BB, 40);

%% extract the boundary and order it
disp('Build boundary..');
BW=bwmorph(BB,'remove',inf);
[y,x] = find(BW);
k = convhull(x,y);
xss = x(k); yss = y(k);
xss(end) = []; yss(end) = [];
[xss,yss] = build_boundary(xss,yss);
[xss,yss] = projection(x,y,xss,yss);

%% snake model refine the extracted boundary
disp('Snake model based boundary-refinement..');
[xs,ys,B] = Segmentation_snake(xss, yss, image_rgb);

%% display the result
BW=bwmorph(B,'remove',inf);
BW = imdilate(BW, strel('disk', 2));
image_r = image_rgb(:,:,1);
image_g = image_rgb(:,:,2);
image_b = image_rgb(:,:,3);
image_r(BW) = 255;
image_g(BW) = 255;
image_b(BW) = 255;
image_rgb(:,:,1) = image_r;
image_rgb(:,:,2) = image_g;
image_rgb(:,:,3) = image_b;
figure; imshow(image_rgb);
end

%% Toolbox built
% build a uniform boundary
function [x, y] = build_boundary(xs, ys)
% input: colum vector
xs_f = [xs(2:end); xs(1)];
ys_f = [ys(2:end); ys(1)];
ds = d(xs, ys, xs_f, ys_f);
d_min = min(ds);
nums = ceil(ds./d_min);
% calculate the weights
weights_f = repmat(1:max(nums), [length(nums),1])./(nums + 1);
weights_f(weights_f>=1) = nan;
weights = 1 - weights_f;
x = weights.*xs + weights_f.*xs_f;
y = weights.*ys + weights_f.*ys_f;
% output row vector
x = reshape(transpose(x), [1, numel(x)]);
y = reshape(transpose(y), [1, numel(y)]);
x(isnan(x)) = [];
y(isnan(y)) = [];
end

% project the rough boundary on initial boundary
function [xb,yb] = projection(x0,y0,xs,ys)
% make sure the projection is feasible
while length(xs) > 0.85*length(x0)
   [xs, ys] = half_boundary_points(xs, ys);
end
xb=[]; yb=[];
xs_f = [xs(2:end), xs(1)];
ys_f = [ys(2:end), ys(1)];
for k=1:length(xs)
    x1 = xs(k); y1 = ys(k);
    x2 = xs_f(k); y2 = ys_f(k);
    [x, y]=projected_points(x1, y1, x2, y2, x0, y0);
    %given not projected points
    if isempty(x)
        continue;
    end
    % select the nearest one as the desired point
    [~, index] = min(abs(((y2-y1).*x-(x2-x1).*y+y1*(x2-x1)-(y2-y1)*x1)./sqrt((y2-y1).^2+(x2-x1).^2)));
    xb = [xb, x(index(1))];
    yb = [yb, y(index(1))];
end
end

%locate the projected points
function [x, y] = projected_points(x1, y1, x2, y2, x_original, y_original)
% tangent vector
n=[y1-y2, x2-x1];
x3=n(1)+x1; y3=n(2)+y1;
if det([x1, y1, 1; x2, y2, 1; x3, y3, 1])<0
    n=-n;
end
% set up the perpenticular vectors
x12=x1+n(1);
y12=y1+n(2);
x22=x2+n(1);
y22=y2+n(2);
% get the projected points
x=[]; y=[];
for k=1:length(x_original)
    if det([x1, y1, 1; x12, y12, 1; x_original(k), y_original(k), 1])<0 &&det([x2, y2, 1; x22, y22, 1; x_original(k), y_original(k), 1])>0
        x(end+1)=x_original(k);
        y(end+1)=y_original(k);
    end
end
end

% dwindle the boundary point to the half
function [x,y] = half_boundary_points(xs,ys)
ind = 1:2:length(xs);
x = xs(ind);
y = ys(ind);
end

% erode and dilate morphologic process
function BB = erode_dilate(B, disk_size)
se=strel('disk', disk_size);
B1=imerode(B,se);
L = bwlabel(B1);
sta = regionprops(L,'Area');
area = cat(1,sta.Area);
index = find(area == max(area));
B2 = ismember(L,index);
BB = imdilate(B2,se);
end

% calculate distance between two points
function distance=d(x1,y1,x2,y2)
distance=sqrt((x1-x2).^2+(y1-y2).^2);
end
