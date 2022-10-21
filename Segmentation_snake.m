function [xs,ys,B]=Segmentation_snake(x,y,image_rgb)
% row vector to colum vector
xss=x';
yss=y';
image=rgb2gray(image_rgb);
%% initial value preparation
%emperical coefficient for external energy
wl=0.22;
we=1.75;
wt=0.27;
%emperical coefficient for internal energy
sigma=1.00;
alpha=0.40;
beta =0.20;
iterations=200;
N=iterations;
%find the image energy
f=external_energy(sigma,wl,we,wt,image);

% half the length of boundary points
[xss, yss] = half_boundary_points(xss,yss);
[xss, yss] = half_boundary_points(xss,yss);
[xs, ys] = half_boundary_points(xss,yss);
disp([num2str(length(xs)), 'points..']);

f_Max=zeros(size(xs));
[m_image,n_image] = size(image);

%% iterations updating boundary
for i=1:N
    disp([num2str(i), 'th iteration..']);
    
    len=length(xs);
    Xs=zeros([len,1]); Ys=zeros([len,1]);
    xs_f = [xs(2:end); xs(1)];
    ys_f = [ys(2:end); ys(1)];
    xs_l = [xs(end); xs(1:end-1)];
    ys_l = [ys(end); ys(1:end-1)];
    for k=1:len
        %initial region points
        x0=xs(k);
        y0=ys(k);
        x1=xs_l(k);
        y1=ys_l(k);
        x2=xs_f(k);
        y2=ys_f(k);
        
        %locate the searching region
        [X1,Y1,X2,Y2,x_area,y_area,boolean]=Locate_searching_area(x0,y0,x1,y1,x2,y2,m_image,n_image);
        
        %xy_area without center point
        x_area=vertcat(x_area,x0);
        y_area=vertcat(y_area,y0);
        
        %xy area missing
        if ~boolean
            Xs(k) = x0;
            Ys(k) = y0;
            continue;
        end
        
        %ensure not repeatition
        if k>1
            x_overlapping=(x_area==Xs(k-1));
            y_overlapping=(y_area==Ys(k-1));
            index_overlapping=logical(x_overlapping.*y_overlapping);
            x_area(index_overlapping)=[];
            y_area(index_overlapping)=[];
        end
        
        %similar value
        similar=similarity_metric(X1,Y1,X2,Y2,image_rgb,image);
        
        %internal energy as a coeffcient
        line_energy_initial=abs(find_line_energy(alpha,beta,x0,y0,x1,y1,x2,y2));
        
        %energy term
        energy_line=abs(find_line_energy(alpha,beta,x_area,y_area,x1,y1,x2,y2));
        energy_image=abs(interp2(f,x_area,y_area));
        
        energy_image_mean=mean(energy_image);
        
        %find the extrem value for all iterations
        energy_image_Max=max(energy_image);
        if f_Max(k)<energy_image_Max
            f_Max(k)=energy_image_Max;
        end
        
        %balance energy stage of energy
        c1=removeExtremValue(energy_line);
        c2=removeExtremValue(energy_image);
        c1_mean=mean(c1);
        c2_mean=mean(c2);
        c1_mean_times=find_magnitude(c1_mean);
        c2_mean_times=find_magnitude(c2_mean);
        balance_energy=c1_mean_times/c2_mean_times;
        
        %energy factors
        gamma1=exp(similar);
        gamma2=line_energy_initial;
        
        gamma=gamma1*gamma2;
        gamma_times=find_magnitude(gamma2);
        
        kappa1=(abs(f_Max(k)-min(energy_image))/energy_image_mean);
        kappa2=exp(abs(1-similar));
        
        kappa=kappa1*kappa2;
        kappa_times=find_magnitude(kappa1);
        if similar==1
            kappa=1;
            kappa_times=1;
        end
        
        %balance energy stage of facter
        gamma=(gamma)/gamma_times;
        kappa=balance_energy*((kappa)/kappa_times);
        
        %energy internal and external
        energy=gamma*energy_line+kappa*energy_image;
        
        %optimaize the points in area consider continous power in the end
        index_Min=find(energy==min(energy));
        Xs(k)=x_area(index_Min(1));
        Ys(k)=y_area(index_Min(1));
    end
    
    %get the values
    xs=Xs;
    ys=Ys;
end

%enlarge
[xss,yss] = build_closed_boundary(xs,ys);
xss = round(xss); yss = round(yss);
idx = sub2ind(size(image),yss,xss);
box = zeros(size(image));
box(idx) = 1;
B = imfill(box,'holes');
end

%% Toolbox built
function [x, y] = build_closed_boundary(xs, ys)
% input: colum vector
xs_f = [xs(2:end); xs(1)];
ys_f = [ys(2:end); ys(1)];
ds = d(xs, ys, xs_f, ys_f);
d_min = 1; % make sure the neighbor points are neighbor pixels
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

% dwindle the boundary point to the half
function [x,y] = half_boundary_points(xs,ys)
ind = 1:2:length(xs);
x = xs(ind);
y = ys(ind);
end

function value=find_magnitude(data)
value = 10^(floor(log10(data)));
end

function f=external_energy(sigma,wl,we,wt,image)
smask = fspecial('gaussian', ceil(3*sigma), sigma);
smth = filter2(smask, image, 'same');

[row, col] = size(image);

%Computing external forces

eline = smth; %eline is simply the image intensities
[grady,gradx] = gradient(smth);
eedge = -1 * sqrt ((gradx .* gradx + grady .* grady)); %eedge is measured by gradient in the image

%masks for taking various derivatives
m1 = [-1 1];
m2 = [-1;1];
m3 = [1 -2 1];
m4 = [1;-2;1];
m5 = [1 -1;-1 1];

cx = conv2(smth,m1,'same');
cy = conv2(smth,m2,'same');
cxx = conv2(smth,m3,'same');
cyy = conv2(smth,m4,'same');
cxy = conv2(smth,m5,'same');

for i = 1:row
    for j= 1:col
        % eterm as deined in Kass et al Snakes paper
        eterm(i,j) = (cyy(i,j)*cx(i,j)*cx(i,j) -2 *cxy(i,j)*cx(i,j)*cy(i,j) + cxx(i,j)*cy(i,j)*cy(i,j))/((1+cx(i,j)*cx(i,j) + cy(i,j)*cy(i,j))^1.5);
    end
end

eext = (wl*eline + we*eedge -wt * eterm); %eext as a weighted sum of eline, eedge and eterm
f=eext;
end

function continus_power=find_line_energy(alpha,beta,x0,y0,x1,y1,x2,y2)
energy1=(x2-x0).^2+(y2-y0).^2;
energy2=(x2-2*x0+x1).^2+(y2-2*y0+y1).^2;
continus_power=alpha.*energy1+beta.*energy2;
end

function energy=removeExtremValue(energy_initial)
energy=energy_initial;
if length(energy)>3
    [id,c]=kmeans(energy_initial,3);
    [cs,index_c_sort]=sort(c);
    index_Medium=index_c_sort(2);
    id_energy=find(id==index_Medium);
    energy=energy_initial(id_energy);
end
end

function similar=similarity_metric(x_left,y_left,x_right,y_right,image_rgb,image)
%output column space ranging from 0 to 1 but also be closed to 1;
%can be improved here;
similar=[];
r=image_rgb(:,:,1);
g=image_rgb(:,:,2);
b=image_rgb(:,:,3);
r=double(r);
g=double(g);
b=double(b);

%sub to ind transformation
index_left=sub2ind(size(image),y_left,x_left);
index_right=sub2ind(size(image),y_right,x_right);

%r value
r_left=r(index_left);
r1=mean(r_left);

r_right=r(index_right);
r2=mean(r_right);
rt=r1/r2;

%g value
g_left=g(index_left);
g1=mean(g_left);

g_right=g(index_right);
g2=mean(g_right);
gt=g1/g2;

%b value
b_left=b(index_left);
b1=mean(b_left);

b_right=b(index_right);
b2=mean(b_right);
bt=b1/b2;
similar(end+1)=3^2/((rt+gt+bt)*(1/rt+1/gt+1/bt));
end

% calculate distance between two points
function distance=d(x1,y1,x2,y2)
distance=sqrt((x1-x2).^2+(y1-y2).^2);
end
