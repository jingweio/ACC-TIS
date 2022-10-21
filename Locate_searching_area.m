function [X1,Y1,X2,Y2,x_area,y_area,boolean]=Locate_searching_area(x0,y0,x01,y01,x02,y02,m_image,n_image)
%searching area
boolean=1;
%output XY column space;
x1=mean([x0,x01]);
y1=mean([y0,y01]);
x2=mean([x0,x02]);
y2=mean([y0,y02]);

% distance=distance/10;
%locate the edge point of the area;
n1=[y1-y0,x0-x1];
n2=[y2-y0,x0-x2];

x11=n1(1)+x1;
y11=n1(2)+y1;
x12=-n1(1)+x1;
y12=-n1(2)+y1;

x21=n2(1)+x2;
y21=n2(2)+y2;
x22=-n2(1)+x2;
y22=-n2(2)+y2;

if det([x1,y1,1;x2,y2,1;x11,y11,1])*det([x1,y1,1;x2,y2,1;x21,y21,1])<0
    tempx=x11;
    tempy=y11;
    x11=x12;
    y11=y12;
    x12=tempx;
    y12=tempy;
end

%determine the boundary of area points
% X=[x11,x12,x21,x22];
% Y=[y11,y12,y21,y22];
% k3=convhull(X,Y);
% X=X(k3);
% Y=Y(k3);

%left parts
X1=[x11;x21;x1;x2];
Y1=[y11;y21;y1;y2];
k1=convhull(X1,Y1);
X1=X1(k1);
Y1=Y1(k1);
X1(end)=[];
Y1(end)=[];

%right parts
X2=[x12;x22;x1;x2];
Y2=[y12;y22;y1;y2];
k2=convhull(X2,Y2);
X2=X2(k2);
Y2=Y2(k2);
X2(end)=[];
Y2(end)=[];


%get a roungh xy by using the max and min of xy
XB1=min(X1):max(X1);
YB1=min(Y1):max(Y1);
XB2=min(X2):max(X2);
YB2=min(Y2):max(Y2);

[XBB1,YBB1]=meshgrid(XB1,YB1);
[XBB2,YBB2]=meshgrid(XB2,YB2);

%change!!!
XBB1=reshape(XBB1,numel(XBB1),1);
YBB1=reshape(YBB1,numel(YBB1),1);
XBB2=reshape(XBB2,numel(XBB2),1);
YBB2=reshape(YBB2,numel(YBB2),1);

%X1 Y1 is the row type
[in,on]=inpolygon(XBB1,YBB1,X1',Y1');
X1=XBB1(in);
Y1=YBB1(in);

%X2 Y2 is the row type
[in,on]=inpolygon(XBB2,YBB2,X2',Y2');
X2=XBB2(in);
Y2=YBB2(in);

[X1,Y1] = correction(X1,Y1,m_image,n_image);
[X2,Y2] = correction(X2,Y2,m_image,n_image);

x_area=vertcat(X1,X2);
y_area=vertcat(Y1,Y2);

%swith XY to integer for saving time
X1=round(X1);
Y1=round(Y1);
[X1,Y1]=removeTheRepeated(X1,Y1);
%in case of out of image
[X1,Y1] = correction(X1,Y1,m_image,n_image);

X2=round(X2);
Y2=round(Y2);
[X2,Y2]=removeTheRepeated(X2,Y2);
[X2,Y2] = correction(X2,Y2,m_image,n_image);

%remove the center point
overlap_x=x_area==x0;
overlap_y=y_area==y0;
overlap=overlap_x.*overlap_y;
x_area(find(overlap))=[];
y_area(find(overlap))=[];

%judge the existence of those area points
if isempty(X1)||isempty(X2)
    boolean=0;
end
end

function [x,y] = correction(xt,yt,m,n)
% [m,n] = size(image);
wrong1 = xt>n|xt<0;
wrong2 = yt>m|yt<0;
wrong = or(wrong1,wrong2);
xt(wrong) = [];
yt(wrong) = [];
x = xt;
y = yt;
end

%remove the repeated points 
function [x,y]=removeTheRepeated(x0,y0)
x0 = reshape(x0, [numel(x0), 1]);
y0 = reshape(y0, [numel(y0), 1]);
XY = [x0,y0];
XY_unique = unique(XY,'rows');
x = XY_unique(:,1);
y = XY_unique(:,2);
end
