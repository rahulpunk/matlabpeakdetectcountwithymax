function P1 = InterX1(L3,varargin)
clc;
clear all;
close all;

filename = 'ss.csv';
hold on
[num,~] = xlsread('ss.csv');
x1 = num(:,1) ; y1 =  num(:,2) ; 
[pks,locs]=findpeaks(-y1);
hd=findpeaks(-y1);
op=text(locs+0.2,-pks,num2str((1:numel(-pks))'));
go=findpeaks(-y1);


%postive peak

[peakValues1, indexes1] = findpeaks(y1);
tValues =( stem(indexes1,peakValues1,'b'));
[pks1,locs1]=findpeaks(y1);
hd1=findpeaks(y1);
op1=text(locs1+0.2,pks1,num2str((1:numel(pks1))'));

go1=findpeaks(y1);

[~,idx4] = max(go1);
maxvalx20 = x1(idx4) ;
t=max(y1);
[~,idx3] = max(y1);
maxvalx = x1(idx3) ;
op=text(locs+0.2,-pks,num2str((1:numel(-pks))'));
go=findpeaks(-y1);



[val44,idx44] = max(pks1) ; 
iwant44 = [go1(idx44) val44]
iwant55 =  [go(idx44-1) y1(idx44-1) ] ; % previous point of max 

iwant66 =  [go(idx44+1) y1(idx44+1)] ; % next point of max

dt1=(iwant55(:,1:2:end-1));
rd2=-dt1;
dt=(iwant66(:,1:2:end-1));
rd1=-dt;

rdia00=t-rd1;
rdia10=rdia00/3;




reco=rdia10+rd1;

%x3=[idx3-500,idx3+500];
%y3=[reco,reco];


x3=[0 500];
y3=[48 48];



L1 = [x1 y1]' ;

L3 = [x3 ; y3] ; 
plot(x1,y1);
hold on

plot(x3,y3);
 

P2 =InterX(L1,L3) ;
 plot(x1,y1,'r')




Ab0=plot(P2(1,:),P2(2,:),'*r');



[maxy,idx30] = max(y1(:));
plot(x1, y1, x1(idx30),y1(idx30),'pr')

legend('y(x)', sprintf('yawot test the Maximum of y = %0.3f',maxy))
axis([xlim    min(y1) max(y1)+1])

%postive peak

[peakValues1, indexes1] = findpeaks(y1);
tValues =( stem(indexes1,peakValues1,'b'));
[pks1,locs1]=findpeaks(y1);
hd1=findpeaks(y1);
op1=text(locs1+0.2,pks1,num2str((1:numel(pks1))'));

go1=findpeaks(y1);
%negative peak
[peakValues2, indexesOfPeaks2] = findpeaks(-y1);


dh=stem(indexesOfPeaks2, y1(indexesOfPeaks2), 'r');

[pks,locs]=findpeaks(-y1);
hd=findpeaks(-y1);

%logic for single peak

t=max(y1);

[maxvaly,idx3] = max(y1);
maxvalx = x1(idx3) ;






%imortant


t1=max(y1);



asaa=(t1-y1:x1:2) ;
answ=max(asaa);


[maxy,idx05] = max(y1(:));


[maxy2,idx2] = max(t1(:));


[df,dr] =findpeaks(-y1);



intersectdiff=diff(P2(1,:,1));
rst1=mean(intersectdiff);
intersectdiff1=diff(P2(1,:,1));

rst2=intersectdiff(:,1:2:end-1);
wiredia1=max(rst2);
open=max(rst1);
rst3=min(intersectdiff1);

gtp=max(rst3);

hold off

  function P = InterX(L1,varargin)
    

    error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';
    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
          
   
   function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
   end
   end
end  
