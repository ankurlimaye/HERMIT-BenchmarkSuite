function [H1,H2,H3]=si(matr,col,ydir,xmin,xmax,ymin,ymax,xlab,ylab,figlab,aspect,sigmin,sigmax) 
%The function SI(matr,col,ydir,xmin,xmax,ymin,ymax,
%                xlab,ylab,figlab,aspect,sigmin,sigmax) 
% makes an color coded image of the matrix 'matr' using the 
% colormap 'col'. (opX means optional at level X).
%
% 'col' (op1) -> 0:cool 1:jet 2:hot 3:gray 4:inv.gray 5:red/white/blue 6:bin
%                7:inv.binary. 8:mod hot  col<0=>-col is used with min=max.
% 'ydir' (op2) -> 1: normal  0:reverse. 
% 'xmin' (op3) -> min. 1. coordinate.
% 'xmax' (op3) -> max. 1. coordinate. Should be > 'xmin'. If not then dx=1.
% 'ymin' (op3) -> min. 2. coordinate.
% 'ymax' (op3) -> max. 2. coordinate. Should be > 'xmin'. If not then dy=1.
% 'xlab' (op4) -> Label for xaxis.
% 'ylab' (op4) -> Label for yaxis.
% 'figlab' (op4) -> Label for figure.
% 'aspect' (op5) -> 0:No aspect correction  1:Aspect correction on.
% 'sigmin' (op6) -> Show figure with values starting from 'sigmin'.
% 'sigmax' (op6) -> Show figure with values ending at 'sigmax'.
%                   If sigmax<=sigmin, these parameters will not be used.
% The function returns three handles. 1:figure  2/3:colorcolumn
% Original Peter Toft, Kim V. Hansen og Lars Risbo Okt 1993.
% Current revision Feb. 23 1996. IMM, Techn. Uni. of Denmark.
 
[cols,rows]=size(matr); 
colmax=15; 
rowmax=15;

if ((exist('sigmin')==0) | (exist('sigmax')==0))
  sigmax=max(max(matr));
  sigmin=min(min(matr));
  mi=sigmin; 
  ma=sigmax; 
else
  if (sigmin<sigmax)
    matr=(matr<=sigmax).*matr+(matr>sigmax).*sigmax;
    matr=(matr>=sigmin).*matr+(matr<sigmin).*sigmin;
  end
  mi=sigmin;
  ma=sigmax;
end
 
if (exist('col')==0) 
  col=1; 
end 

if (col<0)
  col=-col;
  if abs(sigmax)>abs(sigmin)
    sigmax=abs(sigmax);
    sigmin=-sigmax;
  else
    sigmax=abs(sigmin);
    sigmin=-sigmax;
  end
  mi=sigmin;
  ma=sigmax;
end

if (exist('ydir')==0) 
  ydir=1; 
end 
 
if (exist('aspect')==0) 
  aspect=0; 
end 
 
fig=figure;
 
set(fig,'PaperType','a4letter') 
matr2=matr'; 
 
if mi==ma, 
  ma=mi+1; 
end; 
 
matr2 = matr2 - mi; 
maxvalue = ma-mi; 
scale = 63.0 / maxvalue; 
c2=(ma-mi)*63*[0:63]'+mi; 
if col<4 
  matr2 = 1+matr2*scale; 
  c=[64:-1:1]'; 
else 
  matr2 = 64-matr2*scale; 
  c=[1:64]'; 
end 

if col==0, colormap(cool); end 
if col==1, colormap(jet); end 
if col==2, colormap(hot); end 
if col==3, colormap(ones(64,3)-gray); end
if col==4, colormap(gray); end
if col==5,
  rwb=zeros(64,3);
  for i=1:32,
    cl=(i-1)/31;
    rwb(i,1)=cl;
    rwb(i,2)=cl;
    rwb(i,3)=1.0;
  end
  for i=32:64,
    cl=(64-i)/32;
    rwb(i,1)=1.0;
    rwb(i,2)=cl;
    rwb(i,3)=cl;
  end
  colormap(rwb);
end
if (col==6) | (col==7),
  bina=zeros(64,3);
  for i=1:32,
    bina(i,1)=1.0;
    bina(i,2)=1.0;
    bina(i,3)=1.0;
  end
  for i=33:64,
    bina(i,1)=0.0;
    bina(i,2)=0.0;
    bina(i,3)=0.0;
  end
  if col==7,
    bina=1-bina;
  end
  colormap(bina);
end

if (col==8),
  modhot=flipud(hot(64));
  modhot(64,1)=0.4;
  modhot(64,2)=0.4;
  modhot(64,3)=0.4;
  colormap(modhot);
end

if (aspect==1)
  if rows>=cols
    maxyrange=0.75;
    maxxrange=3/4*(0.1+0.65*cols/rows);
  else
    maxxrange=0.75;
    maxyrange=3/4*(0.1+0.65*rows/cols);
  end;
else
  maxxrange=0.75;
  maxyrange=0.75;
end;
 
H1=axes('Position',[0.1 0.125 maxxrange maxyrange]); 
H2=axes('Position',[maxxrange+0.1625 0.125 0.05 0.75]); 
H3=axes('Position',[maxxrange+0.1625 0.125 0.05 0.75]);     

axes(H2),plot(c2); 
axes(H3),image(c); 

if (exist('xmin')*exist('xmax')*exist('ymin')*exist('ymax'))==0 
  xmin=0; 
  xmax=cols-1; 
  ymin=0; 
  ymax=rows-1;
else
  if xmin>=xmax 
    xmin=0;
    xmax=cols-1; 
  end
  end; 
  if ymin>=ymax 
    ymin=0;
    ymax=rows-1; 
  end; 
end; 
 
dx=(xmax-xmin)/(cols-1); 
dy=(ymax-ymin)/(rows-1); 
 
dxh=dx/2; 
dyh=dy/2; 
 
if (cols>colmax) 
  dxh=0; 
end 
 
if (rows>rowmax) 
  dyh=0; 
end 

%axes(H1),image([xmin-dxh,xmax+dxh],[ymin-dyh,ymax+dyh],matr2); 
axes(H1),image([xmin,xmax],[ymin,ymax],matr2); 
H4=get(fig,'CurrentAxes'); 

if (aspect==0) 
  if (cols<=colmax)
    set(H4,'XTick',[xmin:dx:xmax]); 
  end
  if (rows<=rowmax) 
    set(H4,'YTick',[ymin:dy:ymax]); 
  end; 
end;

if (exist('xlab')&exist('ylab')&exist('figlab'))==1 
  xlabel(xlab); 
  ylabel(ylab); 
  title(figlab); 
end 
 
set(H2,'Xtick',[]); 
set(H2,'XtickLabels',[]); 
set(H2,'Ydir','normal'); 
set(H2,'Ylim',[mi ma]); 
set(H3,'Xtick',[]); 
set(H3,'XtickLabels',[]); 
set(H3,'Ytick',[]); 
if ydir==1, 
  set(H1,'Ydir','normal'); 
else 
  set(H1,'Ydir','reverse'); 
end; 
if ((col==3) | (col==4)), 
  set(fig,'InvertHardcopy','off'); 
  set(fig,'Color',[1 1 1]); 
  set(H1,'Xcolor',[0 0 0]) 
  set(H1,'Ycolor',[0 0 0]) 
  set(H1,'Color',[0 0 0]) 
  set(H2,'Xcolor',[0 0 0]) 
  set(H2,'Ycolor',[0 0 0]) 
  set(H3,'Xcolor',[0 0 0]) 
  set(H3,'Ycolor',[0 0 0]) 
  set(get(H1,'Title'),'Color',[0 0 0]) 
end 
if col==3, colormap(gray); end
drawnow; 
