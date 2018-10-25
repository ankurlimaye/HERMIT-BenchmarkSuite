function fifsino(filename,figlabel)
[mat,dx,dy,xmin,ymin]=readfif(filename);
[X Y]=size(mat);
si(mat,4,1,xmin,xmin+(X-1)*dx,ymin,ymin+(Y-1)*dy,'theta','rho',figlabel,0);