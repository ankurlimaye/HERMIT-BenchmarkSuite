function [mat,DeltaX,DeltaY,Xmin,Ymin]=readfif(filename)
%
SizeOfHeader = 236;

fin = fopen(filename,'r','ieee-le');
if fin==-1
  error(sprintf('fif-file %s does not exist',filename));
end  
header = fread(fin,1,'int');

if header~=17737,
  fclose(fin);
  fin = fopen(filename,'r','ieee-be');
  header = fread(fin,1,'int');
  if header~=17737,
    error('Error in fif-file');
  end
end

    
FileName = fread(fin,100,'char');
Dummy = fread(fin,88,'char');
Dummy = fread(fin,1,'int');
XSamples = fread(fin,1,'int');
YSamples = fread(fin,1,'int');
Dummy = fread(fin,1,'int');
Xmin = fread(fin,1,'float');
Ymin = fread(fin,1,'float');
DeltaX = fread(fin,1,'float');
DeltaY = fread(fin,1,'float');
Dummy = fread(fin,1,'float');
Dummy = fread(fin,1,'float');
Dummy = fread(fin,1,'float');

[mat,many] = fread(fin,XSamples*YSamples,'float'); 
if many~=XSamples*YSamples,
  error(sprintf('Could only read %i elements\n',many));
end
mat=reshape(mat,YSamples,XSamples);
mat=mat';
fclose(fin);