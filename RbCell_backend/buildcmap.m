function [cmap]=buildcmap(colors)
% [cmap]=buildcmap(colors)
%
% This function can be used to build your own custom colormaps. Imagine if
% you want to display rainfall distribution map. You want a colormap which
% ideally brings rainfall in mind, which is not achiveved by colormaps such
% as winter, cool or jet and such. A gradient of white to blue will do the
% task, but you might also use a more complex gradient (such as
% white+blue+red or colors='wbr'). This function can be use to build any
% colormap using main colors rgbcmyk. In image processing, w (white) can be
% used as the first color so that in the output, the background (usually
% with 0 values) appears white. In the example of rainfall map, 'wb' will
% produce a rainfall density map where the background (if its DN values are
% 0) will appear as white.
%
% Inputs:
%  colors: string (char) of color codes, any sequence of rgbcmywk
%  representing different colors (such as 'b' for blue) is acceptable. If a
%  gradient of white to blue is needed, colors would be 'wb'; a rainbow of
%  white+blue+red+green would be 'wbrg'.
%
% Example:
%  [cmap]=buildcmap('wygbr');
% %try the output cmap:
% im=imread('cameraman.tif');
% imshow(im), colorbar
% colormap(cmap) %will use the output colormap
%
% First version: 14 Feb. 2013
% sohrabinia.m@gmail.com
%--------------------------------------------------------------------------
plotcolors;
if nargin<1
    colors='wrgbcmyk';
end

if ~ischar(colors)
    error(['Error! colors must be a variable of type char with '...
        'color-names, such as ''r'', ''g'', etc., '...
        'type ''help buildcmap'' for more info']);
end

ncolors=length(colors)-1;


bins=round(255/ncolors);
% diff1=255-bins*ncolors;

vec=zeros(300,3);

switch colors(1)
    case 'w'
        vec(1,:)=1;
    case 'r'
        vec(1,:)=[1 0 0];
    case 'g'
        vec(1,:)=[0 1 0];
    case 'b'
        vec(1,:)=[0 0 1];
    case 'c'
        vec(1,:)=[0 1 1];
    case 'm'
        vec(1,:)=[1 0 1];
    case 'y'
        vec(1,:)=[1 1 0];
    case 'k'
        vec(1,:)=[0 0 0];
    case 'C'
        vec(1,:)=cCyan;
    case 'R'
        vec(1,:)=cRed;
    case 'G'
        vec(1,:)=cGreen;
    case 'Y'
        vec(1,:)=cYellow;
    case 'B'
        vec(1,:)=cBlue;
    case 'O'
        vec(1,:)=cOrange;
    case 'P'
        vec(1,:)=cPurple;
    case 'L'
        vec(1,:)=cLilac;
    case 'p'
        vec(1,:)=cPink;
    case 't'
        vec(1,:)=clGreen;
end

for i=1:ncolors
 beG=(i-1)*bins+1;
 enD=i*bins+1; %beG,enD
 switch colors(i+1)
     case 'w'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD,
     case 'r'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';%colors(i+1),beG,enD
     case 'g'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';%colors(i+1),beG,enD
     case 'b'         
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD
     case 'c'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD
     case 'm'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';
     case 'y'
         vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
     case 'k'
         vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
         
     case 'C'
         vec(beG:enD,1)=linspace(vec(beG,1),cCyan(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),cCyan(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),cCyan(3),bins+1)';%colors(i+1),beG,enD,
     case 'R'
         vec(beG:enD,1)=linspace(vec(beG,1),cRed(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),cRed(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),cRed(3),bins+1)';%colors(i+1),beG,enD
     case 'G'
         vec(beG:enD,1)=linspace(vec(beG,1),cGreen(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),cGreen(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),cGreen(3),bins+1)';%colors(i+1),beG,enD
     case 'Y'         
         vec(beG:enD,1)=linspace(vec(beG,1),cYellow(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),cYellow(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),cYellow(3),bins+1)';%colors(i+1),beG,enD
     case 'B'
         vec(beG:enD,1)=linspace(vec(beG,1),cBlue(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),cBlue(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),cBlue(3),bins+1)';%colors(i+1),beG,enD
     case 'O'
         vec(beG:enD,1)=linspace(vec(beG,1),cOrange(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),cOrange(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),cOrange(3),bins+1)';
     case 'P'
         vec(beG:enD,1)=linspace(vec(beG,1),cPurple(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),cPurple(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),cPurple(3),bins+1)';
     case 'L'
         vec(beG:enD,1)=linspace(vec(beG,1),cLilac(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),cLilac(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),cLilac(3),bins+1)';
     case 'p'
         vec(beG:enD,1)=linspace(vec(beG,1),cPink(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),cPink(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),cPink(3),bins+1)';
     case 't'
         vec(beG:enD,1)=linspace(vec(beG,1),clGreen(1),bins+1)';
         vec(beG:enD,2)=linspace(vec(beG,2),clGreen(2),bins+1)';
         vec(beG:enD,3)=linspace(vec(beG,3),clGreen(3),bins+1)';
  end
end
cmap=vec(1:bins*ncolors,:);
end %end of buildcmap
