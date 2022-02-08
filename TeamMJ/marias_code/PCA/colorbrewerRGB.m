function [colors] = colorbrewerRGB(num, type)
% this funciton taken two inputs:
% NUM = how many color classes to return
% and
% TYPE = nature of your data (sequential, diverging, qualitative)
% it returns a set of colors defined in RGB space
% colors are grabbed from ColorBrewer2.org
% SAMPLE USAGE: 
% [colors] = colorbrewerRGB(3, 'qualitative')

switch type
    case 'qualitative'
        if num > 12
            error('Can return 12 RGB values for type qualitative')
        end
        colormap = [166, 206, 227; 31, 120, 180; 178, 223, 138;...
            51, 160, 44; 251, 154, 153; 227, 26, 28; 253, 191, 111;...
            255, 127, 0; 202, 178, 214; 106, 61, 154; 255, 255, 153;...
            177, 89, 40]/255;
        colors = colormap(1:num,:);    
    case 'diverging'
        if num > 12
            error('Can return 11 RGB values for type diverging')
        end
        colormap = [165, 0, 38; 215, 48, 39; 244, 109, 67; 253, 174, 97;...
            254, 224, 144; 255, 255, 191; 224, 243, 248; 171, 217, 233;...
            116, 173, 209; 69, 117, 180; 49, 54, 149]/255;
        colors = colormap(1:num,:);  
    case 'sequential'
        if num > 9
            error('Can return 9 RGB values for type sequential')
        end
        colormap = [255, 247, 251; 236, 226, 240; 208, 209, 230;...
            166, 189, 219; 103, 169, 207; 54, 144, 192; 2, 129, 138;...
            1, 108, 89; 1, 70, 54]/255;
        colors = colormap(1:num,:); 
end
    
    
end