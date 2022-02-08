function [d1,figh]=RedChannelCorrection(d1)


green_mean = d1.ops.meanImg;
red_mean = d1.ops.meanImg_chan2;
z = d1.ops.meanImg_chan2_corrected;
Imgsize = size(green_mean);

%% corrected red image2
green_mean = green_mean(:);
red_mean =red_mean(:);z=z(:);

[red_mean_corrected,p] = fitmin(green_mean,read_mean,0);

figh=figure;hold on

if isfield(d1,'green')
    scatter(d1.green,d1.red,'.'); 
    scatter(d1.green,d1.red-polyval(p,d1.green),'x');
end

scatter(green_mean,red_mean,'k.');
scatter(green_mean,red_mean_corrected,'rx')
scatter(xData,yData,'rx')
plot(xf([1 end]),polyval(p,xf([1 end])),'r-','LineWidth',2)
scatter(green_mean,z,'bx');
% legend({'originial','corrected','corrected2','fitbaseline'})
d1.ops.meanImg_chan2_corrected2 = reshape(red_mean_corrected,...
    Imgsize(1),Imgsize(2));

figure;
subplot(121)

 maxpicR = prctile(d1.ops.meanImg_chan2_corrected(:), 99.5);%
 maxpicG = prctile(d1.ops.meanImg(:), 99.5);
 img(:,:,1)=d1.ops.meanImg_chan2_corrected/maxpicR;
 img(:,:,2)=d1.ops.meanImg/maxpicG;
 img(:,:,3)=0;
imshow(img)

subplot(122)

 maxpicR = prctile(d1.ops.meanImg_chan2_corrected2(:),99.5);%
 maxpicG = prctile(d1.ops.meanImg(:),99.5);
 img(:,:,1)=d1.ops.meanImg_chan2_corrected2/maxpicR;
 img(:,:,2)=d1.ops.meanImg/maxpicG;
 img(:,:,3)=0;
imshow(img)

%% calculate redcell based on corrected red image2

for k = 1:size(d1.iscell,1)
     linearidx = sub2ind(Imgsize,double(d1.stat{k}.ypix)',double(d1.stat{k}.xpix)');
    d1.redcell2(k,1:3)=[mean(d1.ops.meanImg(linearidx)),mean(d1.ops.meanImg_chan2_corrected(linearidx)),...
     mean(d1.ops.meanImg_chan2_corrected2(linearidx))];
end
%%

iscellidx = find(d1.iscell(:,1)==1&d1.redcell2(:,3)>0); 
%find real cells with positive red signal

x = d1.redcell2(iscellidx,1);%sum(d1.ops.meanImg(linearidx))
y = d1.redcell2(iscellidx,3);%sum(d1.ops.meanImg_chan2_corrected2(linearidx))

[p,S] = polyfit(x,y,1); 
[y_fit,delta] = polyval(p,x,S);

figure;
plot(x,y,'b.'); hold on
plot(x,y_fit,'r-')
plot(x,y_fit+2*delta,'m--',x,y_fit-2*delta,'m--')
title('Linear Fit of Data with 95% Prediction Interval')
legend('Data','Linear Fit','95% Prediction Interval')


isred = iscellidx(y>y_fit+2*delta)%find cells that is likely to be red
d1.redcell2(isred,4)=1;


