% s = '/home/spiegel/Desktop/Programs/SWMF/data/3d__var_1_e20150321-122300-000.out.cdf';
% original = bats('file',s);
% original.clearFields({'b1x','b1y','b1z','e'})
% uni = original.toBatsUni(0.125,{},'xrange',[-40 0],'yrange',[-15 15],'zrange',[-15 15]);

load('/home/spiegel/Desktop/Programs/SWMF/data/uni.mat');

uni.plot('newfigure','slice','xslice',-40,'variable','p','log','colorrange',[2e-2 1e-1],'color','jet','alpha',0.7);
uni.plot('isosurface','variable','bx','level',0,'colorvariable','bz','colorrange',[-2 14],'color','parula', 'colorposition','right','xrange',[-40 -1],'zrange',[-3 3],'alpha',0.7);
uni.plot('slice','yslice',0,'variable','ux','colorrange',[-300 300],'color','autumn','alpha',0.7);
uni.plot('stream','variable','b','start',[-10.38 -2.75 0; -5.38 -2.75 0; -30.38 -2.75 0],'LineWidth',2);
uni.plotEarth;

clear positions;
positions(:,1) = [-18.12:0.125:-15.75];
positions(:,2) = -2.875 + zeros(size(positions,1),1);
positions(:,3) = 0 + zeros(size(positions,1),1);

uni.plot('newfigure','slice','zslice',0,'variable','bz','colorrange',[-2 14],'color','jet','alpha',0.6);
uni.plot('stream','variable','b','start',positions,'LineWidth',1);
%uni.plot('stream','variable','b','start',[-17 -2.875 0],'LineWidth',2,'colorvariable','bx','colorrange',[-10 10],'color','jet');
uni.plot('stream','variable','b','start',[-17 -2.875 0],'LineWidth',2,'colorvariable','bz','colorrange',[-0.1 0.1],'color','jet');
uni.plot('slice','yslice',-2.875,'variable','jxby','colorrange',[-0.01 0.01],'xrange',[-20 -15],'zrange',[-4 4],'yrange',[-3.0 -2.5],'alpha',0.6,'color','jet','colorposition','right','xlim',[-20 -15],'ylim',[-12 0],'zlim',[-4 4]); view([0 0]);

pos(:,1) = [-18.12:0.125:-15.75];
pos(:,2) = -2.875 + zeros(size(pos,1),1);
pos(:,3) = 0.6005 + zeros(size(pos,1),1);
SC = uni.getData(pos);
