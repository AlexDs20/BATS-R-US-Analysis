% s = '/home/spiegel/Desktop/Programs/SWMF/data/3d__var_1_e20150321-122300-000.out.cdf';
% original = bats('file',s);
% original.clearField({'b1x','b1y','b1z','e'})
% uni = original.toBatsUni(0.125,{},'xrange',[-40 0],'yrange',[-15 15],'zrange',[-15 15]);

load('/home/spiegel/Desktop/Programs/SWMF/data/uni.mat');

uni.plot('newfigure','slice','xslice',-40,'variable','p','log','colorrange',[2e-2 1e-1],'color','jet','alpha',0.7);

uni.plot('isosurface','variable','bx','level',0,'colorvariable','bz','colorrange',[-2 14],'color','parula', 'colorposition','right','xrange',[-40 -1],'zrange',[-3 3],'alpha',0.7);

uni.plot('slice','yslice',0,'variable','ux','colorrange',[-300 300],'color','autumn','alpha',0.7);

uni.plot('stream','variable','b','start',[-10.38 -2.75 0; -5.38 -2.75 0; -30.38 -2.75 0],'LineWidth',2);

uni.plotEarth;
