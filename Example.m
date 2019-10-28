% s = '/home/spiegel/Desktop/Programs/SWMF/data/3d__var_1_e20150321-122300-000.out.cdf';
% original = bats('file',s);
% original.clearField({'b1x','b1y','b1z','e'})
% uni = original.toBatsUni(0.125,{},'xrange',[-40 0],'yrange',[-15 15],'zrange',[-15 15]);

load('/home/spiegel/Desktop/Programs/SWMF/data/uni.mat');

uni.plot('newfigure','quiver','zrange',0,'variable','u','increment',3)
