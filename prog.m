s = '/home/spiegel/Desktop/Programs/SWMF/data/3d__var_1_e20150321-122300-000.out.cdf';

b = bats('file',s);

b.clearField({'b1x','b1y','b1z','e'})

uni = b.toBatsUni(0.5,{},'xrange',[-40 -10],'yrange',[-10 10],'zrange',[-5   5]);

%load('../data/test1.mat');
%reduced = uni.reduceDomain('xrange',[-40,0],'yrange',[-15 15],'zrange',[-15 15]);

load('/home/spiegel/Desktop/Programs/SWMF/data/uni.mat')

cut = 11;
X = squeeze(uni.Output.x(:,:,cut));
Y = squeeze(uni.Output.y(:,:,cut));

figure;
p = pcolor(X,Y,squeeze(uni.Output.ux(:,:,cut)));
set(p,'EdgeColor','none');

uni.calc_rhoU;
figure;
p = pcolor(X,Y,squeeze(uni.Derived.rhoUx(:,:,cut)));
set(p,'EdgeColor','none');

uni.calc_j;
figure;
subplot(221); p = pcolor(X,Y,squeeze(uni.Output.jx(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(222); p = pcolor(X,Y,squeeze(uni.Output.jy(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(223); p = pcolor(X,Y,squeeze(uni.Output.jz(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(224); p = pcolor(X,Y,squeeze(uni.Output.j(:,:,cut))); set(p,'EdgeColor','none'); colorbar;

uni.calc_jxb;
figure;
subplot(221); p = pcolor(X,Y,squeeze(uni.Derived.jxbx(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(222); p = pcolor(X,Y,squeeze(uni.Derived.jxby(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(223); p = pcolor(X,Y,squeeze(uni.Derived.jxbz(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(224); p = pcolor(X,Y,squeeze(uni.Derived.jxb(:,:,cut))); set(p,'EdgeColor','none'); colorbar;

uni.calc_temp;
figure;
p = pcolor(X,Y,squeeze(uni.Derived.Temp(:,:,cut))); set(p,'EdgeColor','none');

uni.calc_pb;
figure; p = pcolor(X,Y,squeeze(uni.Derived.Pb(:,:,cut))); set(p,'EdgeColor','none');

uni.calc_beta;
figure; p = pcolor(X,Y,squeeze(uni.Derived.Beta(:,:,cut))); set(p,'EdgeColor','none');

uni.calc_alfven;
figure; p = pcolor(X,Y,squeeze(uni.Derived.Alfven(:,:,cut))); set(p,'EdgeColor','none');

uni.calc_vth;
figure; p = pcolor(X,Y,squeeze(uni.Derived.Vth(:,:,cut))); set(p,'EdgeColor','none');

uni.calc_gyroradius;
figure; p = pcolor(X,Y,squeeze(uni.Derived.Gyroradius(:,:,cut))); set(p,'EdgeColor','none');

uni.calc_plasmafreq;
figure; p = pcolor(X,Y,squeeze(uni.Derived.PlasmaFrequency(:,:,cut))); set(p,'EdgeColor','none');

uni.calc_inertiallength;
figure; p = pcolor(X,Y,squeeze(uni.Derived.InertialLength(:,:,cut))); set(p,'EdgeColor','none');

uni.calc_electronVelocity;
figure;
subplot(221); p = pcolor(X,Y,squeeze(uni.Derived.Vex(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(222); p = pcolor(X,Y,squeeze(uni.Derived.Vey(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(223); p = pcolor(X,Y,squeeze(uni.Derived.Vez(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(224); p = pcolor(X,Y,squeeze(uni.Derived.Ve(:,:,cut))); set(p,'EdgeColor','none'); colorbar;


uni.calc_protonVelocity;
figure;
subplot(221); p = pcolor(X,Y,squeeze(uni.Derived.Vix(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(222); p = pcolor(X,Y,squeeze(uni.Derived.Viy(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(223); p = pcolor(X,Y,squeeze(uni.Derived.Viz(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(224); p = pcolor(X,Y,squeeze(uni.Derived.Vi(:,:,cut))); set(p,'EdgeColor','none'); colorbar;

uni.calc_gradP;
figure;
subplot(221); p = pcolor(X,Y,squeeze(uni.Derived.GradPx(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(222); p = pcolor(X,Y,squeeze(uni.Derived.GradPy(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(223); p = pcolor(X,Y,squeeze(uni.Derived.GradPz(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(224); p = pcolor(X,Y,squeeze(uni.Derived.GradP(:,:,cut))); set(p,'EdgeColor','none'); colorbar;

uni.calc_gradPb;
figure;
subplot(221); p = pcolor(X,Y,squeeze(uni.Derived.GradPbx(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(222); p = pcolor(X,Y,squeeze(uni.Derived.GradPby(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(223); p = pcolor(X,Y,squeeze(uni.Derived.GradPbz(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(224); p = pcolor(X,Y,squeeze(uni.Derived.GradPb(:,:,cut))); set(p,'EdgeColor','none'); colorbar;

uni.calc_divBB;
figure;
subplot(221); p = pcolor(X,Y,squeeze(uni.Derived.DivBBx(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(222); p = pcolor(X,Y,squeeze(uni.Derived.DivBBy(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(223); p = pcolor(X,Y,squeeze(uni.Derived.DivBBz(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(224); p = pcolor(X,Y,squeeze(uni.Derived.DivBB(:,:,cut))); set(p,'EdgeColor','none'); colorbar;

uni.calc_divRhoUU;
figure;
subplot(221); p = pcolor(X,Y,squeeze(uni.Derived.DivRhoUUx(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(222); p = pcolor(X,Y,squeeze(uni.Derived.DivRhoUUy(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(223); p = pcolor(X,Y,squeeze(uni.Derived.DivRhoUUz(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(224); p = pcolor(X,Y,squeeze(uni.Derived.DivRhoUU(:,:,cut))); set(p,'EdgeColor','none'); colorbar;

uni.calc_vorticity;
figure;
subplot(221); p = pcolor(X,Y,squeeze(uni.Derived.Vorticityx(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(222); p = pcolor(X,Y,squeeze(uni.Derived.Vorticityy(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(223); p = pcolor(X,Y,squeeze(uni.Derived.Vorticityz(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
subplot(224); p = pcolor(X,Y,squeeze(uni.Derived.Vorticity(:,:,cut))); set(p,'EdgeColor','none'); colorbar;
