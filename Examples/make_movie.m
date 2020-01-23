% Script to make an animation out of the files that have been converted into a uniform grid (e.g. using cdfToUni script).

% note the last /
path_uni = 'path/to/directory/containg/all/the/batsUni/files/';

filesdir = dir([path_uni,'*.mat']);
files = [];
for i = 1 : numel(filesdir)
  files = [files;filesdir(i).name];
end

fi = '3d__var_1_e20150321-115100-000.out.cdf.mat';
fe = '3d__var_1_e20150321-115100-000.out.cdf.mat';

vid = false;          % if you want to make a video of it, change it to yes

if isempty(fi) & isempty(fe)
  idx_i = 1;
  idx_e = size(files,1);
else
  idx_i = find(strcmp(cellstr(files),fi));
  idx_e = find(strcmp(cellstr(files),fe));
end

cmap1 = multigradient([0 0 0; 1 0 0; 1 1 0]);
cmap2 = multigradient([1 1 1; 0 1 1; 0 0 1]);

for i = idx_i : idx_e               % This loop can be changed into a parfor loop (however when I do that, the saved frames suck... but it works well if you just want to save the image directly.)
  disp(files(i,:));
  data(i) = load([path_uni,files(i,:)]);
  h(i) = ...
  data(i).uni.plot('newfigure',...
          'slice',...
          'zslice',0,...
          'variable','bz',...
          'colorrange',[-2 20],...
          'colorposition','westoutside',...
          'color',cmap1,...
          'alpha',0.9,...
          'xlim',[-40 -5],...
          'ylim',[-15  15] ...
          );

  data(i).uni.plot('quiver',...
          'variable','u', ...
          'xrange',[-30 -10], ...
          'yrange',[-10 10], ...
          'zrange',0.125,...
          'increment',4,...
          'color',cmap2,...
          'colorrange',[0 800],...
          'linewidth',0.5,...
          'HeadAngle',60,...
          'HeadLength',0.25);

  view(h(i).Children(end),[0 90]);

  title(h(i).Children(end),[files(i,21:22),':',files(i,23:24)]);
  set(h(i).Children(end),'xdir','reverse','ydir','reverse');
  drawnow;

  if vid
    I{i} = export_fig(h(i),'-q101');
    F(i) = im2frame(I{i});
  end

  if 0          % if you want to trace a field line at a certain point.
    [x y] = ginput(1);
    z = 0;
    data(i).uni.plot('stream','variable','b','start',[x y z],'LineWidth',2,'color',[0 0 0]);
  end

  if ~isempty(get(groot,'Children'))
    close(h(i));
  end
  data(i).uni=[];
end

if vid
  video = VideoWriter('tail.avi', 'Uncompressed AVI');
  video.FrameRate = 2;
  open(video);
  writeVideo(video, F(idx_i:idx_e));
  close(video);
end
