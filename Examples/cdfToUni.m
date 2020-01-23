% Here is a script to convert cdf files (in Adaptive Grid format)
% into a uniform grid using interpolation (so that you can use the plotting method implemented in this repo).

% Path definitions (note the last /)
path_cdf = 'path/to/CDF/';
path_to_save = 'path/to/CDF/';

% Starting and finishing files
% if empty, all files will be converted
fi = '3d__var_1_e20150321-113000-000.out.cdf';
fe = '3d__var_1_e20150321-121000-000.out.cdf';

% Parameters for the interpolation/resize
clearfields = {'b1x','b1y','b1z','e'};
cellsize = 0.125;
var = {};
xrange = [-40 0];
yrange = [-20 20];
zrange = [-15 15];

% Get the cdf files names
filesdir = dir([path_cdf,'*.cdf']);
files = [];
for i = 1 : numel(filesdir)
  files = [files;filesdir(i).name];
end

if isempty(fi) & isempty(fe)
  idx_i = 1;
  idx_e = size(files,1);
else
  idx_i = find(strcmp(cellstr(files),fi));
  idx_e = find(strcmp(cellstr(files),fe));
end

for i = idx_i : idx_e
  disp(files(i,:));
  F = [path_cdf,files(i,:)];
  file_save = [path_to_save,files(i,:),'.mat'];
  original = bats('file',F);
  original.clearFields(clearfields);
  uni = original.toBatsUni(cellsize,var,'xrange',xrange,'yrange',yrange,'zrange',zrange);
  if 1
    save(file_save,'uni','-v7.3');
  end
  clear uni original F file_save
end
