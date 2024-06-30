% default directory
dirlists = '/home/inwoo/data/Basin/imbie/';

% set exp name for exporting
expname = [dirlists './basinNumbers_8km.exp'];

% get nc file.
nc = [dirlists '/../imbie/basinNumbers_8km.nc'];
x = ncread(nc,'x');
y = ncread(nc,'y');
basinNumber = ncread(nc,'basinNumber');

fprintf('%s\n',nc);

% extend domain for contour.
% 100 values for outline.
basinNumber_ = 100*ones(size(basinNumber)+2); 
basinNumber_(2:end-1,2:end-1) = basinNumber;
dx = unique(diff(x));
dy = unique(diff(y));
x_ = [x(1)-dx; x; x(end)+dx];
y_ = [y(1)-dy; y; y(end)+dy];

% get min and max value for basin.
basinmin = min(basinNumber,[],'all');
basinmax = max(basinNumber,[],'all');
level = [basinmin:basinmax];

% contour_map structure.
contour_map = struct('x',[],'y',[],'level',[],'name',[]);
for i = 1:length(level)
	pos = (basinNumber_ == int64(level(i)));
	tmp = contour(x_,y_,double(pos'),[1]);
	contour_map(i).x = tmp(1,2:end);
	contour_map(i).y = tmp(2,2:end);
	contour_map(i).level = level(i);
	contour_map(i).name = sprintf('%d',level(i));
end
clf; % clear all existed contour maps.

% redraw contour map.
cm = jet;
cr = ceil(linspace(1,64,length(contour_map)));
for i = 1:length(contour_map),
	patch('xdata',contour_map(i).x,'ydata',contour_map(i).y,...
		'edgecolor',cm(cr(i),:),'facecolor','none');
end

fprintf('save %s to %s\n',nc,expname);
expwrite(contour_map,expname);
