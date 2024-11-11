function make_colorwheel(cmap,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate colorwheel files (circular colorscales) based on colormaps
%
% Input (required): 
%     > cmap < colormap to use in string format (e.g. matlab built-in cmaps 
%              like parula, jet, etc.) or external ones like the cyclic 
%              colormaps by F. Crameri which can be downloaded from: 
%              >>> http://www.fabiocrameri.ch/cycliccolourmaps.php <<<
%              To use the latter ones just place the corresponding mat-files
%              (vikO.mat, brocO.mat or corkO.mat) in a directory included in
%              your matlab path
%
% Input (optional): format as string, e.g. 'pdf', 'png' etc., default is
%                   set to pdf when no second parameter is given
%
% Examples: make_colorwheel('parula','pdf')
%           make_colorwheel('jet','png')
%           make_colorwheel('winter','pdf')
%           make_colorwheel('vikO','pdf') => F. Crameri's vikO.mat required
%           make_colorwheel('corkO','pdf') => F. Crameri's corkO.mat required
%           make_colorwheel('brocO','jpeg') => F. Crameri's brocO.mat required
%
% 2019-12-12 -MG-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

if nargin > 2
   error('Too many inputs!') 
end

% use pdf format as default
if nargin > 1
    format = varargin{1};
else
    format='pdf';
end

% check for colormaps
try  
    cmapuse=colormap(cmap);
catch
    cmapin=load([cmap '.mat']);
    cnames=fieldnames(cmapin);
    cmapuse=cmapin.(cnames{1});
end

% interpolate cmaps to go from 0 to 360 degrees (full circle)
if length(cmapuse) ~= 360
    disp('Interpolate...')
    cindex=linspace(1,length(cmapuse),360);
    c1=interp1(1:length(cmapuse),cmapuse(:,1),cindex);
    c2=interp1(1:length(cmapuse),cmapuse(:,2),cindex);
    c3=interp1(1:length(cmapuse),cmapuse(:,3),cindex);
    cmapm=[c1' c2' c3'];
else
    cmapm=cmapuse;
end

cmapa=colormap(cmapm);

% "mirror" colormap
cmapa=flipud(cmapa);

% "rotate" colormap, adjust for your needs 
theta=linspace(0.5*pi,2.5*pi,360);

for ii=1:360 % adjust linewidth for your needs 
    polarplot([theta(ii) theta(ii)],[1.15 2],'linewidth',10,'color',cmapa(ii,:));
    hold on
end

% remove axis, labels etc.
ax=gca; 
ax.ThetaGrid='off';
ax.RGrid='off';
ax.RTickLabel=[]; 
ax.ThetaTickLabel=[];

% required to set background later fully transparent (e.g. for eps format 
% when using the file in Generic Mapping Tools, GMT via the psimage command) 
set(gcf,'color','none');
set(gca,'color','none');

% save plot in selected format
print(['PLOT_colwheel_' cmap],['-d' format])

