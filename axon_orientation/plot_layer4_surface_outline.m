%%
ha1 = figure;
nr5a1 = nrrdread('NR5a1\red\top.nrrd');
h2 = imagesc(nr5a1(:,1:700));
colormap(hot);
caxis([min(nr5a1(:)),max(nr5a1(:))*1.2]);
% set(gca, 'YDir','reverse');
% axis image; axis off;
%%
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
atlas = double(atlas);
atlas1 =zeros(264,228);
for k1 = 1:264
    for k2 = 1:228
        if any(atlas(:,k1,k2))
            firstIndex = find(atlas(:,k1,k2)>0,1);
            atlas1(k1,k2) = atlas(firstIndex,k1,k2);
        end
    end
end
atlas1(:,115:end)=0;
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
root1 = '/997/';
areaName = {'MOp','MOs'};
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
maskPath{12} = '/997/8/567/688/695/315/453/322/353/'; % SSp-n
maskPath{13} = '/997/8/567/688/695/315/453/322/329/'; % bfd
maskPath{14} = '/997/8/567/688/695/315/453/322/337/'; % ll
maskPath{15} = '/997/8/567/688/695/315/453/322/345/'; %m
maskPath{16} = '/997/8/567/688/695/315/453/322/369/'; % ul
maskPath{17} = '/997/8/567/688/695/315/453/322/361/'; % tr
maskPath{18} = '/997/8/567/688/695/315/453/322/182305689/'; % un
maskPath{19} = '/997/8/567/688/695/315/669/385/'; % VISp

% ha1 = figure;
% 
% axx1 = subplot(1,1,1)
scale3 = 5;
hold on;
plotOutline(maskPath(1:3),st,atlas1,[],scale3);
plotOutline(maskPath(4),st,atlas1,[],scale3);
% plotOutline(maskPath(5),st,atlas1,[],scale3);
plotOutline(maskPath(6:11),st,atlas1,[],scale3);
for i = 6:19
    plotOutline(maskPath(i),st,atlas1,[],scale3);
end
axis image; axis off;
print(ha1, 'nr5a1_surface', '-dpdf', '-bestfit', '-painters');