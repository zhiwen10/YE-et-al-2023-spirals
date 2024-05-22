% chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_hStripe_botRow0_ref1.mat';
function chanMap1 = chanMapReorder(Map)
%% reorder probes to 4 shanks
[chanMap, xc, yc, kcoords, NchanTOTdefault] = loadChanMap(Map); % function to load channel map file
%%
if contains(Map,'NPtype24_hStripe')
    xc1 = xc+282;
    yc1 = (yc+15)/15+1;
    mapReal = [xc1,yc1];
    %% find chanMap order to reorder ephys data to 2d map
    map1(1:48,1) = zeros(48,1); % x value
    map1(1:48,2) = 1:48; % y value
    map2 = [];
    for i = 1:8
        mapTemp = map1;
        if mod(i,2) ==0
            mapTemp(:,1) = mapTemp(:,1)+32+250*(i/2-1);
        else
            mapTemp(:,1) = mapTemp(:,1)+250*(i-1)/2;
        end
        map2 = [map2;mapTemp];
    end
    [a,b] = ismember(map2,mapReal,'rows');
    chanMap1 = reshape(b,[48,8]);
elseif contains(Map,'NPtype24_doubleLengthStripe')
    for k  = 1:numel(yc)
        if yc(k)<= 705
            yc1(k) = yc(k)/15+1;
        else yc1(k) = (yc(k)-7.5)/15+1;
        end
    end
    for j = 1:numel(xc)
        if mod(xc(j)-32,250)==0
            xc1(j) = xc(j)-32;
        else xc1(j) = xc(j);
        end
    end     
    mapReal = [xc1',yc1'];
    %% find chanMap order to reorder ephys data to 2d map
    map1(1:96,1) = zeros(96,1); % x value
    map1(1:96,2) = 1:96; % y value
    map2 = [];
    for i = 1:4
        mapTemp = map1;
        mapTemp(:,1) = mapTemp(:,1)+250*(i-1);
        map2 = [map2;mapTemp];
    end
    [a,b] = ismember(map2,mapReal,'rows');
    chanMap1 = reshape(b,[96,4]);
elseif contains(Map,'NPtype24_quadrupleLengthStripe')
    for k  = 1:numel(yc)
        if yc(k)<= 1400
            yc1(k) = yc(k)/15+1;
        else yc1(k) = (yc(k)-15)/15+1;
        end
    end
    for j = 1:numel(xc)
        if mod(xc(j)+32,250)==0
            xc1(j) = xc(j)+32;
        else xc1(j) = xc(j);
        end
    end     
    yc1 = yc1/2+1;
    mapReal = [xc1',yc1'];
    %% find chanMap order to reorder ephys data to 2d map
    map1(1:96,1) = zeros(96,1); % x value
    map1(1:96,2) = 1:96; % y value
    map2 = [];
    for i = 1:4
        mapTemp = map1;
        mapTemp(:,1) = mapTemp(:,1)+250*(i-2);
        map2 = [map2;mapTemp];
    end
    [a,b] = ismember(map2,mapReal,'rows');
    chanMap1 = reshape(b,[96,4]);
elseif contains(Map,'neuropixPhase3B1')
    for k  = 1:numel(yc)
        yc1(k) = yc(k)/20;
        if mod(yc1(k),2) == 0
            yc1(k) = yc1(k)-1;
        end
    end
    for j = 1:numel(xc)
        xc1(j) = (xc(j)-11)/16+1;
    end   
    yc1 = (yc1+1)/2;
    mapReal = [xc1',yc1'];
    %% find chanMap order to reorder ephys data to 2d map
    map1(1:96,1) = zeros(96,1); % x value
    map1(1:96,2) = 1:96; % y value
    map2 = [];
    for i = 1:4
        mapTemp = map1;
        mapTemp(:,1) = mapTemp(:,1)+i;
        map2 = [map2;mapTemp];
    end
    [a,b] = ismember(map2,mapReal,'rows');
    chanMap1 = reshape(b,[96,4]);
    chanMap1(48,2) = 384;
elseif contains(Map,'NPtype21_botRow0')
    for k  = 1:numel(yc)
        yc1(k) = (yc(k)+30)/15;
    end
    for j = 1:numel(xc)
        xc1(j) = (xc(j)+282)/32+1;
    end   
    mapReal = [xc1',yc1'];
    %% find chanMap order to reorder ephys data to 2d map
    chanMap1 = zeros(192,2); % x value
    for i =1:384
        chanMap1(mapReal(i,2),mapReal(i,1)) = i;
    end
end
    