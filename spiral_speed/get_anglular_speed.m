function [angle_offset_all,distance_offset_all] = get_anglular_speed(filteredSpirals3,tracePhase_padded,scale,halfpadding)
pixSize = 3.45/1000/0.6*3;  % mm / pix
for iframe = 1:size(filteredSpirals3,1)
    % iframe = 25;
    frame = filteredSpirals3(iframe,5);
    spiral_radius = (filteredSpirals3(iframe,3))/scale;
    center = round([filteredSpirals3(iframe,2)/scale,filteredSpirals3(iframe,1)/scale]);
    % center = round([(filteredSpirals3(kk,2)-150)/2,(filteredSpirals3(kk,1)-210)/2]);
    center = center+halfpadding;
    frame1 = squeeze(tracePhase_padded(:,:,frame));
    frame2 = squeeze(tracePhase_padded(:,:,frame+1));
    phase_circle1 = []; phase_circle2 = [];
    phase_circle1a = []; phase_circle2a = [];
    phase_circle1b = []; phase_circle2b = [];
    angle_offset = []; distance_offset = [];
    angle = pi/6:pi/6:2*pi; % sample evenly around a circle
    count2 = 1;
    for radius = 1:1:spiral_radius
        for k = 1:numel(angle)
            angle1 = angle(k);
            x(k) = center(2)+round(radius*cos(angle1));
            y(k) = center(1)+round(radius*sin(angle1));
            phase_circle1(k,count2) = frame1(y(k),x(k));
            phase_circle2(k,count2) = frame2(y(k),x(k));
        end
        count2 = count2+1;
    end
    angle_offset = [];
    % phase angle across different radius
    for kk = 1:size(phase_circle1,2)
        phase_circle1a = phase_circle1(:,kk);
        phase_circle2a = phase_circle2(:,kk);
        phase_circle1b(:,kk) = wrapAngle(phase_circle1a);
        phase_circle2b(:,kk) = wrapAngle(phase_circle2a);
    end
    phase_circle1b(abs(phase_circle1b(:))<0.0001) = nan;
    phase_circle1b(abs(phase_circle1b(:)-2*pi)<0.0001) = nan;
    phase_circle2b(abs(phase_circle2b(:))<0.0001) = nan;
    phase_circle2b(abs(phase_circle2b(:)-2*pi)<0.0001) = nan;
    angle_offset = phase_circle2b-phase_circle1b;
    angle_offset = wrapTo2Pi(angle_offset);
    angle_offset(angle_offset(:)>pi) = angle_offset(angle_offset(:)>pi)-pi;
    mean_angle_offset = mean(angle_offset,1);   
    radius_all = 1:1:spiral_radius;
    for kk = 1:numel(radius_all)
        distance_offset(:,kk) = (angle_offset(:,kk)./(2*pi)*(2*pi*radius_all(kk)*pixSize*scale))*35;
    end
    angle_offset_all{iframe} = angle_offset;
    distance_offset_all{iframe} = distance_offset;
end