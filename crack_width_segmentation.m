%  -----------------------------------------------------------------------

% This package contains programs written by Carrasco M, Araya-Letelier G, 
% Vel?zquez R,  % Visconti P., for the implementation of the Image-Based 
% Automated Width Measurement of  % Surface Cracking  described in the article:
%     
%    Carrasco M, Araya-Letelier G, Vel?zquez R, Visconti P. Image-Based 
%    Automated Width Measurement  of Surface Cracking. Sensors. 2021; 21(22):7534. 
%    https://doi.org/10.3390/s21227534 

%  The main function is "crack_width_segmentation.m".
% 
%  Author: Miguel Carrasco. (miguel.carrasco@uai.cl)
%  version: 1.0  01-08-2018
%  version: 1.1  18-10-2018. 
%  version: 1.2  13-09-2021. 



clc;
close all;
clear all;
warning off;

% Image folder
data =dir('data/*.JPG');
n = size(data,1);

figura_analisis = 12; % Chose id image

%=================== PARAMETERS =====================
dm           = 2001; % output dimension
block_points = 10;   % maximum distance per block
step         = 5;    % distance between points
dt           = 42;   % perimenter distance

pct_cut_off = 0.75;  %0.75
tol_distancia= 20;   % inner distance between points
hpts          = 30;  % marker distance

%===================== OTHERS ========================
paint_micro  = 0;    %show output image
paint_lines  = 1;    % show normal lines 
px_to_mm     = [dm-1 300];    % equivale a 2.5 mm con dm=1000;
boot_n       = 1e3;  %numero de valores de bootstrap
kstep        = 10;   

%%
% PHASE I) 
mat_H3 = fspecial('disk',30); 
mat_H1 = fspecial('gaussian',7 );
se     = strel('disk',3);
micro  = strel('square',2);
blurr  = fspecial('gaussian', 8);

points_circle = 32; %16
rad    = linspace(0,360,points_circle);


fixedPoints = [ 0    0
    dm   0
    dm   dm
    0    dm];

% for each image
for i=figura_analisis 
    
    filename = data(i).name
    I=imresize(imread(sprintf('data/%s',filename)),0.5);
    
    ch_R= I(:,:,1);
    ch_G= I(:,:,2);
    ch_B= I(:,:,3);
    
    h_fig= figure;
    h_im= imshow(rgb2gray(I)); hold on; axis on; grid on
    h=impoly;
    input_points= wait(h);
    close(h_fig);
    
    % Homomorphic Transformation
    ch_R_tmp= proyeccion(input_points, dm, ch_R);
    ch_G_tmp= proyeccion(input_points, dm, ch_G);
    ch_B_tmp= proyeccion(input_points, dm, ch_B);
    
    % original dimension
    ipts = [input_points; input_points(1,:)];
    for k=1:4
        dps(k)=  sqrt((ipts(k,1)-ipts(k+1,1))^2 + (ipts(k,2)-ipts(k+1,2))^2 );
    end
    scale = mean(dps)/1000;
    
    J(:,:,1)=ch_R_tmp;
    J(:,:,2)=ch_G_tmp;
    J(:,:,3)=ch_B_tmp;
    
    %%
    % STEP I. L*a*b 
    lab = rgb2lab(J);
    Ilab = imresize(lab, 1/4);
    
    %%
    % STEP II. Perona Malik filtering
    % USE Coherence Filter from Dirk-Jan Kroon (V5)
    JSI = CoherenceFilter(Ilab(:,:,1),struct('T',15,'rho',1,'Scheme','N'));
    JSI_lab = imresize(JSI, 4);
    
    figure, imshow(JSI,[]);
    
    sep_best= lab(:,:,1);
    JS= sep_best;
    JS = JS-min(JS(:));
    JS = JS./max(JS(:));
    data_mean = sort(JSI(:));
 
    x_points= 1:kstep:length(data_mean);
    y_points=interp1(1:length(data_mean), data_mean, 1:kstep:length(data_mean));
    
    area_full=[];
    for ii=1:length(y_points)-1
        area_full(ii)=  y_points(ii)* (length(y_points)-ii);
    end
    
    [~, pos] = max(area_full)
    umbral_full = y_points(pos);
    
    %%
    % STEP III: Binary Segmentation
    bwe = JSI<umbral_full;
    bwe = imresize(bwe,4);
    se = strel('disk',3);
    dst = bwdist(not(bwe));
    im_bwe = imtophat(dst, se);
    
    %%
    % STEP IV. Skeletonization
    bw_skel= bwmorph(im_bwe, 'skel',inf);
    [bwyy,bwxx] = find(bw_skel==1);
    M_sel_coords =  [bwyy,bwxx];
    
    %%
    % STEP V. Kmeans clustering. Iterative step
    
    sw=1;
    clusters=100;
    while(sw)
        [~, M_clus_coords] = kmeans(M_sel_coords,clusters);
        
        data= zeros(1,clusters);
        for pt=1:clusters
            
            point = repmat(M_clus_coords(pt,:), clusters,1);
            tmp =  sort(sqrt(sum((point-M_clus_coords).^2,2)));
            data(pt)= tmp(2);
        end
        distancia =  mode(data);
        
        if(distancia<tol_distancia )
            sw=0;
            fprintf('distance: %1.1f, clusters: %i\n',distancia, clusters);
        else
            clusters= clusters+100;
            fprintf('distance: %1.1f, clusters: %i\n',distancia, clusters);
        end
    end
    
    
    figure, imshow(JS,[]); hold on; axis on; grid on;
    plot(M_clus_coords(:,2), M_clus_coords(:,1), 'gx', 'markersize', 4); drawnow;
        
    
    sts= regionprops(bwe,'all');
    bw_out= zeros(size(bwe));
    
    %%
    % STEP VI. Curve Fitting
    pts=size(M_clus_coords,1);
    jj= ceil(M_clus_coords(:,2));
    ii= ceil(M_clus_coords(:,1));
    
    plot(jj,ii, 'rs','MarkerSize',6,'LineWidth',1);
    ancho_bin = zeros(pts,1);
    orientation = zeros(pts,1);
    
    % for each point detected previusly
    for t=1:pts
        ix= ii(t)-dt:ii(t)+dt-1;
        iy= jj(t)-dt:jj(t)+dt-1;
        try
            if (min(ix)>0 && min(iy)>0)
                section_bw = bwe(ix, iy);
                
                sts_tmp =regionprops(section_bw, 'MinorAxisLength','Orientation','Centroid', 'Area');
                if(size(sts_tmp,1)>1)
                    [~, id_max]= max([sts_tmp.Area]);
                    area_max= sts_tmp(id_max).Area;
                else
                    id_max=1;
                    area_max= sts_tmp(id_max).Area;
                end
                
                if (area_max>10)
                    
                    %%
                    % STEP VII. Curve fitting
                    angle = deg2rad(sts_tmp(id_max).Orientation)+0.0001;
                    m= tan(angle);
                    
                    %recta normal
                    mt= tan(angle-pi/2);
                    
                    lado_b = cos(angle-pi/2)*hpts;
                    xx = linspace(jj(t)-lado_b,jj(t)+lado_b,10);
                    xx = xx-mean(xx);
                    yy  = -m*xx+ii(t);
                    yyt  = -mt*xx+ii(t);
                    xx = xx+ jj(t);
                    
                    plot(xx, yyt, 'm-','LineWidth', 2); drawnow;
                    
                    %%
                    % STEP VIII. Profile width estimation
                    [cx,cy,prf]=improfile(JSI_lab, xx, yyt, 100, 'bilinear');
                    [cx,cy,prfJS]=improfile(JS, xx, yyt, 100, 'bilinear');
                    
                    prf_JS = (prfJS-min(prfJS));
                    prf_JS = prf_JS/max(prf_JS);
                    
                    [~, centerJS] = min(prf);
                    left = 1:centerJS;
                    right = centerJS+1:length(prf_JS);
                    
                    data_left=prf_JS(left);
                    clus_left=kmeans(prf_JS(left),2)-1;
                    pos_left= find (abs(diff(clus_left))==1);
                    distance_left = left(end)-pos_left(end);
                    
                    
                    data_right=prf_JS(right);
                    clus_right=kmeans(prf_JS(right),2)-1;
                    pos_right= find(abs(diff(clus_right))==1);
                    distance_right = pos_right(1);
                    
                   
                    % noise filtering
                    prf = prf(not(isnan(prf)));
                    class =kmeans(prf,2);
                    
                    value_min_cluss =min( [mean(prf(class==1)), mean(prf(class==2))]);
                    value_max_cluss =max( [mean(prf(class==1)), mean(prf(class==2))]); 
                    diff_mean = value_min_cluss/value_max_cluss;
                    
                    [value_min bit] = min([mean(prf(class==1)) mean(prf(class==2))]);
                   
                   
                    if diff_mean<pct_cut_off
                        
                        p1=[cx(centerJS-distance_left) cy(centerJS-distance_left)];
                        p2=[cx(centerJS+distance_right) cy(centerJS+distance_right)];
                        plot(p1(1), p1(2), 'c+', 'MarkerSize',18);
                        plot(p2(1), p2(2), 'c+','MarkerSize',18);
                        ancho_segment = norm(p1-p2);
                        
                        tmp_ancho = sum(class ==bit);
                        ancho_bin(t, :)= ancho_segment;
                        
                        
                        orientation(t)= sts_tmp(id_max).Orientation;
                        
                        if paint_lines
                            lado_b = cos(angle-pi/2)*hpts;
                            xx = linspace(jj(t)-lado_b,jj(t)+lado_b,10);
                            xx = xx-mean(xx);
                            yy  = -mt*xx+ii(t);
                            xx = xx+ jj(t);
                          
                        end
                      
                        
                    else
                        ancho_bin(t, :)=0;
                    end
                    
                end
            end
        catch
            %fprintf('\n out of limits');
        end
        
    end
    
    
    
    back = ancho_bin;
    data_bin= ancho_bin(back>1);
    
    L=data_bin(:);
    [x, bootsam]=bootstrp(boot_n,@mean,L);
    csv_filename = sprintf('%s%s',strtok(filename,'.JPG'), '.csv')
    dlmwrite(csv_filename, L)
    
    prom=mean(data_bin)*px_to_mm(2)/px_to_mm(1)*scale;
    fprintf('\nMean crack width: %3.3f mm', prom);
    
    figure, hist(sort(abs(orientation(abs(orientation)>0))))
    csv_orientation = sprintf('%s%s',strtok(filename,'.JPG'), '_orientation.csv')
    dlmwrite(csv_orientation, orientation)
    
    title('Main orientation')
    xlabel('Angle')
    
end


