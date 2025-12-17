function [Z,Iy] = RSG(Py,Nf,N,clean_method,threshold,epsInterval,eps_def,minptsInterval,trimming,endpoints,subsel,debug)
%% For more details, see https://doi.org/10.3390/math13101626 and https://github.com/aerfangel/RSG
%%
%% The Reference Set Generator (RSG) computes a reference (Z) of size N 
%% from the starting set Py. The specific procedure is different for k=2 
%% and k>2 objectives, however for both cases a filling is performed, 
%% obtaining a set Iy with Nf points. You can select to include the 
%% endpoints in the final set Z by setting endpoints=true. The variable 
%% subsel can be set to 'means', 'medoids', or 'spectral' to select N 
%% points from Iy using either k-means, k-medoids or spectral clustering. 
%% If you want to only compute the filled set Iy, without computing Z, set 
%% trimming = false. The variables epsInterval, eps_def, and
%% minptsInteval, control the component detection. For the component
%% detection, a gridsearch will be performed. If Py consists of one
%% component, we recommend setting epsInterval =
%% [max(max(Py))-min(min(Py)),max(max(Py))-min(min(Py))], eps_def = 1 and
%% minptsInterval = [3,3]. If Py consists of several components (or if 
%% nothing is known about Py), we recommend setting minptsInterval=[2,3],
%% epsInterval=[0.1*d,0.15*d] and eps_def=0.1*d for bi-objective problems,
%% and minptsInterval=[3,4], epsInterval=[0.19*d,0.23*d], and
%% eps_def=0.1*d otherwise, where d is the average pairwise distance 
%% between all the points. Finally, for problems with k>2 objectives, the 
%% variable clean_method and its threshold are parameters for the triangle
%% cleaning. We recommend setting clean_method='long', and to set the
%% values for threshold we recommend the user to have a look at the manual 
%% at: https://github.com/aerfangel/RSG. For biobjective problems, set both 
%% clean_method and threshold as empty arrays.
%%
%% OUTPUT:
%%     Z:  RSG reference set
%%     Iy: filled set
%% INPUT:
%%     Py:             starting set
%%     N:              number of desired points for Z
%%     Nf:             number of desired points for the filling
%%     clean_method:   (for problems with k>2 objectives) it can be 'long', 'cond', 'area', and 'off'. 'long' (recommended), cleans the triangulation using the largest side of the triangles, 'cond' uses the conditional number of the matrix containing the vertices, 'area' uses the area and 'off' skips the cleaning step. We suggest reading the github manual to find the best values.
%%     threshold:      (for problems with k>2 objectives) the threshold used to clean the triangulation. We suggest reading the github manual to find the best values.
%%     epsInterval:    an array [a,b] where a and b are the desired initial and final radius for DBSCAN
%%     eps_def:        the stepsize for epsInterval. DBSCAN will try radius' a, a+eps_def, a+2*eps_def, ..., b
%%     minptsInterval: an array [a,b] where a and b are the desired initial and final minimum number of points for DBSCAN
%%     trimming:       a value of 1 will reduce the filled set Iy to Z, a value of 0 will only compute the filled set
%%     endpoints:      a value of 1 will include the endpoints in Z, a value of 0 will not
%%     subsel:         can be 'means', 'medoids' or 'spectral'. This will select the reduction method to go from Iy to Z

%%%-----------------------V1.0 04/JUNE/2025-----------------------------%%%

if nargin <= 11
    debug = false; % Default value
end


name = [];
basename = [];
show_plots = false;
save_plots = false;
degeneration = false;
comp_points = 'length';
if isempty(clean_method)
    clean_method = 'long'; % Default value for clean_method
end
if isempty(threshold)
    threshold = 3; % Default value for threshold
end

k = size(Py,2);
cluster_iters = 5000;

%----Interpolate--%
Iy = interpolation( Py,Nf,clean_method,threshold,name,basename,epsInterval,eps_def,minptsInterval,comp_points,show_plots,save_plots,degeneration,debug);
if show_plots
    figure
    if size(Iy,2) == 2
        scatter(Iy(:,1),Iy(:,2),'.')
    else
        scatter3(Iy(:,1),Iy(:,2),Iy(:,3),'.')
        view(135,25)
        zlabel('f_3')
    end
    fig = gcf; fig.PaperUnits = 'centimeters'; fig.PaperSize = [9 9]; fig.PaperPosition = [0 0 9 9]; xticklabels(strrep(xticklabels,'-','–')); yticklabels(strrep(yticklabels,'-','–'));
    xlabel('f_1')
    ylabel('f_2')
    % title('Filling')
    if save_plots
        saveas(gcf, [ pwd '/' basename '/Plots/' name '_Filling.fig' ] );
        saveas(gcf, [ pwd '/' basename '/Plots/' name '_Filling.pdf' ] );
        saveas(gcf, [ pwd '/' basename '/Plots/' name '_Filling.png' ] );
    end
end
%-----Trimming----%
if trimming
    if ~endpoints
        if strcmp(subsel,'spectral')
            idx = spectralcluster(Iy,N);
            Z = Iy(idx,:);
        elseif strcmp(subsel,'medoids')
            [~,Z,~,~,~] = kmedoids(Iy,N,'Options',statset('MaxIter',cluster_iters));
        else
            [~,Z] = kmeans(Iy,N,'Options',statset('MaxIter',cluster_iters));
        end
    else
        k = size(Iy,2);
        [~,I] = max(Iy);
        if length(unique(I)) ~= k
            %SOLUCION MIN
            [~,I] = min(Iy);
        end
        if strcmp(subsel,'spectral')
            idx = spectralcluster(Iy,N-k);
            Z = Iy(idx,:);
        elseif strcmp(subsel,'medoids')
            [~,Z,~,~,~] = kmedoids(Iy,N-k,'Options',statset('MaxIter',cluster_iters));
        else
            [~,Z] = kmeans(Iy,N-k,'Options',statset('MaxIter',cluster_iters));
        end
        Z(end+1:end+k,:) = Iy(I,1:k);
%         Z(end+1:end+k,:) = Iy(I,:);
    end
else
    Z = [];
end

if trimming
    if show_plots
        figure
    %     close all
        if k==2
            scatter(Z(:,1),Z(:,2),'d','filled')
            hold on
            scatter(Py(:,1),Py(:,2),'.')
        else
            scatter3(Z(:,1),Z(:,2),Z(:,3),'filled')
            view(135,25)
            hold on
            scatter3(Py(:,1),Py(:,2),Py(:,3),'.')
            zlabel('f_3')
        end
        legend('Z','Py')
        fig = gcf; fig.PaperUnits = 'centimeters'; fig.PaperSize = [9 9]; fig.PaperPosition = [0 0 9 9]; xticklabels(strrep(xticklabels,'-','–')); yticklabels(strrep(yticklabels,'-','–'));
        % title('Final Result (Z) vs Start (Py)')
        xlabel('f_1')
        ylabel('f_2')
        if save_plots
            saveas(gcf, [ pwd '/' basename '/Plots/' name '_RSG_N=' int2str(N) '.fig' ] );
            saveas(gcf, [ pwd '/' basename '/Plots/' name '_RSG_N=' int2str(N) '.pdf' ] );
            saveas(gcf, [ pwd '/' basename '/Plots/' name '_RSG_N=' int2str(N) '.png' ] );
        end
    end
end

end

%%%--------------------------INTERPOLATION------------------------------%%%
function [Fy] = interpolation( Py, Nf, clean_method, threshold, name, basename, epsInterval, eps_def, minptsInterval, comp_points, show_plots, save_plots, degeneration, debug )
    num_obj = size(Py,2);

%--------------------------------------------------------%
%--------------------2D FILLING--------------------------%
%--------------------------------------------------------%
    if num_obj == 2
        %2D Filling
%         epsInterval = [0.25,0.25];
%         eps_def = 0.01;
%         minptsInterval = [2,3];
                
        [~,~,C,num_clusters,~,~] = cluster_gridsearch_(Py,epsInterval,eps_def,minptsInterval);

        if show_plots
            figure
            for i=1:num_clusters
                scatter(C{i}(:,1),C{i}(:,2),'filled');
                hold on
            end
            scatter(Py(:,1),Py(:,2),'.')
            fig = gcf; fig.PaperUnits = 'centimeters'; fig.PaperSize = [9 9]; fig.PaperPosition = [0 0 9 9]; xticklabels(strrep(xticklabels,'-','–')); yticklabels(strrep(yticklabels,'-','–'));
            xlabel('f_1')
            ylabel('f_2')
            % title('Component Detection')
            if save_plots
                saveas(gcf, [ pwd '/' basename '/Plots/' name '_Components.fig' ] );
                saveas(gcf, [ pwd '/' basename '/Plots/' name '_Components.pdf' ] );
                saveas(gcf, [ pwd '/' basename '/Plots/' name '_Components.png' ] );
            end
        end

        %------------------------------------------------INTERPOLATION
        I = {};
        if strcmp(comp_points,'points')
            %Modificación: Proporcional al número de puntos:
            TotalPointsCluster = 0;
            for i=1:num_clusters
                TotalPointsCluster = TotalPointsCluster + size(C{i},1);
            end
            %
        elseif strcmp(comp_points,'length')
            comp_length = zeros(1,num_clusters);
            for comp=1:num_clusters
                X = C{comp};
                [n,~] = size(X);
                Points = sortrows(X,1);
                distances = zeros(1,n-1);
                for i=1:n-1
                    distances(i) = norm(Points(i+1,:)-Points(i,:));
                end
                comp_length(comp) = sum(distances);%SUM OF THE LENGTH SEGMENTS
            end
            total_length = sum(comp_length); 
        end
        
        for i=1:num_clusters
            
            close all

            %Modificación: Proporcional al número de puntos:
            if strcmp(comp_points,'points')
                I{i} = interpolate_tst(C{i},round(Nf*size(C{i},1)/TotalPointsCluster));
            elseif strcmp(comp_points,'components')
                I{i} = interpolate_tst(C{i},round(Nf/num_clusters));
            elseif strcmp(comp_points,'length')
                I{i} = interpolate_tst(C{i},round(Nf*comp_length(i)/total_length));
            end
        end
        comp_idx_tot = [];
        Fy = [];
        for i=1:length(I)
            for j= 1:length(I{i})
                comp_idx_tot(end+1) = i; 
                Fy(end+1,:) = I{i}(j,:);
            end
        end
        
%--------------------------------------------------------%
%--------------------k>3 FILLING-------------------------%
%--------------------------------------------------------%
    elseif num_obj >= 3
        %--------------------COMPONENT DETECTION--------------------
%         epsInterval = [0.25,0.30]; %TESTS!!! MODIFIED VALUE
%         eps_def = 0.05; 
%         minptsInterval = [3,4];

        [~,~,C,num_clusters,~,~] = cluster_gridsearch_(Py,epsInterval,eps_def,minptsInterval);

        if show_plots
            figure
            for i=1:num_clusters
                scatter3(C{i}(:,1),C{i}(:,2),C{i}(:,3),'filled');
                hold on
            end
            scatter3(Py(:,1),Py(:,2),Py(:,3),'.')
            view(135,25)
            fig = gcf; fig.PaperUnits = 'centimeters'; fig.PaperSize = [9 9]; fig.PaperPosition = [0 0 9 9]; xticklabels(strrep(xticklabels,'-','–')); yticklabels(strrep(yticklabels,'-','–'));
            xlabel('f_1')
            ylabel('f_2')
            zlabel('f_3')
            % title('Component Detection')
            if save_plots
                saveas(gcf, [ pwd '/' basename '/Plots/' name '_Components.fig' ] );
                saveas(gcf, [ pwd '/' basename '/Plots/' name '_Components.pdf' ] );
                saveas(gcf, [ pwd '/' basename '/Plots/' name '_Components.png' ] );
            end
        end

        %-----------------------Random INTERPOLATION--------------------
        %Modificación: Proporcional al número de puntos:
        TotalPointsCluster = 0;
        for i=1:num_clusters
            TotalPointsCluster = TotalPointsCluster + size(C{i},1);
        end
%             I{i} = interpolate_tst(C{i},round(Nf*size(C{i},1)/TotalPointsCluster));
        
        %RANDOM INTERPOLATION FOR EACH CLUSTER
        I_rand = cell(1,num_clusters);
        for i = 1:num_clusters
            
%             close all

            %Modificación: Proporcional al número de puntos: (comentar para
            %todos mismo num interpol pts
            interpol_points = round(Nf*size(C{i},1)/TotalPointsCluster);
            
            if size(C{i},1) <= 3
                I_rand{i}=C{i};
            elseif size(unique(C{i}(:,1)),1) <= 3
                I_rand{i}=C{i};
            elseif degeneration
                %DEGENERATION CASE
                
            else
                comp_name = [name '_comp' int2str(i)];
        	    DT = surf_triangulation( C{i}, clean_method, threshold, comp_name, basename, show_plots, save_plots, debug );

                total_volume = 0;
                all_volumes = [];
                for j=1:size(DT,1)
                    %-----Compute Volume-------------%
                    Volume = zeros(num_obj,num_obj-1);
                    for l=2:num_obj
                        Volume(:,l-1) = C{i}(DT(j,l),:) - C{i}(DT(j,1),:);
                    end
                    Volume = det(Volume'*Volume)/factorial(num_obj);
                    all_volumes(end+1) = Volume;
                    total_volume = total_volume + Volume;
                    %-------------------------------------%
                end
%                 I_rand{i} = interpolate_triangle_random(DT.Points(K(1,1),:),DT.Points(K(1,2),:),DT.Points(K(1,3),:),total_area,interpol_points);
                I_rand{i} = interpolate_simplex_random(C{i}(DT(1,:),:),all_volumes(1),total_volume,interpol_points);
                for j=2:size(DT,1)
%                     points = interpolate_triangle_random(DT.Points(K(j,1),:),DT.Points(K(j,2),:),DT.Points(K(j,3),:),total_area,interpol_points);
                    points = interpolate_simplex_random(C{i}(DT(j,:),:),all_volumes(j),total_volume,interpol_points);
                    I_rand{i} = [I_rand{i};points];
                end
            end
        end
        
        %FOR FUNCTION OUTPUT
        I = I_rand;

        Fy = [];
        comp_idx_tot = [];
        for i=1:size(I_rand,2)
            for j= 1:size(I_rand{i},1)
                comp_idx_tot(end+1) = i;
                Fy(end+1,:) = I_rand{i}(j,:);
            end
        end
    end
end
%-----------------------INTERPOLATION SUBFUNCTIONS------------------------%
%----2D Interpolation-----%
function I = interpolate_tst(X,Nf)

    [n,p] = size(X);
    I = zeros(Nf-1,p);
    Points = sortrows(X,1);
    endpoints = [Points(1,:); Points(end,:)];
    distances = zeros(1,n-1);
    points_p_segment = zeros(1,n-1);
    for i=1:n-1
        distances(i) = norm(Points(i+1,:)-Points(i,:));
    end
    % maxdistance = norm(Points(1,:),Points(end,:));
    total_length = sum(distances); %SUM OF THE LENGTH SEGMENTS
    cum_dist = cumsum(distances); %CumSum to see where to put points
    
%     if ~increase_interval
        interpol_length = total_length/(Nf-1);
        dist_left = zeros(1,n);
        for i=1:n-1
            ratio = (distances(i)+dist_left(i))/interpol_length;
            points_p_segment(i) = floor(ratio);
            dist_left(i+1) = (ratio-floor(ratio))*interpol_length;
        end
        count=1;
        for seg=1:n-1
            if points_p_segment(seg)>0
                direction = Points(seg+1,:)-Points(seg,:);
                direction = direction/norm(direction);
                I(count,:) = Points(seg,:)+(interpol_length-dist_left(seg))*direction;
                count=count+1;
                for i=2:points_p_segment(seg)
                    I(count,:) = I(count-1,:)+interpol_length*direction;
                    count=count+1;
                end
            end
        end
%     end 

    I = [endpoints(1,:); I];
    I(end,:) = endpoints(end,:);
    if count<Nf-2
        disp('CUIDADO, MENOS PUNTOS ENCONTRADOS, SET DEGENERADO')
        for i=count:Nf-2
            I(i,:) = endpoints(randi(2),:);
        end
    end
end
%-------------------------%
%----k>3 Interpolation----%
function points = interpolate_simplex_random(P,P_volume,total_volume,total_points)
num_obj = size(P,2);

tri_points = ceil(total_points*P_volume/total_volume); %points for triangle

points = zeros(tri_points,num_obj);
alpha = rand(tri_points,num_obj);
for i=1:tri_points
    alpha(i,:) = alpha(i,:)/norm(alpha(i,:));
    alpha(i,:) = alpha(i,:).^2;
    points(i,:) = sum(P.*alpha(i,:)');
%     points(i,:) = alpha(1,1)*P(1,:)+alpha(1,2)*P(2,:)+alpha(1,3)*P(3,:)+alpha(1,4)*P(4,:);
end

end
%-------------------------%
%--Component Detection----%
function [eps_fin,minpts_fin,C_fin,num_clust_fin,idx_fin,avrg_dist] = cluster_gridsearch_(M,epsInterval,eps_def,minptsInterval)
    num_points = size(M,1);
    Distance = zeros(num_points);
    
    for i=1:num_points
        for j=i:num_points
            Distance(i,j) = norm( M(i,:)-M(j,:) );
            Distance(j,i) = Distance(i,j);
        end
    end
    avrg_dist = sum(sum(Distance)/size(Distance,1))/size(Distance,1);
%     eps = (0.19:0.01:0.23)*avrg_dist; 
    eps = (epsInterval(1):eps_def:epsInterval(2))*avrg_dist;
    num_eps = size(eps,2);
    minpts = minptsInterval(1):minptsInterval(2);
    num_minpts = size(minpts,2);

    % % % % %DEBUG!!! COMENT LATER!! TO TEST COMPONENT DETECTION
    % % % % i=2;
    % % % % j=1;
    % % % % [idx,~] = dbscan(M,eps(i),minpts(j));
    % % % % num_clusters = max(idx);
    % % % % for k=1:num_clusters
    % % % %     scatter3(M(idx==k,1),M(idx==k,2),M(idx==k,3),'.')
    % % % %     hold on
    % % % % end
    % % % % scatter3(M(idx==-1,1),M(idx==-1,2),M(idx==-1,3),'x')


    WLC_min = inf;
%     WLC_max = 0;
    C_fin = {};
    eps_fin = 0;
    minpts_fin = 0;
    num_clust_fin = 0;
    idx_fin = [];
    for i=1:num_eps
%         disp(['Current Epsilon:' int2str(i) '/' int2str(num_eps) ])
        for j=1:num_minpts
            [idx,~] = dbscan(M,eps(i),minpts(j));
            num_clusters = max(idx);
            C = {};
            for k=1:num_clusters
                C{k} = M(find(idx==k),:);
            end
            WLC = WeakestLinkCluster(C);
%             disp([['WLC= ' num2str(WLC) ', eps=' int2str(i) ', minpts=' int2str(minpts(j)) ', ' int2str(num_clusters) ' clusters']])
            if WLC < WLC_min
%             if WLC > WLC_max
                WLC_min = WLC;
%                 WLC_max = WLC;
%                 disp(['NUEVO WLC= ' num2str(WLC) ', eps=' int2str(i) ', minpts=' int2str(minpts(j)) ', ' int2str(num_clusters) ' clusters'])
                eps_fin = eps(i);
                C_fin = C;
                minpts_fin = minpts(j);
                num_clust_fin = num_clusters;
                idx_fin = idx;
            end
        end
    end
end
%-------------------------%
function WLC = WeakestLinkCluster(C)

    num_clust = size(C,2);
%     inter_clust_WLP = zeros(1,num_clust);

    %---This part computes the maximum link intra cluster-----%
    C_sizes = zeros(1,num_clust);
    max_clust = zeros(1,num_clust); %MAXIMUM LINK INTRA CLUSTER
    for i=1:num_clust
        C_sizes(i) = size(C{i},1);
    end
    for i=1:num_clust %compute value over all clusters
%         if C
        for j=1:C_sizes(i)-1 %find maximum dist in cluster
            dist = norm(C{i}(j,:)-C{i}(j+1,:));
            if dist > max_clust(i)
                max_clust(i) = dist;
            end
        end
    end
    %---------------------------------------------------------%

    %----COMPUTING INTRA CLUSTER WLP-----------%
    intra_cluster_WLP = 0;
    for i=1:num_clust-1
        for j=1:C_sizes(i)-1
            for k=j+1:C_sizes(i)
                temp = WeakestLinkPoints(max_clust,C{i}(j,:),C{i}(k,:),C);
                if temp > intra_cluster_WLP
                    intra_cluster_WLP = temp;
                end
            end
        end
    end
    %------------------------------------------%

    %----SHORTEST BETWEEN CLUSTER DISTANCE------%
    inter_clust_WLP = inf;
    for i=1:num_clust-1
        for j=i+1:num_clust
            temp = min_dist_2_clust(C{i},C{j});
            if temp < inter_clust_WLP
                inter_clust_WLP = temp;
            end
        end
    end
    %------------------------------------------%

    WLC = intra_cluster_WLP/inter_clust_WLP;


end
%-------------------------%
function WLP = WeakestLinkPoints(max_clust,x,y,C)
    num_clust = size(C,2);
    WLP_vec = zeros(1,num_clust);
    dist_vec = zeros(1,3);
    for i=1:num_clust
        dist_vec(1) = norm(C{i}(1,:)-x);
        dist_vec(2) = norm(C{i}(end,:)-y);
        dist_vec(3) = max_clust(i);
        WLP_vec(i) = max(dist_vec);
    end
    WLP = min(WLP_vec);
end
%-------------------------%
function min = min_dist_2_clust(X,Y)
    min = inf;
    X_size = size(X,1);
    Y_size = size(Y,1);
    for i=1:X_size
        for j=1:Y_size
            temp = norm(X(i,:)-Y(j,:));
            if temp < min
                min = temp;
            end
        end
    end
end
%-------------------------%
%%%---------------------------------------------------------------------%%%
%%%--------------------CLEAN TRIANGULATION------------------------------%%%
function [ DT_clean ] = rmtriangle( DT, Py, fun_param, threshold, comp_name, basename, show_plots, save_plots,debug )
    
    num_obj = size(Py,2);
    num_bins = 40;

	TA = [];
	for i = 1:size(DT,1)
		TA(i) = fun_param( Py(DT(i,:),1:end) );
	end
    
    if debug
        figure
	    hist = histogram( TA, num_bins );
	    title( 'Histogram');
    end

    if strcmp(threshold,'zero')
        first_zero = find(hist.Values == 0);
        if isempty(first_zero)
            T = inf;
        else
            T = hist.BinEdges(first_zero(1));
        end
    elseif strcmp(threshold,'decreasing')
        first_decreasing = find(diff(hist.Values)<0);
        if isempty(first_decreasing)
            T = inf;
        else
        T = hist.BinEdges(first_decreasing(1)+1);
        end
    else
%         T = mean(TA)*threshold;
        T = threshold;
    end

	
	TA2 = [];
	DT_clean = [];
	for i = 1:size(TA,2)
		% Accept
		if TA(i) <= T
			TA2(end+1) = TA(i);
			DT_clean(end+1,:) = DT(i,:);
		end
	end

    if debug
        figure
	    hist = histogram( TA2, num_bins );
	    title('Histogram Cleaned');
    end
end
%%%---------------------------------------------------------------------%%%
%%%-----------------------PROJECTION VECTOR-----------------------------%%%
function eta = shift_vect(Py,show_plots,save_plots,comp_name,basename)
    
    k = size(Py,2);
    
    if size(Py,1) <k
        eta = ones(1,k);
    else
        [~,I] = max(Py);
        if length(unique(I)) ~= k
            %SOLUCION CLUSTER
%             [~,~,~,~,I] = kmedoids(Py,k,'Options',statset('MaxIter',10));%MODIFICACION
            %SOLUCION RANDOM
%             for i=1:k
%                 index = find(Py(:,i)==Py(I(i),i));
%                 I(i) = max(index); %max
%         %         I(i) = index(randi(length(index))); %rand
%             end
            %SOLUCION MIN
            [~,I] = min(Py);
            if length(unique(I)) ~= k
                %SOLUCION CLUSTER
                [~,~,~,~,I] = kmedoids(Py,k,'Options',statset('MaxIter',10));%MODIFICACION
            end
        end
        
        if show_plots
            figure
            scatter3(Py(:,1),Py(:,2),Py(:,3),'.')
            hold on
            scatter3(Py(I,1),Py(I,2),Py(I,3),'filled')
            title('ENDPOINTS DETECTED')
            fig = gcf; fig.PaperUnits = 'centimeters'; fig.PaperSize = [9 9]; fig.PaperPosition = [0 0 9 9]; xticklabels(strrep(xticklabels,'-','–')); yticklabels(strrep(yticklabels,'-','–'));
            if save_plots
                saveas(gcf, [ pwd '/' basename '/Plots/' comp_name '_endpoints.fig' ] );
                saveas(gcf, [ pwd '/' basename '/Plots/' comp_name '_endpoints.pdf' ] );
                saveas(gcf, [ pwd '/' basename '/Plots/' comp_name '_endpoints.png' ] );
            end
        end

        %         [~,~,~,~,I] = kmedoids(Py,k);%MODIFICACION
        M = zeros(k,k-1);
        for i=1:k-1
            M(:,i) = Py(I(i+1),:) - Py(I(1),:);
        end
        
        [Q,R] = qr(M);
        
        qk = Q(:,end);
        % qk = Q(1,:);
        
        eta = -1*sign(qk(1))*qk;
        % eta = qk;
        
        eta = eta/norm(eta);
    end
end
%%%---------------------------------------------------------------------%%%
%%%-----------------------FILLING FOR k>3-------------------------------%%%
function DT = surf_triangulation(Py, clean_method, threshold, comp_name, basename, show_plots, save_plots, debug )
num_obj = size(Py,2);

%--------------------------------------------------------%
%--------------------k>3 FILLING-------------------------%
%--------------------------------------------------------%
        %CHANGE OF COORDINATES
    normal = shift_vect(Py,show_plots,save_plots,comp_name, basename);
    % normal = [1,1,1]; % ONLY FOR CONV3.4
%     backup = normal; %only for WFG1
%     normal(1) = backup(2);%only for WFG1
%     normal(2) = backup(1);%only for WFG1
%     normal(3) = 0;%only for WFG1
    
    
    if show_plots
        %Plot THE SHIFT DIRECTION
        figure
        s1 = scatter3(Py(:,1),Py(:,2),Py(:,3),'k','.','DisplayName', 'Py');
        hold on
        s2 = scatter3(Py(:,1)+normal(1),Py(:,2)+normal(2),Py(:,3)+normal(3),'y','.','DisplayName', 'Py + \eta');
        for i=1:size(Py,1)
            plot3([Py(i,1),Py(i,1)+normal(1)],[Py(i,2),Py(i,2)+normal(2)],[Py(i,3),Py(i,3)+normal(3)])
        end
        legend([s1,s2])
        % title('Shift direction \eta')
        fig = gcf; fig.PaperUnits = 'centimeters'; fig.PaperSize = [9 9]; fig.PaperPosition = [0 0 9 9]; xticklabels(strrep(xticklabels,'-','–')); yticklabels(strrep(yticklabels,'-','–'));
        if save_plots
            saveas(gcf, [ pwd '/' basename '/Plots/' comp_name '_shift_direction.fig' ] );
            saveas(gcf, [ pwd '/' basename '/Plots/' comp_name '_shift_direction.pdf' ] );
            saveas(gcf, [ pwd '/' basename '/Plots/' comp_name '_shift_direction.png' ] );
        end
    end

    %Projection
    if size(normal,1) == 1
        [Q,R] = qr(normal');
    else
        [Q,R] = qr(normal);
    end
    Q = sign(R(1,1))*Q;
    Ry = Py*Q';
    
    %Triangulation
    if num_obj==3
        DT_raw = delaunay(Ry(:,2),Ry(:,3));
    else
        DT_raw = delaunay(Ry(:,2),Ry(:,3),Ry(:,4));
    end

    if debug
        %Projection:
        figure
        if num_obj==3
            scatter(Ry(:,2),Ry(:,3),'.');
        else
            scatter3(Ry(:,2),Ry(:,3),Ry(:,4),'.');
            view(135,25)
        end
        title('Projection')
    
        %Projection Triangulation:
        figure
        if num_obj==3
            triplot(DT_raw,Ry(:,2),Ry(:,3));
        else
            trisurf(DT_raw,Ry(:,2),Ry(:,3),Ry(:,4));
            view(135,25)
            zlabel('f_3')
        end
        xlabel('f_1')
        ylabel('f_2')
        title('Projection Raw Triangulation')
    
        %Whole Triangulation:
        if num_obj==3
            figure
            trisurf(DT_raw,Ry(:,1),Ry(:,2),Ry(:,3));
            view(135,25)
            xlabel('f_1')
            ylabel('f_2')
            zlabel('f_3')
        else
            figure
            trisurf(DT_raw,Ry(:,1),Ry(:,2),Ry(:,3));
            view(135,25)
            xlabel('f_1')
            ylabel('f_2')
            zlabel('f_3')
            figure
            trisurf(DT_raw,Ry(:,1),Ry(:,2),Ry(:,4));
            view(135,25)
            xlabel('f_1')
            ylabel('f_2')
            zlabel('f_4')
            figure
            trisurf(DT_raw,Ry(:,1),Ry(:,3),Ry(:,4));
            view(135,25)
            xlabel('f_1')
            xlabel('f_3')
            xlabel('f_4')
            figure
            trisurf(DT_raw,Ry(:,2),Ry(:,3),Ry(:,4));
            view(135,25)
            xlabel('f_2')
            xlabel('f_3')
            xlabel('f_4')
        end
        title('Raw Triangulation')
    end

    if strcmp(clean_method,'cond')
        DT = rmtriangle( DT_raw, Ry(:,1:end), @tri_cond, threshold, comp_name, basename, debug, save_plots, debug ); % conditional number
    elseif strcmp(clean_method,'area')
        DT = rmtriangle( DT_raw, Ry(:,1:end), @tri_area, threshold, comp_name, basename, debug, save_plots,debug ); % Area
    elseif strcmp(clean_method,'off')
        DT = DT_raw;
        disp('NO CLEANING');
    else %'long'
        DT = rmtriangle( DT_raw, Ry(:,1:end), @tri_long, threshold, comp_name, basename, debug, save_plots, debug ); % longitud
    end

    if debug
        %Projection Triangulation Clean:
        figure
        if num_obj==3
            triplot(DT,Ry(:,2),Ry(:,3));
        else
            trisurf(DT,Ry(:,2),Ry(:,3),Ry(:,4));
            view(135,25)
            zlabel('f_3')
        end
        xlabel('f_1')
        ylabel('f_2')
        title('Projection Triangulation Clean')
        
        %Clean Triangulation:
        if num_obj == 3
            figure
            trisurf(DT,Py(:,1),Py(:,2),Py(:,3));
            view(135,25)
            xlabel('f_1')
            ylabel('f_2')
            zlabel('f_3')
        else
            figure
            trisurf(DT,Py(:,1),Py(:,2),Py(:,3));
            view(135,25)
            xlabel('f_1')
            ylabel('f_2')
            zlabel('f_3')
            figure
            trisurf(DT,Py(:,1),Py(:,2),Py(:,4));
            view(135,25)
            xlabel('f_1')
            ylabel('f_2')
            zlabel('f_4')
            figure
            trisurf(DT,Py(:,1),Py(:,3),Py(:,4));
            view(135,25)
            xlabel('f_1')
            xlabel('f_3')
            xlabel('f_4')
            figure
            trisurf(DT,Py(:,2),Py(:,3),Py(:,4));
            view(135,25)
            xlabel('f_2')
            xlabel('f_3')
            xlabel('f_4')
        end
        title('Clean Triangulation')
    end

end
%---------------------------SUBFUNCTIONS----------------------------------%
%--Triangulation Cleaning----%
function Volume = tri_area( tridots ) 
    num_obj = size(tridots,2);
    %-----Compute Volume-------------%
    Volume = zeros(num_obj,num_obj-1);
    for l=2:num_obj
        Volume(:,l-1) = tridots(l,:) - tridots(1,:);
    end
    Volume = det(Volume'*Volume)/factorial(num_obj);
    %-------------------------------------%
end
%----------------------------%
function kappa = tri_cond( tridots ) 
    kappa = cond(tridots);
end
%----------------------------%
function A = tri_long( tridots ) 
    if size(tridots,2) == 3
        A = max([...
            norm( tridots(1,:) - tridots(2,:) ),...
            norm( tridots(1,:) - tridots(3,:) ),...
            norm( tridots(2,:) - tridots(3,:) )...
            ]);
    else
        n = size(tridots, 1);  % Number of vertices
        A = 0;  % Initialize max distance
        for i = 1:n-1
            for j = i+1:n
                d = norm(tridots(i,:) - tridots(j,:));
                if d > A
                    A = d;
                end
            end
        end
        % A = max([...
        %     norm( tridots(1,:) - tridots(2,:) ),...
        %     norm( tridots(1,:) - tridots(3,:) ),...
        %     norm( tridots(1,:) - tridots(4,:) ),...
        %     norm( tridots(2,:) - tridots(3,:) )...
        %     norm( tridots(2,:) - tridots(4,:) )...
        %     norm( tridots(3,:) - tridots(4,:) )...
        %     ]);
    end
end
%----------------------------%
%%%---------------------------------------------------------------------%%%
