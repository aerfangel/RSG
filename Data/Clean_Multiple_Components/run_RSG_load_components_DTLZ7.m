%% THIS SCRIPT READS AND CLEANS EACH COMPONENT SEPARATELY FOR DTLZ7.

%------------------------------------------%
clearvars
close all
%------------------------------------------%

%--------GENERATE STARTING POINTS----------%
load('separatedDTLZ7.mat','Comp1','Comp2','Comp3','Comp4')
num_comp = 4;
ref_size = 300;
%------------------------------------------%

%--------RSG COMPONENT PARAMETERS----------%
%--Size of interpolation--%
interpol_size = [700000,700000,700000,700000];
% interpol_size = [100000,100000,100000,100000];
%Other Parameters%
trimming = zeros(1,num_comp);

%--Cleaning--%
clean_method = {'long','long','long','long'};
threshold = [0.02,0.03825,0.03825,0.03825];
for comp=1:num_comp
    epsInterval{comp} = [10.2,10.2]; 
    eps_def{comp} = 0.05;
    minptsInterval{comp} = [2,2];
end
%------------------------------------------%

Py = cell(1,num_comp);
My = cell(1,num_comp);
for comp=1:num_comp
    disp(['comp ' int2str(comp) '/' int2str(num_comp)])
    close all
    %--------RSG PARAMETERS--------------------%
    %--Size of interpolation--%
    interpol_size_comp = interpol_size(comp); 
    %Other Parameters%
    trimming_comp = trimming(comp);
    
    %--Cleaning--%
    clean_method_comp = clean_method{comp};
    threshold_comp = threshold(comp);
    
    %COMPONENT DETECTION PARAMETERS
    epsInterval_comp = epsInterval{comp}; %TESTS!!! MODIFIED VALUE
    eps_def_comp = eps_def{comp};
    minptsInterval_comp = minptsInterval{comp};
    
    %--------GENERATE STARTING POINTS----------%
    eval(['Py{' int2str(comp) '} = Comp' int2str(comp) ';'])
    %------------------------------------------%
    
    %-------------------RSG--------------------%
    [~,My{comp}] = RSG( Py{comp}, interpol_size_comp, ref_size, clean_method_comp, threshold_comp, epsInterval_comp, eps_def_comp, minptsInterval_comp,trimming_comp,false,'means');
    %------------------------------------------%

end

%--------------FINAL TRIMMING--------------%
Iy = [];
for comp=1:num_comp
    Iy = [Iy; My{comp}];
end
max_iters = 5000;
[~,Z] = kmeans(Iy,ref_size,'Options',statset('MaxIter',max_iters));
%------------------------------------------%
