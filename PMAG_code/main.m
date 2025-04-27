clear;
close all;
clc;

addpath('funcs_measures','hybridSelection','funs');
dataname='XXX';
load dataname;

fprintf(fid, 'Result on %s \r\n', dataname);
ds = strcat(datapath,ds.data_file);
for projectionDim= [10:10:100]
    for r = [1.2,1.4,1.6,1.8,2,4,6,8,10,12]
        for cnt=1:20
            count=count+1;
            [result,alpha,runtime,obj] = run_PMAG(ds,projectionDim,r);
            fprintf(fid, 'ACC = %0.4f NMI = %0.4f Purity = %0.4f RunTime= %f r = %0.1f projectionDim=%d \r\n',result(1:3),runtime,r,projectionDim);
        end

    end
end