  
tic
for j=0:99 % timepoints 

    clear output_matrix simulated_data PSD
    
    for i=0:499 % for subjects 

        % open every simulated file
        example=load('/Users/admin/Documents/brainstorm_db/TutorialIntroduction/data/Subject01/S01_AEF_20131218_01_600Hz_resample/data_block001.mat');
        example.F= example.F(:,1:3000); % exclude the last timepoint
        example.Time= example.Time(1:3000); % same reason as above
        example.F=zeros([340,3000]); % clear all existing values and make them as zeros

        simulated_data=readmatrix(['/Users/admin/Downloads/SimulatedData/beta_new/0.5_250_new_new/subject_',num2str(i), '_simulation_6sec_', num2str(j), '.csv']);  % import simulated data files

        example.F(31,:)=(simulated_data'); % transpose and feed in the line 31 as the first MEG channel

        % Save the fields of structure 'example' as *individual variables* in a file
        % called data_block001_copy.mat (basically replace the original file with the new format given by the variable 'example')
        save('/Users/admin/Documents/brainstorm_db/TutorialIntroduction/data/Subject01/S01_AEF_20131218_01_600Hz_resample/data_block001.mat','-struct', 'example') % save filename & variable

        % Input files
        sFiles = {... 
            'Subject01/S01_AEF_20131218_01_600Hz_resample/data_block001.mat'};

        % Process: Power spectrum density (Welch)
        sFiles = bst_process('CallProcess', 'process_psd', sFiles, [], ...
            'timewindow',  [], ...
            'win_length',  3, ...
            'win_overlap', 50, ...
            'units',       'physical', ...  % Physical: U2/Hz
            'sensortypes', 'MLC11', ...
            'win_std',     0, ...
            'edit',        struct(...
                 'Comment',         'Power', ...
                 'TimeBands',       [], ...
                 'Freqs',           [], ...
                 'ClusterFuncTime', 'none', ...
                 'Measure',         'power', ...
                 'Output',          'all', ...
                 'SaveKernel',      0));
             
             
             PSD=load(['/Users/admin/Documents/brainstorm_db/TutorialIntroduction/data/', sFiles.FileName]);
             % to access the certain variable from this directory
             
             output_matrix(i+1,:)= squeeze(PSD.TF); %get rid of the first two dimensions as they're all ones


    end
 
    writematrix(output_matrix, ['/Users/admin/Downloads/SimulatedData/beta_new/0.5_250_new_new_PSD/PSD_', num2str(j), '_simulation_6sec.csv'])
    
    mkdir(['/Users/admin/Downloads/SimulatedData/beta_new/0.5_250_new_new_PSD_mat/',num2str(j),'/'])
    
    files2move = dir('/Users/admin/Documents/brainstorm_db/TutorialIntroduction/data/Subject01/S01_AEF_20131218_01_600Hz_resample/timefreq_psd_*.mat');
    
    for l = 1:length(files2move)
        filename=[files2move(l).folder,'/',files2move(l).name];
  
        copyfile(filename,['/Users/admin/Downloads/SimulatedData/beta_new/0.5_250_new_new_PSD_mat/',num2str(j),'/'])
        delete(filename)
    end
    
    db_reload_database('current')
end

toc

% open all files of time series (reminder matrix rows are subjects and cols
% are frequency)
for k=0:1
    
     simulated_psd(:,:)=readmatrix(['/Users/admin/Downloads/SimulatedData/frequency/0.2_250_PSD/PSD_', num2str(k),'_simulation_6sec.csv']);  


end


plot(PSD.Freqs, squeeze(mean(log(simulated_psd(1,:,:)),3)))

%plot(log(PSD.Freqs), squeeze(mean(log(simulated_psd(7,:,:)),3)))


%%



    