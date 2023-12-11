
function [W,FT] = rate_M(cluster_file,time_file,varargin)

    bin_size = 10;

    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'bin_size'
                bin_size = varargin{ii+1};
        end 
    end 
        
    %% read in data

    spike_clusters = readNPY(cluster_file);
    spike_times = readNPY(time_file);
    
    total_time = max(spike_times); 
    time_bins = 0:bin_size:total_time;
    
    unique_clusters = unique(spike_clusters); 
    num_clusters = length(unique_clusters); 
    spike_matrix = zeros(num_clusters,length(time_bins)-1); 
    
    %% populate spike matrix

    for i = 1:num_clusters  
        curr_cluster = unique_clusters(i); 
        curr_times = spike_times(spike_clusters == curr_cluster); 
    
        spike_counts = histcounts(curr_times, time_bins); 
        spike_matrix(i,:) = spike_counts; 
    end
    
    %% convolution
  
    kernel_width = 50; 
    sigma = 10;  
    
    % Define Gaussian function
    gauss_kernel = fspecial('gaussian', [1, kernel_width], sigma/bin_size);
    
    % Normalize to ensure the area under the kernel is 1
    gauss_kernel = gauss_kernel / sum(gauss_kernel);
    
    W = conv2(spike_matrix,gauss_kernel,'same'); 

    FT = fft(W')

end
    

