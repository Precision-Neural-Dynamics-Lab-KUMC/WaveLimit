function new_waveforms = perform_spike_alignment_for_gpu(waveforms)
% perform_spike_alignment for WaveLimit, v1.3
% Adam Rouse, 4/24/20
% function for aligning spikes based on the peak or trough of the spike
% waveforms
% INPUT: file_path, directory where .mat files created with
% create_spike_mat is located
% OUTPUT: a _aligned.mat file is saved with shifted spike waveforms,
% original time stamps, corrected time stamps based on the alignment shift,
% and the mean waveforms for each original cluster specified by the user


w_pre=8;                       %number of pre-event data points stored before threshold
interp_factor = 4;              %interpolation factor, for 4, new peak/trough can occur at 0.25, .5, .75 samples in addition to whole samples
interp_offset = 5;             %number of samples before and after the alignment point to use for interpolation
max_slide_samples = 2;          %Max number of samples that each waveform can be shifted relative to the alignment point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g_waveforms = gpuArray(waveforms);
  
if size(waveforms,2)>100  %Need a at least 100 spikes on the channel to do the alignment of spike snippets, smaller numbers cause errors and warnings
    
    n_samples = size(waveforms,1);
    %Test the slopes around the threshold time to determine if there was a
    %positve or negative threshold
    n_test_waveforms = min(100, size(waveforms,2));  %Take first 100 waveforms (or all if <100)
    test_slope = diff(waveforms([w_pre w_pre+1],1:n_test_waveforms), [], 1);
    if min(test_slope) > 0
        detection = 'pos';
    else
        detection = 'neg';
    end
    
    %Find the time point in each waveform that has the max/min value
    if strcmpi(detection, 'neg')
        [~, best_i] = min(g_waveforms,[],1);
    else
        [~, best_i] = max(g_waveforms,[],1);
    end
    %Choose the time point with the most max/min's to use as the alignment point
    alignment_index = mode(best_i(best_i > interp_offset & best_i < (n_samples-interp_offset+1)));
    
    
    [xq,yq] = meshgrid(1:size(g_waveforms,2), (alignment_index-interp_offset):(1/interp_factor):(alignment_index+interp_offset) );
    interp_spike = interp2( g_waveforms,xq,yq, 'cubic');
    
    
    %The original aligment point occurs at interp_offset*interp_factor+1,
    %we're going to find the peak/trough that is within 2
    %original samples (max_slide_samples) before or after the alignment point
    if strcmpi(detection, 'neg')
        [~, best_index] = min(interp_spike(((interp_offset*interp_factor+1)-(max_slide_samples*interp_factor)):((interp_offset*interp_factor+1)+max_slide_samples*interp_factor),:),[],1);
    else
        [~, best_index] = max(interp_spike(((interp_offset*interp_factor+1)-(max_slide_samples*interp_factor)):((interp_offset*interp_factor+1)+max_slide_samples*interp_factor),:),[],1);
    end
    %Determine number of samples to shift original spike, 0 is no shift, +2
    %would mean shift spike 2 samples because origninal has a spike peak that is
    %2 samples after the global alignment time
    new_offset = (best_index - (max_slide_samples*interp_factor+1))./interp_factor;
    
    
    %Find waveforms with shifted time, interpolation required when peak/trough
    %does not happen on whole number sample, for example, if original
    %waveform is samples 1:32 and the peak is 1.5 later than the alignment
    %point, then we're going to interpolate to find the points at
    %(3:30)+1.5 or 4.5:1:31.5
    new_waveforms = gpuArray(zeros(n_samples - 2*max_slide_samples, size(waveforms,2)));
    unique_offsets = unique(new_offset);
    for n = 1:length(unique_offsets)
        %TODO see if this can be improved for speed
        if sum(unique_offsets(n)==new_offset)>2 %If it's more than one waveform at the current offset,  we use interp2 on a matrix
            if mod(unique_offsets(n),1)==0  %If offset is an interger, don't have to interpolate
                new_waveforms(:,unique_offsets(n)==new_offset) = g_waveforms(((1+max_slide_samples):(n_samples-max_slide_samples))+unique_offsets(n),unique_offsets(n)==new_offset);
            else
                [xq,yq] = meshgrid(1:sum(unique_offsets(n)==new_offset), ((1+max_slide_samples):(n_samples-max_slide_samples))+unique_offsets(n));
                 new_waveforms(:,unique_offsets(n)==new_offset) = interp2(g_waveforms(:,unique_offsets(n)==new_offset), xq ,yq, 'cubic');
            end
        else
            yq = ((1+max_slide_samples):(n_samples-max_slide_samples))+unique_offsets(n);
            new_waveforms(:,unique_offsets(n)==new_offset) = interp1(waveforms(:,unique_offsets(n)==new_offset), yq, 'V5cubic');
        end
    end 
else
    new_waveforms = waveforms;
end