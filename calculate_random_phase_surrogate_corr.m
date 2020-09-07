% For each of a pair of time series of the same length, generate a set of
% random surrogate time series with the same autocovariance function, by
% using random phase spectra while keeping the amplitude spectrum the same
% as the original time series.
% When calculating correlations between surrogate time series, no surrogate
% is used twice, so that the correlation samples are independent.
% -actual_time_series: 2 x number of time points; the 2 time series
% -correlations_between_surrogates: number_surrogate_pairs x 1
% -surrogate_time_series: number of time points x number_surrogate_pairs x
% 2; all the surrogate time series for each of the two input time series
function surrogate_corr = calculate_random_phase_surrogate_corr(actual_time_series,number_surrogate_pairs,varargin)

incorporate_trial_avg = false;

if nargin == 3
    trial_avg_time_series = varargin{1};
    actual_time_series = actual_time_series - trial_avg_time_series;
    incorporate_trial_avg = true;
end

imaginary_component_tolerance=1e-8;

num_time_points=size(actual_time_series,2);
if mod(num_time_points,2)==1 % odd
    num_repeated_frequencies=(num_time_points-1)/2;
    indices_negative_frequencies=1+num_repeated_frequencies+(1:num_repeated_frequencies);
elseif mod(num_time_points,2)==0 % even
    num_repeated_frequencies=(num_time_points-2)/2;
    indices_negative_frequencies=1+num_repeated_frequencies+1+(1:num_repeated_frequencies);
end
indices_positive_frequencies=1+(1:num_repeated_frequencies);

amplitude_spectrum=cell(size(actual_time_series,1),1);
phase_spectrum=cell(size(actual_time_series,1),1);
for series_i=1:size(actual_time_series,1)
    time_series=actual_time_series(series_i,:)-mean(actual_time_series(series_i,:));
    
    fourier_transform=fft(time_series);
    amplitude_spectrum{series_i}=abs(fourier_transform);
    phase_spectrum{series_i}=angle(fourier_transform);
end


random_phases=-pi+2*pi*rand(2,num_repeated_frequencies,number_surrogate_pairs);

phase_spectrum_mat = vertcat(phase_spectrum{:});
random_phase_spectrum_mat = repmat(phase_spectrum_mat,[1 1 number_surrogate_pairs]);
random_phase_spectrum_mat(:,indices_positive_frequencies,:)=random_phases;
random_phase_spectrum_mat(:,indices_negative_frequencies,:)=-flip(random_phases,2);

amplitude_spectrum_mat = vertcat(amplitude_spectrum{:});
amplitude_spectrum_mat = repmat(amplitude_spectrum_mat,[1 1 number_surrogate_pairs]);

random_phase_fourier_transform_mat=amplitude_spectrum_mat.*exp(1i*random_phase_spectrum_mat);
all_surrogate_mat=ifft(random_phase_fourier_transform_mat,[],2);

if all(imag(all_surrogate_mat)<imaginary_component_tolerance,'all') 
    all_surrogate_mat = real(all_surrogate_mat);
else
    disp('Warning: surrogate data are complex...')
end

if incorporate_trial_avg
    trial_avg_time_series = repmat(trial_avg_time_series,[1 1 number_surrogate_pairs]);
    all_surrogate_mat = all_surrogate_mat + trial_avg_time_series;
end

A = reshape(squeeze(all_surrogate_mat(1,:,:)),num_time_points,number_surrogate_pairs);
B = reshape(squeeze(all_surrogate_mat(2,:,:)),num_time_points,number_surrogate_pairs);

R = corr(A,B);
surrogate_corr = diag(R);

end