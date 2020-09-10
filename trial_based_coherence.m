function cohr = trial_based_coherence(x,y,dT,duration)

K = size(x,1); %Define the number of trials.
N = size(x,2); %Define the number of indices per trial.
Sxx = zeros(K,N); %Create variables to save the spectra.
Syy = zeros(K,N);
Sxy = zeros(K,N);
for k=1:K %Compute the spectra for each trial.
    Sxx(k,:) = 2*dT^2/duration * fft(x(k,:)) .* conj(fft(x(k,:)));
    Syy(k,:) = 2*dT^2/duration * fft(y(k,:)) .* conj(fft(y(k,:)));
    Sxy(k,:) = 2*dT^2/duration * fft(x(k,:)) .* conj(fft(y(k,:)));
end
halfN = floor(N/2);
Sxx = Sxx(:,1:halfN+1); %Ignore negative frequencies.
Syy = Syy(:,1:halfN+1);
Sxy = Sxy(:,1:halfN+1);
Sxx = nanmean(Sxx,1); %Average the spectra across trials.
Syy = nanmean(Syy,1);
Sxy = nanmean(Sxy,1);
cohr = abs(Sxy ./ (sqrt(Sxx) .* sqrt(Syy)));

end