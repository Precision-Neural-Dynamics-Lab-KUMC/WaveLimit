n_wave = 10000;
num_time_pts = 20;
d = -5:0.1:5;
sigma = 1;
sigma2 = 2.5;

waveforms = randn(n_wave,num_time_pts);

% dist = sort( sqrt( sum(waveforms.^2,2) ));
dist = sort(  sum(waveforms.^2,2) );
diff_dist = diff(dist);
[b,a] = butter(1,.001,'low');
%Low-pass to improve estimating the minimum distance rate change
smooth_diff_dist = filter(b,a,diff_dist); 

figure; plot(dist)
figure; plot(diff_dist)
figure; plot(smooth_diff_dist)
smooth_diff_dist(1:100) = NaN;
[~,mode_wave] = min(smooth_diff_dist);


chi2_dist = chi2inv(linspace(0,1,n_wave),num_time_pts);

n_wave_at_mode = n_wave*chi2cdf(num_time_pts-2,num_time_pts);

figure; plot(chi2_dist)

figure; plot(diff(chi2_dist))

[~,model_mode] = min(diff(chi2_dist));


est_25th = round((0.25/chi2cdf(num_time_pts-2,num_time_pts))*mode_wave);

figure;
plot(dist(1:est_25th) )
hold on
plot(chi2_dist(1:est_25th) )

chi2fun = @(x,xdata) chi2inv(x(2).*xdata,num_time_pts)./x(1);  %chi2inv(P,v) p = probability, v = degrees of freedom, in our case the number of time points

 x_values = linspace(0.25./est_25th, 0.25, est_25th);
 y_values = dist(1:est_25th)';

[x,~] = lsqcurvefit(chi2fun, [chi2inv(.25, num_time_pts)./(dist(est_25th)), 1], x_values(1:round(length(x_values)/100):end), y_values(1:round(length(x_values)/100):end), [0 0], [1e10 1e10]);

