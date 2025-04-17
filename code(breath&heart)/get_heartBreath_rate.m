function [breathRate, heartRate] = get_heartBreath_rate(target_rangeProfile, slowtime_fs)

target_num = length(target_rangeProfile(1,:));
breathRate = zeros(1, target_num);
heartRate  = zeros(1, target_num);
useFramesNum = 1200;
for kk = 1:target_num
    target_profile = target_rangeProfile(:,kk);
    dcRemove_ag = angle(target_profile);
    unwrap_dcRemove_ag = unwrap(dcRemove_ag);
    diff_ag = unwrap_dcRemove_ag(1:end-1) - unwrap(dcRemove_ag(2:end));
    agl = diff_ag;

    agl = agl - mean(agl);
    breath_wave = filter(RR_BPF20, agl);
    breath_fft = abs(fftshift(fft(breath_wave)));
    breath_fft = breath_fft(ceil(end/2):end);
    [~, breath_index] = max(breath_fft);
    breath_hz = breath_index * (slowtime_fs / 2 / length(breath_fft));
    breathRate(kk) = ceil(breath_hz * 60);

    heart_wave = filter(HR_BPF20, agl);
    heart_fft = abs(fftshift(fft(heart_wave)));
    heart_fft = heart_fft(ceil(end/2):end);
    [~, heart_index] = max(heart_fft);
    heart_hz = heart_index * (slowtime_fs / 2 / length(heart_fft));
    heartRate(kk) = ceil(heart_hz * 60);

    xAxis = linspace(0, slowtime_fs/2, length(breath_fft)) * 60;
    t_axis = linspace(0, 1/slowtime_fs * length(target_profile), length(target_profile)-1 );
    
    figure;
    subplot(221); plot(t_axis, breath_wave); title('breath'); xlabel('time(s)'); ylabel('A'); grid on;
    subplot(223); plot(xAxis, breath_fft); title('breath'); xlabel('f'); ylabel('A'); grid on;
    subplot(222); plot(t_axis, heart_wave); title('heart'); xlabel('time(s)'); ylabel('A'); grid on;
    subplot(224); plot(xAxis, heart_fft); title('heart'); xlabel('f'); ylabel('A'); grid on;
end
