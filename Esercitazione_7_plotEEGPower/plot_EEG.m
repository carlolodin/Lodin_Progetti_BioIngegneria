function plot_EEG(names, data, sampling_rate)
    % Plot EEG power
    % names: cell array of channel names
    % data: matrix of EEG data (channels x time)
    % sampling_rate: sampling rate in Hz

    % Check input dimensions
    if size(data, 1) ~= length(names)
        error('Number of channels in data does not match number of names provided.');
    end

    power = zeros(size(data, 1), floor(size(data, 2)/sampling_rate));
    
    figure;
    hold on;
    for channel = 1:size(data, 1)
        % Compute power for each channel
        num_segments = floor(size(data, 2) / sampling_rate);
        for seg = 1:num_segments
            start_idx = (seg - 1) * sampling_rate + 1;
            end_idx = seg * sampling_rate;
            segment_data = data(channel, start_idx:end_idx);
            power(channel, seg) = mean(segment_data .^ 2); % Power calculation
        end
        subplot(3, 3, channel);
        stem(1:num_segments, power(channel, :), 'filled');
        title(['Power of ' names{channel}]);
        xlabel('Time (seconds)');
        ylabel('Power ($\mu$V$^2$)', 'Interpreter', 'latex');
        grid on;
    end
end