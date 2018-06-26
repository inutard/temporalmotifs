filenames = { ...
    'email-Eu-core-temporal', ...
    'sx-mathoverflow', ...
    'sx-superuser', ...
    'CollegeMsg'};
% 'email-Eu-core-temporal-data.txt'
% 'CollegeMsg-data.txt'
% 'sx-mathoverflow-data.txt'
% 'sx-superuser-data.txt'

for i = 1:length(filenames)
    filename = filenames{i};
    
    subplot(length(filenames), 4, 4 * i - 3);
    fn = strcat(filename, '-edge-data.txt');
    fileID = fopen(fn, 'r');
    A = fscanf(fileID, '%f');
    plot(cumsum(A));
    title(strcat('#edges vs time. ', filename));
    
    subplot(length(filenames), 4, 4 * i - 2);
    fn = strcat(filename, '-data2.txt');
    fileID = fopen(fn, 'r');
    A = fscanf(fileID, '%f');
    plot(cumsum(A));
    title('#motifs vs time');
    
    subplot(length(filenames), 4, 4 * i - 1);
    fn = strcat(filename, '.txt');
    data = dlmread(fn);
    times = data(:,3);
    times = times - times(1);
    times = mod(times, 24 * 60 * 60) / 60 / 60; 
    histogram(times);
    title('#edges vs time of day');
    xlabel('time (h)');
    xlim([0 24]);
    
    subplot(length(filenames), 4, 4 * i);
    fn = strcat(filename, '.txt');
    data = dlmread(fn);
    times = data(:,3);
    times = (times - times(1)) / 60 / 60; 
    histogram(times);
    title('#edges vs time (hist)');
    xlabel('time (h)');
end

