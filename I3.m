filename1 = ("./output/ACF.dat");
data1 = readmatrix(filename1);
f1 = figure;
plot(data1(:,2));
samplenum = 8192;
xlim([0, samplenum])

filename2 = ('./output/PSD_re.dat');
data2 = readmatrix(filename2);
f2 = figure;
plot(data2(:, 2));
xlim([samplenum samplenum*2])
%ylim([0 1e10])