% Analyze the areaBetweenCurves for all files
fileList = dir('Data'); % list of the data folder files
nFiles = length(fileList); % total files

results = [];

for i = 1:nFiles
    filename = fileList(i).name;
    if length(filename) > 3
        if strcmp(filename(end-4:end), '.json') == true
            disp(filename)
        else continue
        end
    else continue
    end

   filenameShort = filename(1:end-5);
   [G2N, G2D, D2N] = areaBetweenCurves(filenameShort);
   sampleID = filenameShort;
   results2 = table(string(sampleID), G2N, G2D, D2N);
   results = [results; results2];
    %results2 = analyzeFile(filename(1:end-5));
    %results = [results; results2];
end
close all
%writetable(results, 'tabulatedResults.xlsx');