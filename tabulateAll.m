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

   
    results2 = analyzeFile(filename(1:end-5));
    results = [results; results2];
end

writetable(results, 'tabulatedResults.xlsx');