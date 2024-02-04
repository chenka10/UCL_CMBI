function WriteLineToCSV(data, filename)    
    fileID = fopen(filename, 'a');
    
    fprintf(fileID, '%s', num2str(data{1}));
    for i = 2:numel(data)
        fprintf(fileID, ',%s', num2str(data{i},5));
    end
    
    fprintf(fileID, '\n');

    fclose(fileID);
end