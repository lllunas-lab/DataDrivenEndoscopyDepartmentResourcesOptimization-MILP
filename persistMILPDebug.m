function persistMILPDebug(TdebugMILP, outDir, baseFilename)

    arguments
        TdebugMILP table
        outDir {mustBeText}
        baseFilename (1,:) char
    end

    if ~isfolder(outDir)
        mkdir(outDir);
    end

    % Force char path (IMPORTANT)
    basePath = fullfile(char(outDir), baseFilename);

    % MAT
    save([basePath '.mat'], 'TdebugMILP', '-v7.3');

    % JSON
    try
        jsonText = jsonencode(TdebugMILP);
        fid = fopen([basePath '.json'],'w');
        fwrite(fid, jsonText, 'char');
        fclose(fid);
    catch ME
        warning('JSON export failed: %s', ME.message);
    end

    % CSV (scalar summary)
    varNames = TdebugMILP.Properties.VariableNames;
    isScalarNumeric = false(size(varNames));
    for i = 1:numel(varNames)
        col = TdebugMILP.(varNames{i});
        isScalarNumeric(i) = isnumeric(col) || islogical(col);
    end

    Tcsv = TdebugMILP(:, isScalarNumeric);
    writetable(Tcsv, [basePath '_summary.xlsx']);
end
