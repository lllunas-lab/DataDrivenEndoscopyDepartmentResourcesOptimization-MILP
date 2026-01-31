function TdebugMILP = afegeixMILPDebug( ...
    TdebugMILP, bucle, nP, nR, nT, fval, exitflag, output, elapsed)
%AFEGEIXMILPDEBUG Logs selected MILP debug information
%
% Columns:
%   bucle, nP, nR, nT, fval, exitflag, elapsed
%   + one column per field of output (stored as cell)

    % --- Base scalar fields ---
    Tbase = table( ...
        bucle, nP, nR, nT, fval, exitflag, elapsed, ...
        'VariableNames', ...
        {'bucle','nP','nR','nT','fval','exitflag','elapsed'} ...
    );

    % --- Output fields (force cell storage) ---
    outFields = fieldnames(output);
    nFields   = numel(outFields);

    outData = cell(1, nFields);
    for k = 1:nFields
        outData{k} = {output.(outFields{k})};
    end

    Tout = cell2table(outData, 'VariableNames', outFields');

    % --- Full row ---
    Trow = [Tbase Tout];

    % --- Append or initialise ---
    if isempty(TdebugMILP)
        TdebugMILP = Trow;
    else
        TdebugMILP = [TdebugMILP; Trow];
    end
end
