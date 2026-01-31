function [dies_seleccio, nDiesReal, nextDay] = generaDiesLaborables( ...
    data_ini, nDies, excluded_dates)
% generaDiesLaborables
%   Laborables: dilluns–divendres (weekday 2..6)
%
%   Outputs:
%     dies_seleccio : array datetime amb dies laborables sense exclosos
%     nDiesReal     : nDies menys exclosos dins del rang
%     data_fi       : data posterior a l'últim dia de dies_seleccio

    if nargin < 3 || isempty(excluded_dates)
        excluded_dates = datetime.empty(0,1);
    end

    % --- Normalització d'entrades ---
    data_ini = dateshift(datetime(data_ini), "start", "day");
    excluded_dates = dateshift(datetime(excluded_dates), "start", "day");

    % --- Dilluns a divendres segons MATLAB ---
    DIAS_LABORABLES = 2:6;

    % --- Generació de nDies laborables consecutius ---
    dies = datetime.empty(0,1);
    d = data_ini;

    while numel(dies) < nDies
        if ismember(weekday(d), DIAS_LABORABLES)
            dies(end+1,1) = d; %#ok<AGROW>
        end
        d = d + days(1);
    end

    % --- Exclusió de dates ---
    if isempty(excluded_dates)
        dies_seleccio = dies;
    else
        dies_seleccio = dies(~ismember(dies, excluded_dates));
    end

    % --- GARANTIA: dies_seleccio sempre és array columna ---
    dies_seleccio = reshape(dies_seleccio, [], 1);

    nDiesReal = numel(dies_seleccio);

    % --- Data posterior a l'última data seleccionada ---
    if nDiesReal == 0
        nextDay = data_ini;
    else
        nextDay = dies_seleccio(end) + days(1);
    end
end
