function persistReservesVsDemandaCSV( ...
    outDir, csvFilename, bucle, ...
    TdescInData, temps_Reserva, temps_Reserva_Acumulat)
%PERSISTRESERVESVSDEMANDACSV
% Persists reserve vs demand summary into a CSV file (append mode)
%
% Inputs:
%   outDir                 : output directory (char or string)
%   csvFilename            : CSV file name (e.g. 'reserves_vs_demanda.csv')
%   bucle                  : loop index (double)
%   TdescInData            : table with field Tipo_Prueba
%   temps_Reserva_Acumulat : accumulated reserved time (minutes)

    %% --- Ensure directory exists ---
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    csvPath = fullfile(char(outDir), csvFilename);

    %% --- Count procedures ---
    tipus_desc = TdescInData.Tipo_Prueba;

    nGastros  = sum(tipus_desc == "Gastroscopia");
    nColonos  = sum(tipus_desc == "Colonoscopia");

    %% --- Required time (HU) ---
    dDescFecha = zeros(height(TdescInData), 1);
    dDescFecha(tipus_desc == "Gastroscopia") = 30;
    dDescFecha(tipus_desc == "Colonoscopia") = 60;

    tempsNecessariHU = sum(dDescFecha);

    %% --- Deficit / surplus ---
    tempsReservaDeficit = temps_Reserva_Acumulat - tempsNecessariHU;

    %% --- Create one-row table ---
    Trow = table( ...
        bucle, nGastros, nColonos, ...
        tempsNecessariHU, temps_Reserva_Acumulat, tempsReservaDeficit, ...
        'VariableNames', { ...
            'bucle', ...
            'nGastros', ...
            'nColonos', ...
            'tempsNecessariHU', ...
            'tempsReservaAcumulat', ...
            'tempsReservaDeficit' ...
        });

    %     tipus_desc = TdescInData.Tipo_Prueba;  % "Gastroscopia" o "Colonoscopia" o "Otros"
    % nUH_bucle_Fecha = height (TdescInData);
    % dDescFecha = zeros(1, nUH_bucle_Fecha);         % durada en minuts
    % dDescFecha(tipus_desc=="Gastroscopia")  = 30;
    % dDescFecha(tipus_desc=="Colonoscopia")  = 60;

    % temps_Necessari_HU=sum(dDescFecha);
    % temps_Reserva_deficit = temps_Reserva_Acumulat - temps_Necessari_HU;
    fprintf('================ RESERVES VS DEMANDA REAL ================\n');
    fprintf('Per realitzar totes les intervencions %d H/U es necessiten %d minuts\n', ...
        nGastros+nColonos, tempsNecessariHU);
    fprintf(['Es disposa de reserva de %d minuts reservats per bloc,\n ' ...
        '        que amb un total de %d blocs,\n suposa un temps acumulat de %d ' ...
        '        minuts reservats,\n i una diferència de %d minuts.\n'], ...
        temps_Reserva, bucle, temps_Reserva_Acumulat, tempsReservaDeficit);


    %% --- Write / append to CSV ---
    if isfile(csvPath)
        writetable(Trow, csvPath, 'WriteMode', 'append');
    else
        writetable(Trow, csvPath);
    end
end
