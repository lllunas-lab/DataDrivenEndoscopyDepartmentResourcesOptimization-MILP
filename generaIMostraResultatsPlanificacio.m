function [Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat,...
    temps_Reserva, miniAvail, miniReservaType, nMiniAvailDisponible] = ...
    generaIMostraResultatsPlanificacio( ...
    offset_hour, nDies, slotRangeDia, Delta, slotsDia, dies_seleccio, ...
    nP, nR, d, tipus, T2p, ...
    avail, reserveType, ...
    x_opt, u_opt, idle_opt, delta_opt, mu, ...
    Total_Files, nP_total, nUH_total, ...
    nomDirResults, inP_ini, ...
    last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva)


 % function [Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva] = ...
 %    generaIMostraResultatsPlanificacio( ...
 %        % Bàsics / temps
 %        offset_hour, nDies, slotRangeDia, Delta, slotsDia, dies_seleccio, ...
 %        % Dimensions i dades
 %        nP, nR, d, tipus, T2p, ...
 %        % Disponibilitat / reserves
 %        avail, reserveType, ...
 %        % Solució MILP
 %        x_opt, u_opt, idle_opt, delta_opt, mu, ...
 %        % Debug/prints
 %        Total_Files, nP_total, nUH_total, ...
 %        % Persistència
 %        nomDirResults, inP_ini, ...
 %        % Estat acumulat (modificable)
 %        last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva, ...
 %        % Callbacks
 %        nextRowResFcn)

%GENERAI MOSTRA RESULTATS PLANIFICACIO
% Encapsula:
%  - Presentació de l'assignació per dia (lectura de x_opt)
%  - Construcció de la taula Res (procediments)
%  - Afegir reserves (avail/reserveType) i intervals lliures (NOT_USED)
%  - Estadístiques i visualització (matriu per dia + llistat detallat)
%  - Persistència a CSV (parcial i acumulat)
%
% INPUT (struct "in") - camps necessaris:
%   in.nDies_real
%   in.slotRangeDiaSafe           % function handle: @(dia) ...
%   in.nP, in.nR
%   in.x_opt                      % (nP x nR x nT)
%   in.d                          % durades procediments (min)
%   in.tipus                      % cellstr o string array indexable per i
%   in.Delta                      % minuts per slot
%   in.slotsDia
%   in.dies_seleccio              % datetime array (1..nDies_real)
%   in.avail                       % (nR x nT_avail) 1=disponible, 0=reservat
%   in.reserveType                 % (nR x nT_avail) string: tipus reserva o ""
%   in.T2p                         % taula amb camps: Prioridad, ID, FECHA, HORA
%   in.nextRowRes                  % function handle: @(Res,fila,growBlock) ...
%
%   (per estadístiques/prints finals)
%   in.Total_Files, in.nP_total, in.nUH_total
%   in.u_opt, in.idle_opt, in.delta_opt, in.mu
%
%   (persistència)
%   in.nomDirResults               % string o char
%   in.inP_ini                     % per nom fitxer parcial
%
%   (acumuladors / variables que es modifiquen)
%   in.last_date_non_prog_proc_filter  % datetime (es pot inicialitzar a NaT)
%   in.temps_Reserva_Acumulat          % numeric
%   in.temps_Reserva                   % numeric
%
% OUTPUT:
%   Res  : taula de resultats (ordenada)
%   out  : struct amb:
%          out.assignat, out.assignatsPerDia, out.slotsDisponiblesPerDia,
%          out.slotsLliuresPerDia, out.totalAssignats, out.totalNoAssignats,
%          out.last_date_non_prog_proc_filter (ACTUALITZAT),
%          out.temps_Reserva_Acumulat (ACTUALITZAT),
%          out.nomFitxerParcial, out.nomFitxerTotal



    %% --------- 1) Presentació assignació per dia ----------
    for dia = 1:nDies
        fprintf('=== DIA %d ===\n', dia);
        tRange = slotRangeDia(dia);

        for i = 1:nP
            found = false;
            for k = 1:nR
                for t = tRange
                    if x_opt(i,k,t) > 0.5
                        found = true;
                        slotLocal = t - tRange(1);
                        startMin  = slotLocal * Delta;
                        endMin    = startMin + d(i);
                        fprintf('  Proc %2d (%s, %2d min): Sala %d, de %3d a %3d min\n', ...
                            i, string(tipus(i)), d(i), k, startMin, endMin);
                    end
                end
            end
            if ~found && dia==nDies
                fprintf('  Proc %2d NO assignat en els %d dies\n', i, nDies);
            end
        end
        fprintf('\n');
    end

    %% generar resultats per simular de nou estadística %%
    %% 11. Construcció de la taula de resultats Res (procediments)
    assignat           = false(1, nP);            % si el procediment i s'ha assignat a algun dia
    assignatsPerDia    = zeros(1, nDies);    % # procediments assignats per dia
    slotsPerDia        = nR * slotsDia;           % capacitat total de slots per dia: nT
    slotsDisponiblesPerDia = zeros(1, nDies);
    for dia = 1:nDies
        tRange = slotRangeDia(dia);
        slotsDisponiblesPerDia(dia) = sum(avail(:, tRange), 'all');  % només els slots avail=1
    end
    slotsLliuresPerDia = slotsDisponiblesPerDia;  % inicialment, tots disponibles són lliures

    % --- 1) Mida inicial i mida del bloc d'ampliació ---
    initRows  = max(nP, 400);   % com a mínim nP, però posa un mínim raonable
    growBlock = 50;            % o 5000 segons volum

    % --- 2) Crear Res preallocada amb initRows files ---
    Res = table( ...
        NaN(initRows,1), ...                        % ID_PROC (NaN per missing)
        NaN(initRows,1), ...                        % SALA
        NaT(initRows,1), ...                        % FECHA
        duration(zeros(initRows,1),0,0), ...        % HORA
        duration(zeros(initRows,1),0,0), ...        % HORA_FI
        strings(initRows,1), ...                    % TIPUS
        strings(initRows,1), ...                    % PRIORITAT
        NaN(initRows,1), ...                        % ID
        NaT(initRows,1), ...                        % FECHA_ORG
        duration(zeros(initRows,1),0,0), ...        % HORA_ORG
        duration(zeros(initRows,1),0,0), ...        % DUR_RES_MIN
        duration(zeros(initRows,1),0,0), ...        % DUR_IDLE_MIN
        duration(zeros(initRows,1),0,0), ...        % DUR_PROC
        'VariableNames', {'ID_PROC','SALA','FECHA','HORA','HORA_FI','TIPUS','PRIORITAT', ...
        'ID','FECHA_ORG','HORA_ORG','DUR_RES_MIN','DUR_IDLE_MIN','DUR_PROC'} );

    fila = 0;

    for dia = 1:nDies
        tRange = slotRangeDia(dia);

        for i = 1:nP
            for k = 1:nR
                for t = tRange
                    if x_opt(i,k,t) > 0.5
                        % Marquem com assignat i actualitzem indicadors
                        if ~assignat(i)
                            assignat(i) = true;
                            assignatsPerDia(dia) = assignatsPerDia(dia) + 1;
                            slotsUsats = ceil(d(i) / Delta);
                            slotsLliuresPerDia(dia) = slotsLliuresPerDia(dia) - slotsUsats;
                        end

                        % Càlcul d'hores
                        slotLocal = t - tRange(1);      % 0..slotsDia-1
                        startMin  = slotLocal * Delta;  % minuts des de l'inici del dia
                        endMin    = startMin + d(i);

                        horaIni = minutes(startMin) + offset_hour;
                        horaFi  = minutes(endMin)  + offset_hour;
                        duradaEvent = horaFi - horaIni;

                        % Afegim fila a la taula Res
                        [Res, fila] = nextRowRes(Res, fila, growBlock);

                        Res.ID_PROC(fila)    = i;
                        Res.SALA(fila)       = k;
                        Res.FECHA(fila)      = dies_seleccio(dia);
                        Res.HORA(fila)       = horaIni;
                        Res.HORA_FI(fila)    = horaFi;
                        Res.TIPUS(fila)      = string(tipus(i));
                        Res.PRIORITAT(fila)  = string(T2p.Prioridad(i));
                        Res.ID(fila)         = T2p.ID(i);
                        Res.FECHA_ORG(fila)  = T2p.FECHA(i);
                        Res.HORA_ORG(fila)   = T2p.HORA(i);
                        Res.DUR_RES_MIN(fila)= duration(0,0,0);
                        Res.DUR_IDLE_MIN(fila)= duration(0,0,0);
                        Res.DUR_PROC(fila)   = duradaEvent;

                        if (last_date_non_prog_proc_filter < dies_seleccio(dia))
                            last_date_non_prog_proc_filter = dies_seleccio(dia);
                        end
                    end
                end
            end
        end
    end

    %% 12. Afegir a Res també els slots reservats per Urgències i Hospitalitzacions

    % IMPORTANT: usem la mida REAL d'avail/reserveType (evitem index out of bounds)
    nT_avail = size(avail, 2);
    nDies_avail = ceil(nT_avail / slotsDia);
    slotRangeDiaAvail = @(dia) ((dia-1)*slotsDia + 1) : min(dia*slotsDia, nT_avail);

    for dia = 1:nDies_avail
        tRange = slotRangeDiaAvail(dia);

        % Protecció per si dies_seleccio és més curt (clamp)
        diaFecha = min(dia, numel(dies_seleccio));

        for k = 1:nR
            t = tRange(1);
            while t <= tRange(end)

                % Hi ha reserva en aquest slot?
                if avail(k,t) == 0 && reserveType(k,t) ~= ""
                    tipusReserva = reserveType(k,t);

                    % Agrupem slots contigus amb la mateixa reserva
                    tStart = t;
                    while t <= tRange(end) && avail(k,t) == 0 && reserveType(k,t) == tipusReserva
                        t = t + 1;
                    end
                    tEnd = t - 1;

                    % Càlcul d'hores (minuts des de l'inici del dia)
                    slotLocalStart = tStart - tRange(1);      % 0..slotsDia-1
                    startMin       = slotLocalStart * Delta;
                    duracioMin     = (tEnd - tStart + 1) * Delta;
                    endMin         = startMin + duracioMin;

                    horaIni = minutes(startMin) + offset_hour;
                    horaFi  = minutes(endMin)  + offset_hour;
                    duradaEvent = horaFi - horaIni;

                    % Afegim fila COMPLETA a Res (evita warnings)
                    % Nova fila RESERVA
                    [Res, fila] = nextRowRes(Res, fila, growBlock);

                    Res.ID_PROC(fila)      = NaN;
                    Res.SALA(fila)         = k;
                    Res.FECHA(fila)        = dies_seleccio(diaFecha);
                    Res.HORA(fila)         = horaIni;
                    Res.HORA_FI(fila)      = horaFi;
                    Res.TIPUS(fila)        = "RESERVA";
                    Res.PRIORITAT(fila)    = "RESERVA";
                    Res.ID(fila)           = NaN;                     % millor que 0
                    Res.FECHA_ORG(fila)    = dies_seleccio(diaFecha);
                    Res.HORA_ORG(fila)     = horaIni;
                    Res.DUR_RES_MIN(fila)  = duradaEvent;
                    Res.DUR_IDLE_MIN(fila) = duration(0,0,0);
                    Res.DUR_PROC(fila)     = duration(0,0,0);

                else
                    t = t + 1;
                end
            end
        end
    end

    %% 12b. Afegir a Res també els slots Lliures

    % 1) Matriu d'ocupació per sala i slot global
    nT = size(x_opt, 3);
    ocupat = zeros(nR, nT);

    for i = 1:nP
        slotsUsats = ceil(d(i) / Delta);
        for k = 1:nR
            for t = 1:nT
                if x_opt(i,k,t) > 0.5
                    tEnd = min(t + slotsUsats - 1, nT);
                    ocupat(k, t:tEnd) = 1;
                end
            end
        end
    end

    % 2) Afegim intervals LLIURES a Res (avail=1 i ocupat=0)
    for dia = 1:nDies
        tRange = slotRangeDia(dia);
        for k = 1:nR
            freeMask = (avail(k, tRange) == 1) & (ocupat(k, tRange) == 0);

            % Detectar runs (intervals) de 1s a freeMask
            dmask = diff([0 freeMask 0]);
            starts = find(dmask == 1);
            ends   = find(dmask == -1) - 1;

            for r = 1:length(starts)
                tStartLocal = starts(r) - 1;        % 0..slotsDia-1
                tEndLocal   = ends(r);              % 1..slotsDia

                startMin = tStartLocal * Delta;
                endMin   = tEndLocal * Delta;       % fi en minuts

                horaIni = minutes(startMin) + offset_hour;
                horaFi  = minutes(endMin)  + offset_hour;
                duradaEvent = horaFi - horaIni;

                [Res, fila] = nextRowRes(Res, fila, growBlock);
                Res.ID_PROC(fila)      = NaN;
                Res.SALA(fila)         = k;
                Res.FECHA(fila)        = dies_seleccio(dia);
                Res.HORA(fila)         = horaIni;
                Res.HORA_FI(fila)      = horaFi;
                Res.TIPUS(fila)        = "NOT_USED";
                Res.PRIORITAT(fila)    = "NOT_USED";
                Res.ID(fila)           = NaN;
                Res.FECHA_ORG(fila)    = dies_seleccio(dia);
                Res.HORA_ORG(fila)     = horaIni;
                Res.DUR_RES_MIN(fila)  = duration(0,0,0);
                Res.DUR_IDLE_MIN(fila) = duradaEvent;
                Res.DUR_PROC(fila)     = duration(0,0,0);
            end
        end
    end

    %% retornem miniavail,minireserva i si hi ha slots lliures per tal de fer el minimilp
    % | Condició             | miniAvail | miniReservaType |
    % | -------------------- | --------: | --------------- |
    % | `avail=0`            |         0 | `"RESERVA"`     |
    % | `avail=1 & ocupat=1` |         0 | `"PROGRAM"`     |
    % | `avail=1 & ocupat=0` |         1 | `"NOT_USED"`    |

    % Mateixa mida que avail
    miniAvail       = false(nR, nT);
    miniReservaType = strings(nR, nT);

    % 1) Base: tot el que NO està disponible (avail=0) és RESERVA / bloqueig
    miniReservaType(avail == 0) = "RESERVA";

    % 2) Slots ocupats pel programa (assignats) -> PROGRAM (això té prioritat)
    miniReservaType(ocupat == 1) = "PROGRAM";

    % 3) Slots lliures reals (equivalent a freeMask global)
    miniAvail = (avail == 1) & (ocupat == 0);
    miniReservaType(miniAvail) = "NOT_USED";

    nMiniAvailDisponible = double(sum(miniAvail(:)));

    Res = Res(1:fila, :);

    %% 13. Mostrar taula i estadístiques

    % Mostrar la taula de planificació
    disp(Res);

    % Estadístiques globals dels indicadors
    totalAssignats   = sum(assignat);
    totalNoAssignats = nP - totalAssignats;

    %% 14. Visualització simple per dia (matriu sala x temps)
    for dia = 1:nDies
        tRange = slotRangeDia(dia);
        schedule = zeros(nR, slotsDia);
        for k = 1:nR
            for t = tRange
                for i = 1:nP
                    if x_opt(i,k,t) > 0.5
                        schedule(k, t - tRange(1) + 1) = i;
                    end
                end
            end
        end
        fprintf('Matriu DIA %d (files = sales, columnes = slots de %d min):\n', ...
            dia, Delta);
        disp(schedule);
    end

    %% 15. Planificació detallada per pantalla

    % Ordenem la taula per FECHA i HORA
    Res = sortrows(Res, {'FECHA','HORA','SALA'});

    fprintf('\n========= PLANIFICACIÓ DETALLADA =========\n');
    for f = 1:height(Res)
        fprintf('ID %4g | Sala %d | %s | %s - %s | %-12s | %s\n', ...
            Res.ID_PROC(f), ...
            Res.SALA(f), ...
            string(Res.FECHA(f), 'yyyy-MM-dd'), ...
            char(Res.HORA(f)), ...
            char(Res.HORA_FI(f)), ...
            Res.TIPUS(f), ...
            Res.PRIORITAT(f));
    end
    fprintf('==========================================\n\n');

    fprintf('\n================ RESUM GLOBAL ================\n');
    fprintf('Dades problema: Total Files In: %d, Total proc: %d Total', ...
        Total_Files, nP_total, nUH_total);
    fprintf('Total procediments           : %d\n', nP);
    fprintf('Assignats a algun dia        : %d\n', totalAssignats);
    fprintf('NO assignats                 : %d\n', totalNoAssignats);
    fprintf('Percent NO assignats         : %.1f%%%%\n', 100 * totalNoAssignats / nP);

    fprintf('\n--- Detall per dia ---\n');
    for dia = 1:nDies
        slotsTotalsDia = nR * slotsDia;
        slotsReservatsDia = slotsTotalsDia - slotsDisponiblesPerDia(dia);

        fprintf(['Dia %d: %2d procediments assignats, slots lliures = %3d de %3d ' ...
            '(disponibles=%3d, reservats=%3d)\n'], ...
            dia, assignatsPerDia(dia), slotsLliuresPerDia(dia), ...
            slotsTotalsDia, slotsDisponiblesPerDia(dia), slotsReservatsDia);
    end

    fprintf('==============================================\n\n');

    fprintf('=== Estadístiques per sala (sobre tots els dies) ===\n');
    for k = 1:nR
        fprintf('Sala %d: ús = %.1f min, oci = %.1f min, delta = %.1f (mu=%.1f)\n', ...
            k, u_opt(k), idle_opt(k), delta_opt(k), mu);
    end

    %% 16. Estadística extra: diferència entre hospitalitzacions/urgències reals i reserves fixades

    % Comptem hospitalitzacions i urgències reals (les descartades del
    % model): nHosp, nUrg
    temps_Reserva_Acumulat = temps_Reserva_Acumulat + temps_Reserva;

    %% Bloc 17. Persistència de resultats

    %% 17a. Guardar la taula mini avail en un fitxer
    Tmini = table( ...
        miniAvail(:), ...
        miniReservaType(:), ...
        'VariableNames', {'MINI_AVAIL', 'MINI_RESERVA_TYPE'} );
    nomFitxerParcial = string(nomDirResults) + filesep + "miniAvail_debug" + string(inP_ini) + ".xlsx";
    writetable(Tmini, nomFitxerParcial);
    fprintf('Fitxer CSV generat correctament: %s\n', nomFitxerParcial);

    %% 17b. Guardar la taula Res en un fitxer CSV (parcial)
    nomFitxerParcial = string(nomDirResults) + filesep + "resultat_planificacio" + string(inP_ini) + ".xlsx";
    writetable(Res, nomFitxerParcial);
    fprintf('Fitxer CSV generat correctament: %s\n', nomFitxerParcial);

    %% 17c. Guardar la taula Res en un fitxer CSV (afegint-hi resultats)
    nomFitxerTotal = string(nomDirResults) + filesep + "resultat_planificacio_total.xlsx";
    if isfile(nomFitxerTotal)
        % El fitxer ja existeix → afegim files SENSE repetir capçaleres
        writetable(Res, nomFitxerTotal, ...
            'WriteMode', 'append', ...
            'WriteVariableNames', false);
    else
        % Primer cop → creem el fitxer amb capçaleres
        writetable(Res, nomFitxerTotal);
    end

    fprintf('Resultats afegits al fitxer CSV: %s\n', nomFitxerTotal);


end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------ LOCAL FUNCTIONS ------------------------ %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Res, fila] = nextRowRes(Res, fila, growBlock)
%NEXTROWRES Incrementa fila i amplia Res si cal (creixement per blocs).
%
%   [Res, fila] = nextRowRes(Res, fila, growBlock)

    fila = fila + 1;

    if fila > height(Res)
        extra = growBlock;

        % Ampliar duplicant una "plantilla" de files
        Res = [Res; Res(1:extra,:)];

        % Netejar les noves files (opcional però recomanat)
        idxNew = (height(Res)-extra+1) : height(Res);

        Res.ID_PROC(idxNew)    = NaN;
        Res.SALA(idxNew)       = NaN;
        Res.FECHA(idxNew)      = NaT;
        Res.HORA(idxNew)       = duration(0,0,0);
        Res.HORA_FI(idxNew)    = duration(0,0,0);
        Res.TIPUS(idxNew)      = "";
        Res.PRIORITAT(idxNew)  = "";
        Res.ID(idxNew)         = NaN;
        Res.FECHA_ORG(idxNew)  = NaT;
        Res.HORA_ORG(idxNew)   = duration(0,0,0);
        Res.DUR_RES_MIN(idxNew)= duration(0,0,0);
        Res.DUR_IDLE_MIN(idxNew)=duration(0,0,0);
        Res.DUR_PROC(idxNew)   = duration(0,0,0);
    end
end

