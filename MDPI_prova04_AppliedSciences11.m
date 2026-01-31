%% MILP per a planificació d'endoscopies (versió simple)
% - Variables x_{i,k,t} (un procediment i comença a sala k, slot temporal t)
% - Indexació "natural MATLAB": idx(i,k,t) = i + (k-1)*nP + (t-1)*nP*nR
% - Malla temporal de 15 minuts (Delta)
% - Objectiu: prioritzar puntualitat (minimitzar hora de començament)
% - Slots reservats (urgents / hospitalitzats) com a temps bloquejat

% Amb el model anterior teníem:

% #1. Només variables x_{i,k,t} (és a dir, "només decideixo on i quan comença cada procediment").
% #2. No solapament garantit per sala: un cop un procediment comença, el model compta tots els slots que ocupa i no deixa que un altre s'hi encavalqui.
% #3. Reserves (URG/HOSP) modelades com a avail(k,t) = 0 i ub(idx(i,k,t)) = 0, de manera que cap procediment "REG" pot ocupar aquests slots.
% #4. Puntualitat: la funció objectiu minimitza la suma dels temps de començament (t-1)*delta, així que el solver ...
%     intentarà posar els procediments els més aviat possible, respectant bloquejos i no solapament.

% Ara afegim:

% Afegim variables contínues:

% #7. u_k: minuts d'ús total de cada sala 
% #8. idle_k: minuts d'oci total de cada sala 
% #9. delta_k: desviació d'ús respecte d'un objectiu mu (balanceig entre sales)

% 4. Nova funció objectiu (multi-terme):
% amb pesos w_start, w_idle, w_bal (pots jugar-hi).

% selecció de nP per un valor de bloc inP (v2)

% v3: iteració per blocs
% Aborta si el temps disponible és menor que el temps total dels
% procediments

clear; clc;
close all;

%% VARIABLES DE LOG I RESULTATS
tClockStart = datetime('now');
tMain = tic;

sClockStart =  char(tClockStart, 'HH:mm:ss');
disp(['Start: ' sClockStart]);

%% 1. Paràmetres bàsics del problema
little=false;
if little
    num_days_desired = 1;            % nombre de dies
    default_inP = 20;
    num_rooms = 2;
else
    num_days_desired = 5;            % nombre de dies
    default_inP = 70;
    num_rooms = 4;
end


% Using datetime (recommended)
ts = datetime('now','Format','yyyyMMdd_HHmmss');   % e.g. 20251227_143025
logFilename = "myLog_" + string(ts) + ".txt";
logDirname = "myLog_" + string(ts);
logProblem = "_B"+string(num_days_desired)+"inP"+string(default_inP)+"NR"+string(num_rooms);
% logFilename is a string: "myLog_20251227_143025.txt"

% Crear directori de resultats (si no existeix)
nomDirResults = "Results_" + string(logDirname)+string(logProblem);

if ~isfolder(nomDirResults)
    mkdir(nomDirResults);
end

% Camí complet del fitxer de log
logPath = fullfile(nomDirResults, logFilename);

% Activar diary
diary(logPath);
diary on




%% CARREGA FITXER ENTRADA

if little
    filename = 'Dades_COL_GAS_LittleOK.mat';
else
    filename = 'Dades_COL_GAS.mat';
end

% Comprovem si el fitxer 'Dades_COL_GAS.mat' existeix abans de carregar-lo
if isfile(filename)   
    S = load(filename);   % carrega en una struct
    noms = fieldnames(S);
    T = S.(noms{1});      % assigna el contingut a T
else  
    error('El fitxer %s no existeix.', filename);
end

% Ordenem per data de sol·licitud
T = sortrows(T, {'FECHA', 'HORA'});

if ~ismember('ID', T.Properties.VariableNames)
    T.ID = (1:height(T)).';
end
T = movevars(T, 'ID', 'Before', 1);

% Si little=true, forcem prioritat "Urgente" als registres indicats
% Afegim urgents i hospitalitzados a efectes de debug
if little
    idx = [5 8 12];
    % Evitem índexs fora de rang
    idx = idx(idx >= 1 & idx <= height(T));
    T.Prioridad(idx) = "Urgente";

    idx = [3 6];
    % Evitem índexs fora de rang
    idx = idx(idx >= 1 & idx <= height(T));
    T.Prioridad(idx) = "Hospitalizado";
end
% END CARREGA DE DADES

%% Declaració lastMilp, última volta del bucle
lastMilp = false;
missatgeLastMilp = ['Set LastMilp = true, el valor que estem agafant de procesos\n' ...
                'en aquest bucle, deixarà el total de procediments tractats.' ...
                'nP_Total %d < inP_ini %d + inP %d\n'];

%% Info bàsica de les dades
Total_Dies = numel(unique(T.FECHA)); % Nombre de dies recollits
Total_Files = size(T,1);



%% PARÀMETRES SIMULACIÓ
excluded_dates = datetime.empty(0,1);
% excluded_dates = ["2026-01-06", "2026-01-01"];
start_date = "2024-02-01";
Delta      = 15;                    % minuts per slot
hores_torn = 10;                    % hores per dia (p.ex. 8:00-16:00)
slotsDia   = hores_torn*60/Delta;   % slots per dia (8h -> 32 slots)
mins_colono=60;
mins_gastro=30;
mins_other=60;
offset_hour = hours(9);   % offset de +9 hores
miniMILP = false;
Tdebug = [];
TdebugMILP = [];


[dies_seleccio, nDies, next_start_date] = generaDiesLaborables( ...
    start_date, num_days_desired, excluded_dates);
[nT, nR, nSlots, nMinsTot] = calculaDimensions(nDies, slotsDia, Delta, num_rooms);

ini_date_non_prog_proc_filter = start_date;
last_date_non_prog_proc_filter = next_start_date - days(1);

%% Impressió de dades bàsiques del problema
fprintf('\n=============================================\n');
fprintf('         INFORMACIÓ BÀSICA DEL PROBLEMA      \n');
fprintf('=============================================\n');

% Dades d'origen
fprintf('Total de dies al fitxer (FECHA):     %d\n', Total_Dies);
fprintf('Total de procediments al fitxer (FECHA):     %d\n', Total_Files);

% Paràmetres temporals
fprintf('\n--- Paràmetres temporals ---\n');
fprintf('Delta (minuts per slot):             %d min\n', Delta);
fprintf('Hores per torn(hores_torn):          %d h\n', hores_torn);
fprintf('Slots per dia per sala(slotsDia):    %d\n', slotsDia);
fprintf(['Nombre de dies A PLANIFICAR per defecte en cada bloc(num_days_desired):  %d\n'], num_days_desired);
fprintf('Total slots disponibles nDies*slotsDia (nT):         %d\n', nT);

% Recursos
fprintf('\n--- Recursos ---\n');
fprintf('Nombre de sales (nR):                 %d\n', nR);
fprintf('Nombre slots total = slotsDia*nDies(nT)*nSales(nR) %d\n', nSlots);
fprintf('Total Minuts disponibles SENSE RESERVES %d\n', nMinsTot);
fprintf('\n');
fprintf('\n=============================================\n\n');


%% Valors de Prioridad que vols conservar per a la simulació
valsPrioridad = ["Screening"  "Preferente"  "Regular"]; 

% "Hospitalizado" i "Urgente" ja tenen els seus espais
valsDescarte = ["Urgente" "Hospitalizado"];

% Màscara lògica sobre el camp T.Prioridad (sense filtrar encara per dies)
maskPrioridadSim = ismember(T.Prioridad, valsPrioridad);
maskDescarteSim = ismember(T.Prioridad, valsDescarte);

% Subtaula amb només aquestes prioritats (tots els dies)
Tprio = T(maskPrioridadSim, :);

% Subtaula amb només descartats (tots els dies)
Tdesc = T(maskDescarteSim, :);

% Guardar la taula Tdesc en un fitxer CSV 
nomFitxer = nomDirResults + filesep + "descartats" + ".xlsx";
writetable(Tdesc, nomFitxer);
fprintf('Fitxer CSV Tdesc generat correctament: %s\n', nomFitxer);

% Nombre total de procediments candidats abans de filtrar per dies
nP_total = size(Tprio, 1);
nUH_total = size (Tdesc, 1);

fprintf('\n=============================================\n');
fprintf('     RESULTATS DE SELECCIÓ DE PRIORITATS \n');
fprintf('=============================================\n\n');

%% Resum inicial abans de filtrar per dies
fprintf('--- Resum inicial (Tprio, només prioritats de simulació) ---\n');
fprintf('Total de registres %d\n', Total_Files);
fprintf('Total de registres amb prioritats [%s]: %d\n', strjoin(valsPrioridad, ', '), nP_total);
fprintf('Total de registres Descartats (urgent i hospitalitzat): %d\n', nUH_total);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% b) selecció per numero de procediments
% Ara, SELECCIÓ PER NOMBRE DE PROCEDIMENTS (inP) SOBRE Tprio
inP = default_inP;   % Nombre de procediments que vols simular
inP_ini=1;
remanent=0;

if(inP_ini+inP > nP_total)
    inP=nP_total - inP_ini;
end


%% SELECCIÓ BLOCS DE DATES DE PROGRAMACIÓ
if ~isempty(dies_seleccio)
    fmt = 'yyyy-MM-dd';
    fprintf( ...
        'PROGRAMACIO BLOC INICIAL \nDia Inici: %s, dia Fi: %s, total dies bloc %d, nDies %d\n', ...
        string(dies_seleccio(1), fmt), ...
        string(dies_seleccio(end), fmt), ...
        numel(dies_seleccio), nDies);
else
    fprintf('Bloc sense dies seleccionats\n');
end


% Slots disponibles i reserves (URG/HOSP)
[avail, reserveType] = buildAvailReserves(nDies, nR, nT, slotsDia, dies_seleccio);

% Helper per calcular rangs de t dins de cada dia
slotRangeDia = @(dia) ((dia-1)*slotsDia + 1) : (dia*slotsDia);

temps_Disponible_4_sales=sum(sum(avail))*Delta;
temps_Reserva = nMinsTot - temps_Disponible_4_sales;
slots_Reserva = temps_Reserva / Delta;
temps_Reserva_Acumulat = 0; % inicialitzem contador per acumular temps de reserva dels bucles
fprintf ('Amb les condicions dissenyades entre les nR %d sales disposen de %d minuts TEÒRICS \n', nR, nMinsTot);
fprintf ('Amb les condicions dissenyades entre les nR %d sales disposen de %d minuts REALS\n', nR, temps_Disponible_4_sales);
fprintf ('Suposa una reserva de %d minuts, i %d slots reservats en les nR %d sales\n', temps_Reserva, slots_Reserva, nR);

%% INICIALITZEM BLOC BUCLE PRINCIPAL
bucle=1;
first_iter = false;
lastID = min(T.ID);

tmp = temps_Disponible_4_sales/mins_colono;
if (default_inP>tmp)
    fprintf (['Correcció default_inP de usuari, per ajustar a condicions de sala: ' ...
        'default_inP usuari %d, default_inP:%d\n'], default_inP, tmp);
    default_inP = tmp;
end


while inP > 0

    % Reinicialitzem estat del bloc de selecció
    first_iter = true;
    remanent = 0;
    tBloc = tic;

    while true
        %% Dies seleccionats
        % Ens assegurem de no demanar més del que hi ha
        % Taula final de simulació: només prioritats desitjades i els primers nP procediments
        %T2p = Tprio(inP_ini:inP_ini+nP-1, :);
        % Mantén el teu càlcul (no el toquem)

        nP = min(inP, size(Tprio,1));

        % Creació robusta de T2p
        [T2p, idx_fi] = buildProceduresTableToProcess(Tprio, inP_ini, nP, nP_total);
        debugPrintSelecProcediments(inP, nP, inP_ini, idx_fi, bucle);
        dies_T2p = getInputDataDaysRange (T2p);

        %% Construcció del vector de durades d (últim pas)
        [d, temps_Necessari, temps_Disponible_4_sales, remanent, remanent_prev, tipus] = ...
            compute_d_ProcedureDurationsAndRemainingTime( ...
            T2p, nP,mins_gastro, mins_colono, avail, Delta, remanent);

        % impressió temps disponible vegeu: avail
        debugPrintTempsDisponibilitat(temps_Necessari, temps_Disponible_4_sales, ...
            remanent, nR);

        % SI ENS PASSEM DE TEMPS → TORNAR A L'ESTAT ANTERIOR  --> Ojo que
        if remanent < 0
            if ~first_iter
                % Recuperem l'estat del pas anterior
                nP        = nP_prev;
                inP       = inP_prev;
                T2p       = T2p_prev;
                d         = d_prev;
                tipus     = tipus_prev;
                remanent  = remanent_prev;
                fprintf('Recuperem estat anterior \n');
            end
            break;
        end

        %% I provem amb un inP més gran
        if(inP_ini + inP > nP_total)
            fprintf(missatgeLastMilp, nP_total, inP_ini, inP);
            fprintf('\iniciem MILP amb aquest últim slot\n');
            lastMilp = true;
            Tdebug = afegeixDebugBloc(Tdebug, bucle, nP, inP_ini, nT, nR, ...
                nSlots, nMinsTot, tipus, nDies, dies_seleccio, nHosp, nUrg, ...
                dies_T2p,lastID, maxID, nP_total, Total_Files, Total_Dies, ...
                lastMilp, default_inP, "set LastMilp true: inP_ini + inP > nP_total");
            break;
            % surt del bucle pq ja ha trobat el màxim valor de procediments a tractar
        else
            inP = inP + 1;
        end

        %% Si aquest pas és vàlid, el guardem com a "pas anterior"
        nP_prev       = nP;
        inP_prev      = inP;
        T2p_prev      = T2p;
        d_prev        = d;
        tipus_prev    = tipus;
        remanent_prev = remanent;
        first_iter     = false;


    end % while


    %% BLOC : Recuperem els procediments descartats
    maxID = max(T2p.ID);
    T2d = Tdesc(Tdesc.ID >= lastID & Tdesc.ID <= maxID, :);
    nHosp = sum(strcmp(T2d.Prioridad, 'Hospitalizado'));
    nUrg  = sum(strcmp(T2d.Prioridad, 'Urgente'));

    % imprimim condicions de crida del MILP
    Tdebug = debugPrintFileMILPConditions(bucle, nP, inP_ini, nT, nR, nSlots, ...
        nMinsTot, tipus, nDies, dies_seleccio, nHosp, nUrg, T2p, Tdebug, ...
        lastID, maxID, nP_total, Total_Files, Total_Dies, ...
        lastMilp, default_inP);

    fprintf('\n------------ CRIDEM FUNCIÓ MILP -----------------------\n');
    fprintf('BUCLE %d\n\n', bucle);
    % Variables decidides pel MILP → x_opt, u_opt, idle_opt, delta_opt, mu
    % Variables de context → avail, Delta, d, nP, nR, nT
    % Variables d'estat → exitflag, output, fval
    % Si el proves i veus que:
    % (1) concentra massa en una sala → puja w_bal
    % (2) t'omple massa tard el dia → puja w_start
    % (3) et deixa molt oci perquè no hi ha prou demanda → baixa w_idle o fes mu més petit

    sol = executarMILP(nP, nR, nT, avail, Delta, d);

    if sol.exitflag ~= 1
        fprintf('Sortint del bucle: MILP no resolt correctament (exitflag=%d)\n', sol.exitflag);
        break   % <-- ARA sí, perquè estàs dins del while
    end
    fprintf('------------ SORTIDA FUNCIÓ MILP -----------------------\n');
    TdebugMILP = afegeixMILPDebug(TdebugMILP, bucle, nP, nR, nT, sol.fval, ...
        sol.exitflag, sol.output, sol.elapsed);
    %% recuperar la sortida del MILP
    % --- Desempaquetar "sol" a variables locals (com abans del refactor) -----
    [ exitflag, output, fval, ...
        x_opt, u_opt, idle_opt, delta_opt, mu, ...
        nP, nR, nT, Delta, d, avail, ...
        nAssignats, theta ] = unpackSolMiniMILP(sol);


    %% 9. Mostrar l'horari resultant
    %nAssignat = 0;
    [nT_real, nDies, slotRangeDia] = fixTimeInconsistency( ...
        nT, x_opt, slotsDia, nDies, dies_seleccio, slotRangeDia);

    [Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat,...
        temps_Reserva, avail_mini, reserveType_mini, nAvailDisponible_mini] = ...
        generaIMostraResultatsPlanificacio( ...
        offset_hour, nDies, slotRangeDia, Delta, slotsDia, dies_seleccio, ...
        nP, nR, d, tipus, T2p, avail, reserveType, ...
        x_opt, u_opt, idle_opt, delta_opt, mu, ...
        Total_Files, nP_total, nUH_total, nomDirResults, inP_ini, ...
        last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva);

    if miniMILP
        if nAvailDisponible_mini > 1 % qualsevol de les proves que executem com a mínim són 2 slots
            tempStatusMainMilp = packMiniMilpRunState( ...
                exitflag, output, fval, x_opt, u_opt, idle_opt, delta_opt, mu, ...
                nP, nR, nT, Delta, d, avail, nAssignats, theta, ...
                nT_real, nDies, slotRangeDia, ...
                Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva, ...
                reserveType, 0, ...
                slotsDia, dies_seleccio, offset_hour, tipus, T2p, reserveType, ...
                Total_Files, nP_total, nUH_total, nomDirResults, inP_ini);
    
            % minimilpbucle
            default_inP_mini = 1;
            lastMilp_mini = false;
            first_iter_mini = true;
            bucle_mini=1;
            % %T2p_mini, avail_mini, reserveType_mini, nAvailDisponible_mini
            % TO DEBUG WITH DATA MINI
            % Tprio.Tipo_Prueba(22,1) = "Gastroscopia";
            % Tprio.Tipo_Prueba(23,1) = "Gastroscopia";
            % avail_mini(1,39) = 1;
            % avail_mini(1,40) = 1;
            % reserveType_mini(1,39) = "NOT_USED";
            % reserveType_mini(1,40) = "NOT_USED";
    
    
            while ~lastMilp_mini
                if first_iter_mini
                    [inP_ini, inP] = actualitzaIniP_ini(inP_ini, nP, default_inP_mini, nP_total);
                    iniP_ini_mini = inP_ini;
                else
                    first_iter_mini = false;
                    if(inP_ini+bucle> nP_total) %he consumit tots els procediments
                        fprintf('Ajusto tamany bloc als procediments restant %d \n', inP);
                        inP=nP_total - inP_ini;
                    end
                end
                % Creació robusta de T2p
                [T2p_mini, idx_fi] = buildProceduresTableToProcess(Tprio, inP_ini, bucle_mini, nP_total);
                debugPrintSelecProcediments(inP, nP, inP_ini, idx_fi, bucle_mini);
                
                [d, temps_Necessari, temps_Disponible_4_sales, remanent, remanent_prev, tipus] = ...
                compute_d_ProcedureDurationsAndRemainingTime(T2p_mini, nP ,mins_gastro, ...
                mins_colono, avail_mini, Delta, remanent);
    
                debugPrintTempsDisponibilitat(temps_Necessari, temps_Disponible_4_sales, ...
                    remanent, nR);
    
                Tdebug = debugPrintFileMILPConditions(bucle_mini, bucle_mini, inP_ini, nT, nR, nSlots, ...
                    nMinsTot, tipus, nDies, dies_seleccio, 0, 0, T2p_mini, Tdebug, ...
                    lastID, maxID, nP_total, Total_Files, Total_Dies, ...
                    lastMilp_mini, default_inP_mini);
    
                sol = executarMILP(nP, nR, nT, avail_mini, Delta, d);
    
                if sol.exitflag ~= 1
                    fprintf('Sortint del bucle: miniMILP no resolt correctament (exitflag=%d)\n', sol.exitflag);
                    lastMilp_mini = true;
                end
                if lastMilp_mini %genero resultats quan ha acabat el bucle
                    if ~first_iter_mini
                        %recuperem pas anterior
                        [exitflag, output, fval, x_opt, u_opt, idle_opt, delta_opt, mu, ...
                            nP, nR, nT, Delta, d, avail_mini, nAssignats, theta, ...
                            nT_real, nDies, slotRangeDia, ...
                            Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva, ...
                            reserveType_mini, nAvailDisponible_mini, ...
                            slotsDia, dies_seleccio, offset_hour, tipus, T2p_mini, reserveType_mini, ...
                            Total_Files, nP_total, nUH_total, nomDirResults, inP_ini] = ...
                            unpackMiniMilpRunState(tempStatus);
                        fprintf('\n HEM POGUT ASSIGNAR AMB EL MILP %d\n', bucle_mini-1);
        
                        [Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat,...
                            temps_Reserva, avail_mini, reserveType_mini, nAvailDisponible_mini] = ...
                            generaIMostraResultatsPlanificacio( ...
                            offset_hour, nDies, slotRangeDia, Delta, slotsDia, dies_seleccio, ...
                            bucle_mini-1, nR, d, tipus, T2p_mini, avail_mini, reserveType_mini, ...
                            x_opt, u_opt, idle_opt, delta_opt, mu, ...
                            Total_Files, nP_total, nUH_total, nomDirResults, inP_ini, ...
                            last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva);
                    else
                        fprintf('\n hem intentat fer miniMILP sense èxit \n');
                    end
                    % recuperem estat main principal
                    [exitflag, output, fval, x_opt, u_opt, idle_opt, delta_opt, mu, ...
                        nP, nR, nT, Delta, d, avail, nAssignats, theta, ...
                        nT_real, nDies, slotRangeDia, ...
                        Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva, ...
                        reserveType, nAvailDisponible, ...
                        slotsDia, dies_seleccio, offset_hour, tipus, T2p, reserveType, ...
                        Total_Files, nP_total, nUH_total, nomDirResults, inP_ini] = unpackMiniMilpRunState(tempStatusMainMilp);
                    % com que hem processat files, hem d'actualitzar el
                    % punt d'inici per revisar la matriu
                    [inP_ini, inP] = actualitzaIniP_ini(iniP_ini_mini, bucle_mini-1, default_inP, nP_total);                
    
                else
                    [exitflag, output, fval,x_opt, u_opt, idle_opt, delta_opt, mu, ...
                        nP, nR, nT, Delta, d, avail_mini, nAssignats, theta ] = unpackSolMiniMILP(sol);
    
                    [nT_real, nDies, slotRangeDia] = fixTimeInconsistency( ...
                        nT, x_opt, slotsDia, nDies, dies_seleccio, slotRangeDia);
    
                    tempStatus = packMiniMilpRunState( ...
                        exitflag, output, fval, x_opt, u_opt, idle_opt, delta_opt, mu, ...
                        nP, nR, nT, Delta, d, avail_mini, nAssignats, theta, ...
                        nT_real, nDies, slotRangeDia, ...
                        Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva, ...
                        reserveType_mini, nAvailDisponible_mini, ...
                        slotsDia, dies_seleccio, offset_hour, tipus, T2p_mini, reserveType_mini, ...
                        Total_Files, nP_total, nUH_total, nomDirResults, inP_ini);
                end
                bucle_mini = bucle_mini + 1; 
                first_iter_mini = false;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% BLOC 18 PREPARACIÓ SEGÜENT BUCLE

    fprintf('----------- PREPARACIO SEGUENT BULCE---------------\n');
    fprintf('En el bucle ACTUAL %d:\nPer un total de procediments (np) %d \n', bucle, nP);
    fprintf('Iniciant en el procediment %d \n', inP_ini); 
    fprintf('Finalitzant en el procediment %d \n', inP_ini+nP-1);
    fprintf('nP_total %d \n', nP_total);
 
    if lastMilp 
        fprintf('\nTotal de registres processats, lastMilp TRUE, finalitzem procés\n\n');
        break;
    end

    bucle=bucle+1;
    first_iter = false;
    lastID = maxID;

    [dies_seleccio, nDies, next_start_date] = generaDiesLaborables( ...
        next_start_date, num_days_desired, excluded_dates);
    [nT, nR, nSlots, nMinsTot] = calculaDimensions(nDies, slotsDia, Delta, num_rooms);

    % Muntem les reserves
    [avail, reserveType] = buildAvailReserves(nDies, nR, nT, slotsDia, dies_seleccio);
    [inP_ini, inP] = actualitzaIniP_ini(inP_ini, nP, default_inP, nP_total);


    fprintf('---------------------------------------------\n---------------------------------------------\n');
    fprintf('En el bucle següent %d:\nPer un total de procediments inicials (inp) %d \n', bucle, inP);
    fprintf('Iniciant en el procediment %d \n', inP_ini); 
    fprintf('nP_total %d \n', nP_total);
    if ~isempty(dies_seleccio)
        fmt = 'yyyy-MM-dd';
        fprintf( ...
            'Dia Inici: %s, dia Fi: %s, total dies bloc %d, nDies %d\n', ...
            string(dies_seleccio(1), fmt), ...
            string(dies_seleccio(end), fmt), ...
            numel(dies_seleccio), nDies);
    else
        fprintf('Bloc sense dies seleccionats\n');
    end

    nDiesCheck = numel(dies_seleccio);

    if nDiesCheck ~= nDies
        warning( ...
            'generaDiesLaborables:nDiesMismatch', ...
            ['nDies (%d) diferent de dies realment seleccionats (%d). ' ...
            'Possibles exclusions dins del rang.'], ...
            nDies, nDiesCheck);
    end


    fprintf('---------------------------------------------\n---------------------------------------------\n');
   
    elapsed = toc(tBloc);
    elapsed_hms = seconds(elapsed);
    fprintf('Temps transcorregut pel bloc: %s\n', string(elapsed_hms, 'hh:mm:ss'));
    fprintf('------------------FI BLOC--------------------\n---------------------------------------------\n');
    fprintf('---------------------------------------------\n---------------------------------------------\n');
    fprintf('---------------------------------------------\n---------------------------------------------\n');

end % While bucle principal

Tdebug = afegeixDebugBloc( Tdebug, bucle,nP, inP_ini, nT, nR, nSlots, nMinsTot, ...
    tipus, nDies, dies_seleccio, nHosp, nUrg, ...
    dies_T2p,lastID, maxID, nP_total, Total_Files, Total_Dies, ...
    lastMilp, default_inP, "sortit de MILP");

%% Construcció del vector de durades d (últim pas)
tipus_desc = Tdesc.Tipo_Prueba;  % "Gastroscopia" o "Colonoscopia" o "Otros"
nUH_bucle = size (T2d, 1);
dDesc = zeros(1, nUH_bucle);         % durada en minuts
dDesc(tipus=="Gastroscopia")  = 30;
dDesc(tipus=="Colonoscopia")  = 60;

temps_Necessari_HU=sum(dDesc);
temps_Reserva_deficit = temps_Reserva_Acumulat - temps_Necessari_HU;
fprintf('================ RESERVES VS DEMANDA REAL ================\n');
fprintf('Per realitzar totes les intervencions %d H/U es necessiten %d minuts\n', height(Tdesc), temps_Necessari_HU);
fprintf('Es disposa de reserva de %d minuts reservats per bloc,\n que amb un total de %d blocs,\n suposa un temps acumulat de %d minuts reservats,\n i un dèficit de %d minuts.\n', ...
    temps_Reserva, bucle, temps_Reserva_Acumulat, temps_Reserva_deficit);

disp('ATENCIÓ PER VALIDAR EL TEMPS DE LES URGÈNCIES I HOSPITALITZACIONS CALDRIA VALIDAR QUE LA DATA DE REALITZACIÓ CAU DINS LA DATA DE PLANIFICACIÓ');
ini_date_non_prog_proc_filter = datetime(ini_date_non_prog_proc_filter);
maxFecha = max(Tdesc.FECHA);
last_date_non_prog_proc_filter  = next_start_date - days(1);
if last_date_non_prog_proc_filter > maxFecha
    last_date_non_prog_proc_filter  = maxFecha;
end


TdescInData = Tdesc(Tdesc.FECHA >= ini_date_non_prog_proc_filter & Tdesc.FECHA <= last_date_non_prog_proc_filter, :);
fprintf('Data inici programació: %s\n', string(ini_date_non_prog_proc_filter, 'yyyy-MM-dd hh:mm:ss'));
fprintf('Data última programació: %s\n', string(last_date_non_prog_proc_filter, 'yyyy-MM-dd hh:mm:ss'));

nomFitxer = nomDirResults + filesep + "descartats_filtre_data" + ".xlsx";
idx = Tdesc.FECHA >= ini_date_non_prog_proc_filter & Tdesc.FECHA <= last_date_non_prog_proc_filter;
writetable(Tdesc(idx, :), nomFitxer);
fprintf('Fitxer CSV Tdesc generat correctament: %s\n', nomFitxer);

nomFitxer = nomDirResults + filesep + "tauladebug" + ".xlsx";
writetable(Tdebug, nomFitxer);
fprintf('Fitxer CSV Tdebug generat correctament: %s\n', nomFitxer);

nomFitxer = 'tauladebugMILP';
persistMILPDebug(TdebugMILP, nomDirResults, nomFitxer);
fprintf('Fitxers TdebugMILP generat correctament: %s\n', nomFitxer);

nomFitxer = 'reserves_vs_demanda.xlsx';
persistReservesVsDemandaCSV( ...
    nomDirResults, ...
    nomFitxer, ...
    bucle, ...
    TdescInData, ...
    temps_Reserva, ...
    temps_Reserva_Acumulat);
fprintf('Fitxers ReservesVsDemanda generat correctament: %s\n', nomFitxer);
fprintf('==========================================================\n');


elapsed=toc(tMain);
elapsed_hms = seconds(elapsed);

fprintf('Temps transcorregut MILP: %s\n', datestr(elapsed_hms, 'HH:MM:SS'));


% - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [inP_ini, inP] = actualitzaIniP_ini(inP_ini, nP, default_inP, nP_total)

    inP_ini = inP_ini+nP; % actualizo el valor del procediment inicial pq segueixi després de l'últim procediment processat
    inP=default_inP; %reinicialitzo el valor de la temptativa de procediments al valor inicial per defecte
    if(inP_ini+inP> nP_total) %he consumit tots els procediments
        fprintf('Ajusto tamany bloc als procediments restant %d \n', inP);
        inP=nP_total - inP_ini;
    end
end


function [ ...
    exitflag, output, fval, ...
    x_opt, u_opt, idle_opt, delta_opt, mu, ...
    nP, nR, nT, Delta, d, avail, ...
    nAssignats, theta ] = unpackSolMiniMILP(sol)
%UNPACKSOLMINIMILP Desempaqueta l'struct sol del MiniMILP
% Retorna variables individuals (no struct)

    % --- Estat del solver ---
    exitflag = sol.exitflag;
    output   = sol.output;
    fval     = sol.fval;

    % --- Solució ---
    x_opt     = sol.x_opt;
    u_opt     = sol.u_opt;
    idle_opt  = sol.idle_opt;
    delta_opt = sol.delta_opt;
    mu        = sol.mu;

    % --- Metadades (opc.) ---
    nP    = getfield_or_empty(sol, 'nP');
    nR    = getfield_or_empty(sol, 'nR');
    nT    = getfield_or_empty(sol, 'nT');
    Delta = getfield_or_empty(sol, 'Delta');
    d     = getfield_or_empty(sol, 'd');
    avail = getfield_or_empty(sol, 'avail');

    % --- Derivats (opc.) ---
    nAssignats = getfield_or_empty(sol, 'nAssignats');
    theta      = getfield_or_empty(sol, 'theta');
end

% -------------------------------------------------------------------------
function v = getfield_or_empty(s, fname)
% Retorna [] si el camp no existeix
    if isfield(s, fname)
        v = s.(fname);
    else
        v = [];
    end
end


function [nT_real, nDies, slotRangeDia] = fixTimeInconsistency( ...
    nT, x_opt, slotsDia, nDies_real, dies_seleccio, slotRangeDia)
%FIXTIMEINCONSISTENCY Ajusta rangs temporals si nT no quadra amb x_opt
%
% Inputs:
%   nT            : nT esperat (metadada)
%   x_opt         : solució MILP (3a dimensió = temps)
%   slotsDia      : nombre de slots per dia
%   nDies_real    : dies real (si el tens calculat fora)
%   dies_seleccio : array de datetimes seleccionats
%   slotRangeDia  : handle actual (es retorna actualitzat si cal)
%
% Outputs:
%   nT_real       : nT efectiu (size(x_opt,3))
%   nDies         : dies calculats a partir de nT_real i slotsDia
%   slotRangeDia  : handle actualitzat per truncar a nT_real

    nT_real = nT;
    nDies   = ceil(nT / slotsDia);

    % Si falta slotRangeDia, el construïm per defecte
    if nargin < 6 || isempty(slotRangeDia)
        slotRangeDia = @(dia) ((dia-1)*slotsDia + 1) : min(dia*slotsDia, nT);
    end

    % --- Inconsistència ---
    nT_x = size(x_opt, 3);
    if nT ~= nT_x
        warning(['Inconsistència temporal: nT=%d però size(x_opt,3)=%d. ' ...
                 'Ajustant impressió a nT_real.'], nT, nT_x);

        nT_real = nT_x;
        slotRangeDia = @(dia) ((dia-1)*slotsDia + 1) : min(dia*slotsDia, nT_real);
        nDies = ceil(nT_real / slotsDia);

        fprintf('\nDies real (%d dies):\n i nDies (dies_selecció) %d, nDies %d \n', ...
            nDies_real, numel(dies_seleccio), nDies);
    end
end

function debugPrintSelecProcediments(inP, nP, inP_ini, idx_fi, bucle)
%DEBUGPRINTSELECPROCEDIMENTS Imprimeix el resum de selecció de procediments

    fprintf('\n--- Selecció per nombre de procediments ---\n');
    fprintf('Procediments demanats (inP): %d\n', inP);
    fprintf('Procediments finalment seleccionats (nP): %d\n', nP);
    fprintf('Inici selecció: %d\n', inP_ini);
    fprintf('Fi selecció: %d\n', idx_fi);
    fprintf('En el Bucle: %d\n', bucle);

end


function [T2p, idx_fi] = buildProceduresTableToProcess( ...
    Tprio, inP_ini, nP, nP_total)
%BUILDPROCEDURESTABLETOPROCESS
% Construeix la taula de procediments a processar (T2p) a partir de Tprio
%
% Inputs:
%   Tprio     : taula ordenada per prioritat
%   inP_ini   : índex inicial
%   nP        : nombre de procediments a seleccionar
%   nP_total  : nombre total de files a Tprio
%
% Outputs:
%   T2p       : subtaula amb els procediments a processar
%   idx_fi    : índex final real utilitzat ([] si fora límits)
% agafo els nP possibles de la taula de procediments ja filtrada
% Tprio, no considera els urgents/hospi.
% per construir la taula d'urgents/hospi que es descarten, ens cal
% agafar l'id del Tprio (idx_fi) i seleccionar tots els Tdes que
% tingun un id inferior però superior al bucle anterior


    % Cas fora de límits
    if inP_ini < 1 || inP_ini > nP_total
        T2p = Tprio([],:);   % taula buida amb mateixes columnes
        idx_fi = [];
        fprintf('\n--- FORA LIMITS TAULA ---\n');
        return;
    end

    % Cas normal
    idx_fi = min(inP_ini + nP - 1, nP_total);
    T2p    = Tprio(inP_ini:idx_fi, :);
end

function [d, temps_Necessari, temps_Disponible, remanent, remanent_prev, tipus] = ...
    compute_d_ProcedureDurationsAndRemainingTime( ...
        T2p, nP, mins_gastro, mins_colono,  avail, Delta, remanent)
%COMPUTEPROCEDUREDURATIONSANDREMAININGTIME
% Calcula durades per procediment i temps disponible / remanent
%
% Inputs:
%   T2p          : taula de procediments a processar
%   nP           : nombre de procediments
%   mins_gastro  : durada gastroscòpia (min)
%   mins_colono  : durada colonoscòpia (min)
%   avail        : matriu disponibilitat (nR x nT)
%   Delta        : minuts per slot
%   remanent     : remanent anterior
%
% Outputs:
%   d                 : vector durades per procediment (1 x nP)
%   temps_Necessari   : minuts totals necessaris
%   temps_Disponible  : minuts totals disponibles (totes les sales)
%   remanent          : nou remanent
%   remanent_prev     : remanent anterior

    % --- Tipus de procediment ---
    tipus = T2p.Tipo_Prueba;   % "Gastroscopia" | "Colonoscopia" | "Otros"

    % --- Durades ---
    d = zeros(1, nP);
    d(tipus == "Gastroscopia") = mins_gastro;
    d(tipus == "Colonoscopia") = mins_colono;
    % Els "Otros" queden a 0 (mateixa lògica que l'original)

    % --- Temps necessari ---
    temps_Necessari = sum(d);

    % --- Temps disponible ---
    temps_Disponible = sum(avail, 'all') * Delta;

    % --- Remanent ---
    if remanent == 0
        remanent = temps_Disponible;
    end

    remanent_prev = remanent;
    remanent      = temps_Disponible - temps_Necessari;
end


function debugPrintTempsDisponibilitat( ...
    temps_Necessari, ...
    temps_Disponible_4_sales, ...
    remanent, ...
    nR)
%DEBUGPRINTTEMPSDISPONIBILITAT
% Imprimeix el resum de temps necessari vs disponible

    disp(['Per realitzar totes les intervencions es necessiten ' ...
          num2str(temps_Necessari) ' minuts']);

    disp(['Amb les condicions dissenyades entre les 4 sales disposen de ' ...
          num2str(temps_Disponible_4_sales) ' minuts']);

    disp(['En els nDies acumulats hi ha una Diferència de : ' ...
          num2str(remanent) ' minuts']);

    disp(['Son ' ...
          num2str((temps_Disponible_4_sales - temps_Necessari) / nR) ...
          ' min per sala']);
end


function S = packMiniMilpRunState( ...
    exitflag, output, fval, x_opt, u_opt, idle_opt, delta_opt, mu, ...
    nP, nR, nT, Delta, d, avail_mini, nAssignats, theta, ...
    nT_real, nDies, slotRangeDia, ...
    Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva, ...
    reservaAvail_mini, nMiniAvailDisponible, ...
    slotsDia, dies_seleccio, offset_hour, tipus, T2p_mini, reserveType_mini, ...
    Total_Files, nP_total, nUH_total, nomDirResults, inP_ini)

%PACKMINIMILPRUNSTATE Empaqueta estat d'un run MiniMILP en una struct

    S = struct();

    % --- Solver / solució ---
    S.exitflag  = exitflag;
    S.output    = output;
    S.fval      = fval;
    S.x_opt     = x_opt;
    S.u_opt     = u_opt;
    S.idle_opt  = idle_opt;
    S.delta_opt = delta_opt;
    S.mu        = mu;

    % --- Metadades MILP ---
    S.nP    = nP;
    S.nR    = nR;
    S.nT    = nT;
    S.Delta = Delta;
    S.d     = d;

    % --- Disponibilitat mini + derivats ---
    S.avail_mini         = avail_mini;
    S.reservaAvail_mini  = reservaAvail_mini;
    S.nMiniAvailDisponible = nMiniAvailDisponible;

    % --- Camps opcionals del sol ---
    S.nAssignats = nAssignats;
    S.theta      = theta;

    % --- Temps/Calendari ajustat ---
    S.nT_real      = nT_real;
    S.nDies        = nDies;
    S.slotRangeDia = slotRangeDia;

    % --- Resultats planificació ---
    S.Res                         = Res;
    S.last_date_non_prog_proc_filter = last_date_non_prog_proc_filter;
    S.temps_Reserva_Acumulat      = temps_Reserva_Acumulat;
    S.temps_Reserva               = temps_Reserva;

    % --- Inputs rellevants del run (per reproduir/debug) ---
    S.slotsDia         = slotsDia;
    S.dies_seleccio    = dies_seleccio;
    S.offset_hour      = offset_hour;
    S.tipus            = tipus;
    S.T2p_mini         = T2p_mini;
    S.reserveType_mini = reserveType_mini;

    S.Total_Files   = Total_Files;
    S.nP_total      = nP_total;
    S.nUH_total     = nUH_total;
    S.nomDirResults = nomDirResults;
    S.inP_ini       = inP_ini;
end
function [ ...
    exitflag, output, fval, x_opt, u_opt, idle_opt, delta_opt, mu, ...
    nP, nR, nT, Delta, d, avail_mini, nAssignats, theta, ...
    nT_real, nDies, slotRangeDia, ...
    Res, last_date_non_prog_proc_filter, temps_Reserva_Acumulat, temps_Reserva, ...
    reservaAvail_mini, nMiniAvailDisponible, ...
    slotsDia, dies_seleccio, offset_hour, tipus, T2p_mini, reserveType_mini, ...
    Total_Files, nP_total, nUH_total, nomDirResults, inP_ini] = ...
    unpackMiniMilpRunState(S)

%UNPACKMINIMILPRUNSTATE Desempaqueta estat d'un run MiniMILP des d'una struct

    % --- Solver / solució ---
    exitflag  = S.exitflag;
    output    = S.output;
    fval      = S.fval;
    x_opt     = S.x_opt;
    u_opt     = S.u_opt;
    idle_opt  = S.idle_opt;
    delta_opt = S.delta_opt;
    mu        = S.mu;

    % --- Metadades MILP ---
    nP    = S.nP;
    nR    = S.nR;
    nT    = S.nT;
    Delta = S.Delta;
    d     = S.d;

    % --- Disponibilitat mini + derivats ---
    avail_mini          = S.avail_mini;
    reservaAvail_mini   = S.reservaAvail_mini;
    nMiniAvailDisponible = S.nMiniAvailDisponible;

    % --- Camps opcionals del sol ---
    nAssignats = S.nAssignats;
    theta      = S.theta;

    % --- Temps/Calendari ajustat ---
    nT_real      = S.nT_real;
    nDies        = S.nDies;
    slotRangeDia = S.slotRangeDia;

    % --- Resultats planificació ---
    Res                          = S.Res;
    last_date_non_prog_proc_filter = S.last_date_non_prog_proc_filter;
    temps_Reserva_Acumulat       = S.temps_Reserva_Acumulat;
    temps_Reserva                = S.temps_Reserva;

    % --- Inputs rellevants del run (per reproduir/debug) ---
    slotsDia         = S.slotsDia;
    dies_seleccio    = S.dies_seleccio;
    offset_hour      = S.offset_hour;
    tipus            = S.tipus;
    T2p_mini         = S.T2p_mini;
    reserveType_mini = S.reserveType_mini;

    Total_Files   = S.Total_Files;
    nP_total      = S.nP_total;
    nUH_total     = S.nUH_total;
    nomDirResults = S.nomDirResults;
    inP_ini       = S.inP_ini;
end
