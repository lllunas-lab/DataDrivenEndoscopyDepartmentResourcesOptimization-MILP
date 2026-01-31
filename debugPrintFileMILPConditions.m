function Tdebug = debugPrintFileMILPConditions( ...
    bucle, nP, inP_ini, ...
    nT, nR, nSlots, nMinsTot, ...
    tipus, ...
    nDies, dies_seleccio, ...
    nHosp, nUrg, ...
    T2p, ...
    Tdebug, ...
    lastID, maxID, nP_total, Total_Files, Total_Dies, ...
    lastMilp, default_inP)
%DEBUGPRINTFILEMILPCONDITIONS
% Funció de debug: imprimeix condicions d'entrada al MILP i registra el bloc a Tdebug.
% No retorna res (Tdebug es passa per referència? NO: a MATLAB no hi ha referència).
%
% IMPORTANT: si vols conservar el nou Tdebug, fes que la funció el retorni,
% o bé fes que Tdebug sigui global / handle. Ara mateix, aquesta crida
% genera un Tdebug nou però es perd fora d'aquesta funció.

    fprintf('------------ A PUNT PEL MILP -----------------------\n');
    fprintf('Al Bucle %d \n', bucle);
    fprintf('Per un total de procediments (np) %d \n', nP);
    fprintf('Iniciant en el procediment %d \n', inP_ini);
    fprintf('Finalitzant en el procediment %d \n', inP_ini+nP-1);
    fprintf('[nT %d, nR %d, nSlots %d, nMinsTot %d]\n', nT, nR, nSlots, nMinsTot);

    % --- Resum tipus procediments ---
    [nGastro, nColono, nOtros] = countProcedureTypes(tipus);

    fprintf('\n--- Resum de tipus de procediments ---\n');
    fprintf('Gastroscopies (30min) : %d\n', nGastro);
    fprintf('Colonoscopies (60min) : %d\n', nColono);
    fprintf('Otros (?min)          : %d\n', nOtros);

    % --- Dies seleccionats per programar ---
    printDiesSeleccio(dies_seleccio, nDies);

    fprintf('TOTAL Porcediments que es processaran(nP)            : %d\n', nP);
    fprintf('\nProcediments NO considerats en aquest model (perquè tenen el seu espai propi),\n');
    fprintf('només dins dels dies seleccionats:\n');
    fprintf('   Hospitalització (Hospitalizado): %d\n', nHosp);
    fprintf('   Urgents (Urgente):               %d\n', nUrg);
    fprintf('   TOTAL NO considerats:            %d\n', nHosp + nUrg);
    fprintf('\n');

    % --- Dies originals dels procediments seleccionats (T2p) ---
    dies_T2p = unique(T2p.FECHA);
    printDiesT2p(dies_T2p);

    % --- Registrar a Tdebug (side-effect: només útil si Tdebug es gestiona fora) ---
    %#ok<NASGU>  % (evita warning si no uses Tdebug després)
    Tdebug = afegeixDebugBloc( ...
        Tdebug, bucle, ...
        nP, inP_ini, ...
        nT, nR, nSlots, nMinsTot, ...
        tipus, ...
        nDies, dies_seleccio, ...
        nHosp, nUrg, ...
        dies_T2p, lastID, maxID, nP_total, Total_Files, Total_Dies, ...
        lastMilp, default_inP, "bucle");
end

% -------------------------------------------------------------------------
function [nGastro, nColono, nOtros] = countProcedureTypes(tipus)
    nGastro = sum(tipus == "Gastroscopia");
    nColono = sum(tipus == "Colonoscopia");
    nOtros  = sum(tipus == "Otros");
end

% -------------------------------------------------------------------------
function printDiesSeleccio(dies_seleccio, nDies)
    if numel(dies_seleccio) ~= nDies
        fprintf('\n ERROR CAS DIFERENT: Dies seleccionats PER PROGRAMAR (dies):%d i nDies %d\n', ...
            numel(dies_seleccio), nDies);
        return;
    end

    fprintf('\nDies seleccionats PER PROGRAMAR \n');
    if isempty(dies_seleccio)
        fprintf('Bloc sense dies seleccionats\n');
        return;
    end

    fmt = 'yyyy-MM-dd';
    fprintf('Dia Inici: %s, dia Fi: %s, total dies bloc %d, nDies %d\n', ...
        string(dies_seleccio(1), fmt), ...
        string(dies_seleccio(end), fmt), ...
        numel(dies_seleccio), nDies);
end

% -------------------------------------------------------------------------
function printDiesT2p(dies_T2p)
    fprintf('\n--- Dies originals dels procediments seleccionats (T2p) ---\n');
    for i = 1:numel(dies_T2p)
        % millor que datestr: imprimeix datetime com a string amb format
        fprintf('   - %s', string(dies_T2p(i), 'yyyy-MM-dd'));
    end
    fprintf('\nDies dels procediments in: %d\n', numel(dies_T2p));
end
