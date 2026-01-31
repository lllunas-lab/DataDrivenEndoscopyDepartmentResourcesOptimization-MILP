%% ----------------------- FUNCIONS AUXILIARS -------------------%%

% Exemple 1 de reserves (esquemàtic, per veure el mecanisme):
%  - Sala 3: primeres 4 franges (1h) reservades a urgents -> bloquejades
%  - Sala 1: franges 9-12 (2h, de 2:00 a 4:00) reservades a hospitalitzats
%  - Sala 2: franges 9-10 (30') també reservades
%
% URG a la Sala 3 (slots 1..4)
% avail(3, 1:4) = 0;
%
% HOSP a Sala 1 (slots 9..12)
% avail(1, 9:12) = 0;
%
% HOSP parcial a Sala 2 (slots 9..10)
% avail(2, 9:10) = 0;


% Exemple de reserves (esquemàtic, per veure el mecanisme):

% URGÈNCIES
% Dilluns calen 3 endoscopies de reserva
% Suposem temps màxim: 60mins o 4 slots (60/15)
% Per tant dilluns calen 12 Slots reservats

% Resta de dies calen 2 endoscopies de reserva
% Suposem temps màxim: 60mins o 4 slots (60/15)
% Per tant dilluns calen 8 Slots reservats

% Sala 1: PER ARA FEM TOTES LES RESERVES EN AQUESTA SALA


%Sala 2, 3 i 4 iguals
%avail(2:4, 1:8) = 0;            % Dia 1 les 2 primeres hores
%avail(2:4, slotsDia+(1:4)) = 0; % Dia 2 les 2 primeres hores

function [avail, reserveType] = buildAvailReserves(nDies, nR, nT, slotsDia, dies_programacio)

    % Inicialització:
    % 1 = disponible, 0 = reservat
    avail = ones(nR, nT);
    reserveType = strings(nR, nT);

    % --- Robustesa bàsica ---
    if isempty(dies_programacio)
        warning('dies_programacio és buit: no s''aplica cap reserva.');
        return;
    end

    % Assegurem consistència: nDies ve donat per dies_programacio
    nDiesReal = numel(dies_programacio);
    if nDiesReal ~= nDies
        warning('nDies (%d) difereix de numel(dies_programacio) (%d). S''usarà dies_programacio.', nDies, nDiesReal);
        nDies = nDiesReal;
    end

    % --- Aplicació de reserves per cada dia segons weekday ---
    for d = 1:nDies
        base = (d-1) * slotsDia;   % offset dins de l'eix temps (1..nT)

        wd = weekday(dies_programacio(d)); % 1=dg,2=dl,...,7=ds

        if wd == 2
            % =========================
            % DILLUNS
            % =========================
            % Sala 1: 12 slots (3h)
            idx1 = base + (1:12);
            idx1 = idx1(idx1 >= 1 & idx1 <= nT);
            avail(1, idx1) = 0;
            reserveType(1, idx1) = "RESERVA";

            % Sala 2: 4 slots (1h)
            idx2 = base + (1:4);
            idx2 = idx2(idx2 >= 1 & idx2 <= nT);
            if nR >= 2
                avail(2, idx2) = 0;
                reserveType(2, idx2) = "RESERVA";
            end

        elseif wd >= 3 && wd <= 6
            % =========================
            % DIMARTS..DIVENDRES (patró "dia 2" i successius)
            % =========================
            % Sala 1: 8 slots (2h)
            idx1 = base + (1:8);
            idx1 = idx1(idx1 >= 1 & idx1 <= nT);
            avail(1, idx1) = 0;
            reserveType(1, idx1) = "RESERVA";

            % Sala 2: 4 slots (1h)
            idx2 = base + (1:4);
            idx2 = idx2(idx2 >= 1 & idx2 <= nT);
            if nR >= 2
                avail(2, idx2) = 0;
                reserveType(2, idx2) = "RESERVA";
            end

        else
            % Caps de setmana o valors inesperats
            % Si preferiu, podeu fer warning aquí:
            % warning('Data %s no és laborable (weekday=%d): no s''apliquen reserves.', string(dies_programacio(d)), wd);
        end
    end
end
