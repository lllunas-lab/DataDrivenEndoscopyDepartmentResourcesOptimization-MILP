function sol = executarMILP (nP, nR, nT, avail, Delta, d)
    original = false;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Variables:
    % - Variables:
    %   x_{i,k,t} (binàries): el procediment i comença a sala k, slot temporal t
    %   u_k       (contínua): minuts d'ús total de la sala k
    %   idle_k    (contínua): minuts d'oci total de la sala k
    %   delta_k   (contínua): desviació d'ús de la sala k respecte a mu
    %
    % - Indexació natural MATLAB per x:
    %   idx_(i,k,t) = i + (k-1)*nP + (t-1)*nP*nR
    %
    
    % - Objectiu: prioritzar puntualitat + minimitzar oci + balancejar ús entre sales

    %% sensibilitat
    % (1) concentra massa en una sala → puja w_bal 
    % (2) t'omple massa tard el dia → puja w_start
    % (3) et deixa molt oci perquè no hi ha prou demanda → baixa w_idle o fes mu més petit

    
    %% 3. Indexació de les variables x_{i,k,t}
    % Vector de variables z: , de mida nP*nR*nT
    % Primer totes les x_{i,k,t}, després u_k, delta_k i idle_k
    fprintf('------------ INICI CÀLCUL VARIABLES PROCÉS -----------------------\n');
    nX = nP * nR * nT;       % nombre de variables x
    fprintf('\n\n nX(nombre variables) = %d\n nP(procediments) %d* \n nR(sales) %d* \n nT(slots) %d \n', nX, nP, nR, nT);
    theta = ceil(d / Delta);   % nombre de slots per procediment

    % 0j0 nVars ara és nX 0j0 0j0 0j0
    offsetU     = nX;               % u_k comença a nX+1
    offsetDelta = nX + nR;          % delta_k comença a nX+nR+1
    offsetIdle  = nX + 2*nR;        % idle_k comença a nX+2*nR+1
    
    nVars = nX + 3*nR;       % Nombre total
    % nVars = nP * nR * nT; ESBORRAR
    
    
    %% Indexació variables
    % Funció idx_ coherent amb la indexació lineal de MATLAB:
    % idx(i,k,t) er retorna la posició a la variable vectoritzada 
    idx_ = @(i,k,t) i + (k-1)*nP + (t-1)*nP*nR;
    upos    = @(k)      offsetU     + k;              % per u_k
    deltpos = @(k)      offsetDelta + k;              % per delta_k
    idlepos = @(k)      offsetIdle  + k;              % per idle_k
    
    
    %% 4. Tipus de variable, límits i funció objectiu (punctualitat)
    intcon = 1:nX;              % només les x són enteres
    lb = zeros(nVars,1);        % límit inferior
    ub = inf(nVars,1);          % límit superior
    ub(1:nX,1)=1;               % límit superior (binàries); x ∈ [0,1]; la resta [0, inf)
    fprintf("%% 4. Tipus de variable, límits i funció objectiu (punctualitat)");
    if original
        % Bloquegem els slots no disponibles: x_{i,k,t} no pot ser 1 si avail(k,t)=0
    
        for k = 1:nR
            for t = 1:nT
                if avail(k,t) == 0
                    for i = 1:nP
                        ub(idx_(i,k,t)) = 0;  % força x(i,k,t)=0
                    end
                end
            end
        end

    else %% [NOU] 4b. Límits superiors per x: T_i + no-spill amb avail
        for i = 1:nP
            tMaxStart = nT - theta(i) + 1;   % últim slot de començament vàlid (T_i)

            if tMaxStart < 1
                error('El procediment %d (theta=%d) no cap en l''horitzó nT=%d.', i, theta(i), nT);
            end

            % 1) Fora de T_i: t > tMaxStart  -> ub = 0 per a totes les sales
            if tMaxStart < nT
                tForbidden = (tMaxStart+1):nT;  % slots de començament que "surten" de l'horitzó
                for k = 1:nR
                    ub(idx_(i,k,tForbidden)) = 0;
                end
            end

            % 2) No-spill respecte a avail(k,·):
            %    si algun slot [t, t+theta(i)-1] està bloquejat, el començament t no és vàlid
            for k = 1:nR
                for t = 1:tMaxStart
                    span = t:(t+theta(i)-1);              % slots que ocuparia el procediment
                    if any(avail(k, span) == 0)
                        ub(idx_(i,k,t)) = 0;              % força x_{i,k,t} = 0
                    end
                end
            end
        end
    end
    % FI [NOU] 4b. Límits superiors per x: T_i + no-spill amb avail
    fprintf('% NOU: Càlcul capacitat per sala (minuts disponibles) i u-objectiu mu');
    % NOU: Càlcul capacitat per sala (minuts disponibles) i u-objectiu mu
    cap_k = zeros(nR,1);
    for k = 1:nR
        cap_k(k) = sum(avail(k,:)) * Delta;   % slots lliures * 15'
    end
    
    total_demand = sum(d);
    mu = min( total_demand/nR, mean(cap_k) ); % objectiu d'ús per sala
    
    
    %____ES POT EXPERIMENTAR_________________________________________________
    % Pesos objectiu
    w_start = 1.0;    % prioritzar començar aviat
    w_idle  = 0.01;   % penalitzar oci
    w_bal   = 0.1;    % penalitzar desequilibri entre sales
    %________________________________________________________________________
    
    
    %% Vector f
    f = zeros(nVars,1);
    % Funció objectiu f'*z
    
    %________________________________________________________________________
    % Puntualitat: minimitzar la suma de temps de començament.
    % Terme de puntualitat: w_start * sum (t-1)*Delta * x_{i,k,t}
    
    for t = 1:nT
        cost_t = (t-1)*Delta;   % minuts des de l'inici del torn
        for k = 1:nR
            for i = 1:nP
                f(idx_(i,k,t)) = f( idx_(i,k,t) ) + w_start*cost_t; % Igual que l'anterior si w_start=1
            end
        end
    end
    
    % plot(f,'.') % Control plot
    %________________________________________________________________________
    
    %________________________________________________________________________
    % Terme d'oci total: w_idle * sum idle_k
    for k = 1:nR
        f( idlepos(k) ) = f( idlepos(k) ) + w_idle;
    end
    %________________________________________________________________________
    
    %________________________________________________________________________
    % Terme de balanceig: w_bal * sum delta_k
    for k = 1:nR
        f( deltpos(k) ) = f( deltpos(k) ) + w_bal;
    end
    %________________________________________________________________________
    
    
    
    
    
    
    
    %% 5. Restriccions d'igualtat: cada procediment una única vegada (Aeq, beq)
    % 5.1. Cada procediment exactament un cop:
    %      sum_{k,t} x_{i,k,t} = 1
    Aeq = sparse(nP, nVars);
    beq = ones(nP,1);
    
    %% [NOU] 5.1. Assignació única: sum_{k,t in T_i} x_{i,k,t} = 1

    if original
        % VELL Vectorització parcial: per cada i, totes les seves variables
        % estan en posicions i:nP:nVars
        for i = 1:nP
            Aeq(i, i:nP:nX) = 1;
        end
    else %% [NOU] 5.1. Assignació única: sum_{k,t in T_i} x_{i,k,t} = 1
        for i = 1:nP
            tMaxStart = nT - theta(i) + 1;      % límit superior de T_i
            tRange_i  = 1:tMaxStart;            % T_i

            for k = 1:nR
                cols = idx_(i, k, tRange_i);    % totes les x_{i,k,t} amb t ∈ T_i
                Aeq(i, cols) = 1;
            end
        end
    end
    % FI [NOU] 5.1. Assignació única:
    
    
    % 5.2. Definició d'u_k:
    %      u_k = sum_i d_i * sum_t x_{i,k,t}
    % ->   u_k - sum_i d_i * sum_t x_{i,k,t} = 0
    
    Aeq_u  = sparse(nR, nVars);
    beq_u  = zeros(nR,1);
    
    %% [NOU] 5.2. Definició d'u_k amb t ∈ T_i
    % VELL
    if original
        for k = 1:nR
            Aeq_u(k, upos(k)) = 1;
            for i = 1:nP
                for t = 1:nT
                    Aeq_u(k, idx_(i,k,t)) = Aeq_u(k, idx_(i,k,t)) - d(i);
                end
            end
        end
    else
        for k = 1:nR
            Aeq_u(k, upos(k)) = 1;      % coeficient de u_k

            for i = 1:nP
                tMaxStart = nT - theta(i) + 1;   % T_i
                tRange_i  = 1:tMaxStart;

                cols = idx_(i, k, tRange_i);
                Aeq_u(k, cols) = Aeq_u(k, cols) - d(i);
            end
        end
    end
    %FI [NOU] 5.2. Definició d'u_k amb t ∈ T_i
    
    
    % 5.3. Relació ús + oci = capacitat: u_k + idle_k = cap_k
    Aeq_cap = sparse(nR, nVars);
    beq_cap = cap_k;
    
    for k = 1:nR
        Aeq_cap(k, upos(k))   = 1;
        Aeq_cap(k, idlepos(k)) = 1;
    end
    
    % Combinar totes les Condicions que s'han de complir amb igualtat 
    Aeq = [Aeq; Aeq_u; Aeq_cap];
    beq = [beq; beq_u; beq_cap];
    
    %% 6. Restriccions d'inequació: no solapament a cada sala i temps
    % Per a cada sala k i franja temporal t:
    %   sum_i sum_{tau in [t-theta(i)+1, t]} x_{i,k,tau} <= 1
    %
    % Versió semi-vectoritzada: bucles sobre k i t, i assignació vectorial sobre tau.
    
    nConstrOverlap = nR * nT;
    A_overlap = spalloc(nConstrOverlap, nVars, nConstrOverlap * nP * 2); % estimació
    b_overlap = ones(nConstrOverlap, 1);
    
    row = 0;
    %% [NOU] 6.1. No solapament: tau ∈ T_i ∩ [max(1,t-theta(i)+1), t]

    % VELL
    if original
        for k = 1:nR
            for t = 1:nT
                row = row + 1;
                % Si la sala k en el temps t està bloquejada (avail=0),
                % podem ometre parcialment (igualment ja hem posat ub=0)
                for i = 1:nP
                    span = theta(i);
                    tau_min = max(1, t - span + 1);
                    tau_max = t;
                    if tau_min <= tau_max
                        taus = tau_min:tau_max;
                        cols = idx_(i, k, taus);
                        A_overlap(row, cols) = A_overlap(row, cols) + 1;
                    end
                end
            end
        end
    else
        % [NOU] 6.1 No solapament amb T_i
        nConstrOverlap = nR * nT;
        A_overlap = spalloc(nConstrOverlap, nVars, nConstrOverlap * nP * 2);
        b_overlap = ones(nConstrOverlap, 1);

        row = 0;
        for k = 1:nR
            for t = 1:nT
                row = row + 1;
                for i = 1:nP
                    span    = theta(i);
                    tau_min = max(1, t - span + 1);
                    tau_max = min(t, nT - span + 1);   % tall amb T_i

                    if tau_min <= tau_max
                        taus = tau_min:tau_max;
                        cols = idx_(i, k, taus);
                        A_overlap(row, cols) = A_overlap(row, cols) + 1;
                    end
                end
            end
        end
    end
    
    % 6.2. Balanceig: |u_k - mu| <= delta_k
    %      u_k - mu <= delta_k  -> u_k - delta_k  <= mu
    %      mu - u_k <= delta_k  -> -u_k - delta_k <= -mu
    A_bal = sparse(2*nR, nVars);
    b_bal = zeros(2*nR,1);
    
    
    
    for k = 1:nR
        % u_k - delta_k <= mu
        A_bal(2*k-1, upos(k))    =  1;
        A_bal(2*k-1, deltpos(k)) = -1;
        b_bal(2*k-1)             =  mu;
    
        % -u_k - delta_k <= -mu
        A_bal(2*k, upos(k))      = -1;
        A_bal(2*k, deltpos(k))   = -1;
        b_bal(2*k)               = -mu;
    end
    
    % Combinar totes les inequacions
    A = [A_overlap; A_bal];
    b = [b_overlap; b_bal];
    
    
    %% 7. Crida a intlinprog
    %Options del Pere: 
    % options = optimoptions('intlinprog', ...
    %    'Display','iter', ...
    %    'Heuristics','advanced', ...
    %    'CutGeneration','basic');
    
    % Options Legacy
    % FIX: evita el backend HiGHS (bug intern amb "optimstatus" en alguns casos)
    %options = optimoptions('intlinprog', ...
    %'Algorithm','legacy', ...
    %'Display','iter', ...
    %'Heuristics','advanced', ...
    %'CutGeneration','basic');


    % Options per entendre el trigger de l'error
    options = optimoptions('intlinprog', ...
       'Display','iter');   % sense Heuristics ni CutGeneration

    disp('__________INI____________')
    tMILP = tic;
    [z_opt, fval, exitflag, output] = intlinprog( ...
        f, intcon, A, b, Aeq, beq, lb, ub, options);
    elapsed=toc(tMILP);
    disp('____________FINI____________')
    fprintf('Elapsed time MILP: %.6f seconds\n', elapsed);
    elapsed_hms = seconds(elapsed);
    fprintf('Temps transcorregut pel MILP: %s\n', string(elapsed_hms, 'hh:mm:ss'));

    disp(output)
    fieldnames(output)

    ver('optim')
    options = optimoptions('intlinprog'); 
    disp(options.Algorithm)
    

    if exitflag ~= 1
        warning('MILP no òptim (exitflag = %d).', exitflag);
        fprintf('Atenció: intlinprog no ha retornat exitflag = 1 (òptim trobat).\n');
        % Empaquetar mínim retorn
        sol = struct();
        sol.exitflag = exitflag;
        sol.output   = output;
        sol.fval     = fval;

        return
    end

    
    % Reconstruïm el tensor x(i,k,t) directament amb reshape
    % x només són les primeres nX variables
    x_vec = z_opt(1:nX);
    x_opt = reshape(x_vec, [nP nR nT]);
    
    u_opt    = z_opt(offsetU+1     : offsetU+nR);
    delta_opt= z_opt(offsetDelta+1 : offsetDelta+nR);
    idle_opt = z_opt(offsetIdle+1  : offsetIdle+nR);
    
    
    %% 9. Mostrar l'horari resultant
    nAssignat = 0;

    if nT ~= size(x_opt,3)
        warning('Inconsistència temporal: nT=%d però size(x_opt,3)=%d. Ajustant impressió a nT_real.', ...
            nT, size(x_opt,3));
        nT_real = size(x_opt, 3);
        slotRangeDia = @(dia) ((dia-1)*slotsDia + 1) : min(dia*slotsDia, nT_real);
        nDies = ceil(nT_real / slotsDia);
        fprintf('\nDies real (%d dies):\n i nDies (dies_selecció) %d, nDies %d \n',nDies_real, numel(dies_seleccio), nDies);   
    end

    %% --- Empaquetar resultats en "sol" (retorn de la funció) ----------------
    sol = struct();

    % Estat del solver
    sol.exitflag = exitflag;     % 1 = èxit, altres = no òptim / infeasible / etc.
    sol.output   = output;       % info detallada del solver
    sol.fval     = fval;         % valor objectiu

    % Solució (decisions)
    sol.x_opt     = x_opt;       % [nR x nP x nT] assignació sala-proc-slot
    sol.u_opt     = u_opt;       % [nP x 1] start time en minuts o en unitats (segons ho defineixis)
    sol.idle_opt  = idle_opt;    % [nR x 1] idle per sala (si ho tens així)
    sol.delta_opt = delta_opt;   % [nR x 1] desviació/càrrega (segons model)
    sol.mu        = mu;          % [nR x 1] càrrega o ús (segons model)

    % Metadades mínimes útils (no és "in", és per traçabilitat)
    sol.nP    = nP;
    sol.nR    = nR;
    sol.nT    = nT;
    sol.Delta = Delta;
    sol.d     = d;
    sol.avail = avail;   % context de disponibilitat utilitzat pel MILP
    sol.elapsed = elapsed;


    % (Opcional) checks/derivats útils
    % Percentatge assignat (si x_opt és 0/1):
    try
        sol.nAssignats = sum(any(any(x_opt,1),3));  % compte procediments amb alguna assignació
    catch
        % si la dimensió no encaixa, ignorem
    end

end