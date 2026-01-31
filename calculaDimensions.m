function [nT, nR, nSlots, nMinsTot] = calculaDimensions(nDies, slotsDia, Delta, sales)
% calculaDimensions
%   Calcula les dimensions globals del problema de planificació
%
%   Inputs:
%     nDies     : nombre de dies a programar
%     slotsDia  : nombre de slots per dia
%     Delta     : durada d'un slot (minuts)
%
%   Outputs:
%     nT        : slots globals (total slots per sala en els nDies previstos de programar)
%     nR        : nombre de sales (ROOM 1..4)
%     nSlots    : Total slots en les nR Sales
%     nMinsTot  : Total minuts en les nR Sales

    nT         = nDies * slotsDia;      % slots globals (total slots per sala en els nDies previstos de programar)
    nR         = sales;                     % nombre de sales (ROOM 1..4)
    nSlots     = nT * nR;               % Total slots en les nR Sales
    nMinsTot   = nSlots * Delta;        % Total minuts en les nR Sales
end
