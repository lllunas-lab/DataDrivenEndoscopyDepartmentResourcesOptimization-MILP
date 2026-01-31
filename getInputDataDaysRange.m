function dies_T2p = getInputDataDaysRange (T2p)
    dies_T2p = unique(T2p.FECHA);   % dies presents a T2p
    fprintf('\n--- Dies originals dels procediments seleccionats (T2p) ---\n');
    for i = 1:numel(dies_T2p)
        fprintf('   - %s', datestr(dies_T2p(i), 'yyyy-mm-dd'));
    end
    fprintf('\nDies dels procediments in: %d\n', numel(dies_T2p));
end