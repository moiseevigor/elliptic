function m = nome2m(q)
%NOME2M  Inverse of Moiseev's nomeq:  q  ->  m  (0<q<1).

    assert(all(abs(q)<1 & q>0),'q must satisfy 0<q<1.');
    f = @(m) nomeq(m) - q;          % nomeq is forward map m -> q
    m = arrayfun(@(qq) fzero(f, [1e-8 1-1e-8]), q);   % vectorised
end

