function S1 = S1solver(S2,muk, alpha, gamma, k1, k3)
S10 = 0;
S1 = fzero(@(S1)S1eq(S1, S2,muk, alpha, gamma, k1, k3),S10);
end

function eq = S1eq(S1, S2,muk, alpha, gamma, k1, k3)
eq = S1 + muk.*norm(S2).*(alpha + k3*S1)./(norm(alpha + k3*S1).*(1 + (gamma./ (alpha + (k3 - k1).*S1) )^2)^(1/2));
end