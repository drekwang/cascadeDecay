function dpdt = odefun(t, p, kIn, kOut)
    dpdtIn = kIn*p;
    dpdtOut = zeros(length(p),1);
    for i = 1:length(p)
        dpdtOut(i, 1) = -1*sum(kOut(i,:)*p(i,1));
    end
    dpdt = dpdtIn + dpdtOut;
end