function y = fd(Ei,Ef,T)
y = zeros(size(Ei));
kb = 8.617333262E-5;
b = 1/(kb*T);
if(T==0)
   y = heaviside(Ef - Ei);
else
   for ie = 1 : length(Ei)
      y(ie) = (1+exp((Ei(ie)-Ef)*b))^(-1);
   end
end
