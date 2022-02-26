function [U] = U_potenziale(q)
%   U_potenziale Data configurazione calcola il potenziale totale di lennard jones
%   (caso 2d)
%   

sigma=1;
epsilon=1;
U=0;
for i=1:length(q)
    for j=i+1:length(q)
       U=U+4*epsilon*(sigma/norm(q(:,1)-q(:,j))^12-(sigma/norm(q(:,1)-q(:,j))^6));
    end
end

end