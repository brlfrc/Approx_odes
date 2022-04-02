function [U] = U_potenziale_box(q,x_box,y_box,U0)
%   U_potenziale Data configurazione calcola il potenziale totale di lennard jones
%   (caso 2d)
%   

sigma=1;
epsilon=1;
U=0;
for i=1:length(q)
    for j=i+1:length(q)
       U=U+4*epsilon*(sigma/norm(q(:,i)-q(:,j))^12-(sigma/norm(q(:,i)-q(:,j))^6));
    end
end

for a=1:length(q)
    if abs(q(1,a))> x_box/2
       U=U+U0*(q(1,a)-x_box/2)*(x_box/2+q(1,a));
    end
    if abs(q(2,a))> y_box/2
       U=U+U0*(q(2,a)-y_box/2)*(y_box/2+q(2,a));
    end
end

end