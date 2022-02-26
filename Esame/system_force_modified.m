function [f] = system_force_modified(z)
    
    if(size(z,2)~=1)
        z=z';
    end

    i=1;
    for j=1:length(z)/4
        q(:,j)=z(i:i+1,1);
        i=i+2;
    end
    for j=1:length(z)/4
        p(:,j)=z(i:i+1,1);
        i=i+2;
    end

    f=0;
    i=1;
    for a=1:length(q)
        f(i:i+1,1)=p(:,a);
        i=i+2;
    end
    for a=1:length(q)
        f(i:i+1,1)=forza_modified(q,a);
        i=i+2;
    end
end