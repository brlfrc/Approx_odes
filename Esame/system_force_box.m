function [f] = system_force_box (z,x_box,y_box)
    
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
        if abs(q(1,a))> x_box/2
            if (q(1,a)*p(1,a)>0)
                f(i,1)=-p(1,a);
            end
        end
        if abs(q(2,a))> y_box/2
            if (q(2,a)*p(2,a)>0)
                f(i+1,1)=-p(2,a);
            end
        end
        i=i+2;
    end
    for a=1:length(q)
        f(i:i+1,1)=forza(q,a);
        if abs(q(1,a))> x_box/2
            if(q(1,a)>0)
                f(i,1)=-500;
            else
                f(i,1)=500;
            end
        end
        if abs(q(2,a))> y_box/2
            if(q(2,a)>0)
                f(i+1,1)=-500;
            else
                f(i+1,1)=500;
            end
        end
        i=i+2;
    end
end