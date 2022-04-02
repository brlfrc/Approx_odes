function [f,count] = forza(q,a)  
    f=0;
    temp=q(:,a);
    q(:,a)=q(:,1);
    q(:,1)=temp;
    count=0;
    for j=2:length(q)
        count=count+1;
        f=f+4*(12/norm(q(:,1)-q(:,j))^13-6/norm(q(:,1)-q(:,j))^7)*(q(:,1)-q(:,j))/norm(q(:,1)-q(:,j));
    end
end
