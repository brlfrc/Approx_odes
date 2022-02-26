classdef Molecola
   properties
      x; y; z;
      p_x; p_y; p_z;
      delta;
   end

   methods
      function obj = Molecola(x,y,z,p_x, p_y, p_z, delta)
        if nargin == 0
           obj.x=0; obj.y=0; obj.z=0;
           obj.p_x=0; obj.p_y=0; obj.p_z=0;
           obj.delta=0;
        elseif nargin == 3
            obj.x=x; obj.y=y; obj.z=z;
            obj.p_x=0; obj.p_y=0; obj.p_z=0;
            obj.delta=1;
        elseif nargin == 6
            obj.x=x; obj.y=y; obj.z=z;
            obj.p_x=p_x; obj.p_y=p_y; obj.p_z=p_z;
            obj.delta=1;
        else
            obj.x=x; obj.y=y; obj.z=z;
            obj.p_x=p_x; obj.p_y=p_y; obj.p_z=p_z;
            obj.delta=delta;
        end
      end
      
      function r = distance(obj1,obj2)
         r = norm(obj1.return_Position(obj1)-obj2.return_Position(obj2));
      end

      function r = distance(obj1,obj2)
         r = norm(obj1.return_Position(obj1)-obj2.return_Position(obj2));
      end
   end

   methods (Static)
    function pos = return_Position(obj)
         pos = [obj.x; obj.y; obj.z];
      end
   
   
   end
end