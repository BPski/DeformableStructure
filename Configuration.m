classdef Configuration
    %a representation for a homogenerous transformation
    %gives some convenient shorthand
    
    properties
        delta %the displacements to be considered
    end
    
    methods
        function obj = Configuration(val)
            %can be initialized either with the displacements or with a
            %matrix
           
            
            if nargin == 0 %initialize to identity
                obj.delta = [0;0;0;0;0;0];
            else
                [m,n] = size(val);
                if (m==6 && n==1) || (m==1 && n==6) %specified a delta
                    obj.delta = reshape(val,6,1);
                elseif m == 4 && n == 4 %specified a matrix
                    %check entries
                    if any([0,0,0,1]~=val(4,:)) %doesn't have the right shape
                        error("if specifying a matrix it is expected that the matrix looks like [R(3x3), p(3x1); 0(1x3), 1]")
                    end
                    %minus sign because the R and extractAngles seem to be
                    %opposite each other
                    obj.delta = [-extractAngles(val(1:3,1:3));val(1:3,4)];
                else
                    error("do not recognize how to make configuration from the entry, either need 6x1 or 4x4");
                end
            end
        end
        
        function R = R(obj)
            %extract rotation
            R = Rx(obj.delta(1))*Ry(obj.delta(2))*Rz(obj.delta(3));
        end
        
        function p = p(obj)
            %extract position
            p = obj.delta(4:6);
        end
        
        function M = M(obj)
            %the matrix representation
            M = [obj.R,obj.p;0,0,0,1];
        end
        
        function display(obj)
            display(obj.M);
        end
        
        function Ad = Ad(obj)
            %the adjoint transform
            Ad = [obj.R,zeros(3);skew(obj.p)*obj.R,obj.R];
        end
        
        function res = mtimes(obj,r)
            %overloads the matrix multiplication
            %the kinds of multiplication to cover are
                %3x1 vector, in this case return the 3x1 vector
                %anything else is regular matrix multiplcation
            if all(size(r) == [3,1])
                v = obj.M*[r;1]; %could have issues if the vector is not "free"
                res = v(1:3);
            elseif isa(r,'Configuration') %retain configuration between configurations
                res = Configuration(obj.M*r.M);
            else
                res = obj.M*r;
            end
        end
        
        function obj = set.delta(obj,val)
            if all(size(val)==[6,1]) || all(size(val)==[1,6])
                obj.delta = reshape(val,6,1);
            else
                error("delta needs to be set as a 6x1 vector [3x angles, 3x displacements]");
            end
        end
    end
end
    
        