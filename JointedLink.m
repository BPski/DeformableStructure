classdef JointedLink < Link
    %the pseudo link that forms when the truss, joint, truss structure is
    %approximated as a link by saying that the joint is much softer
    
    properties
        B
        beta0
        joint_stiffness
    end
    
    properties (Dependent)
        stiffness
    end
    
    methods
        function obj = JointedLink(B,beta0,joint_stiffness)          
            %assuming b_up=b_down=B
%             if length(top) ~= 3 || length(bottom) ~= 3
%                 %make sure defined in 3D
%                 error('link expects positions to be 3D');
%             end
            
            init_length = sqrt(2*B^2*(1-cos(beta0)));
            obj@Link(init_length,joint_stiffness);
            obj.B = B;
            obj.beta0 = beta0;
            obj.joint_stiffness = joint_stiffness;
            obj.top = [0;0;obj.init_length];
            obj.bottom = [0;0;0];
            
        end
        
        
%         function length = get.length(obj)
%             length = norm(obj.top-obj.bottom);
%         end
%         
%         function obj = set.top(obj,val)
%             if length(val) ~= 3
%                 error('top position of link needs to be in 3D')
%             end
%             obj.top = reshape(val,3,1);
%         end
%         
%         function obj = set.bottom(obj,val)
%             if length(val) ~= 3
%                 error('bottom position of link needs to be in 3D')
%             end
%             obj.bottom = reshape(val,3,1);
%         end
%         
%         function strain = get.strain(obj)
%             strain = (obj.length-obj.init_length)/obj.init_length;
%         end
        
        function stiffness = get.stiffness(obj)
            %the actual stiffness
            b = obj.length();
            b0 = obj.init_length;
            B = obj.B;
            beta = obj.beta0;
            KJ = obj.joint_stiffness;
            del = b-b0;
            if b ~= obj.init_length
                %stiffness = obj.joint_stiffness*obj.init_length^2*((acos((obj.b_up0^2+obj.b_down0^2-b^2)/(2*obj.b_up0*obj.b_down0))-obj.beta0)/(b-obj.init_length))^2;
                %stiffness = obj.joint_stiffness*obj.init_length*(acos(cos(obj.beta0)-(del^2+2*del*obj.init_length)/(2*obj.b_up0*obj.b_down0))-obj.beta0)/del;
                stiffness = KJ*b0/B*(acos(cos(beta)-(del^2+2*del*b0)/(2*B^2))-beta)/(del*sqrt(1-b0^2/(4*B^2)-(del^2+2*del*b0)/(4*B^2)));
            else
                %stiffness = obj.joint_stiffness*obj.init_length^4/(obj.b_up0*obj.b_down0*sin(obj.beta0))^2;
                %stiffness = obj.joint_stiffness*obj.init_length^2/(obj.b_up0*obj.b_down0*sin(obj.beta0));
                stiffness = KJ*4*b0^2/(4*B^2-b0^2);
            end
        end
        
        function energy = energy(obj)
%             %have to approximate the energy due to the stiffness being
%             %nonlinear
%             %using an integration of a 3rd order taylor expansion
%             del = obj.length-obj.init_length;
%             b0 = obj.init_length;
%             b_up = obj.b_up0;
%             b_down = obj.b_down0;
%             KJ = obj.joint_stiffness;
%             beta = obj.beta0;
%             energy = -KJ*(del^2*(b0/(b_up*b_down*sin(beta)))/2+del^3*(1/(2*b_up*b_down*sin(beta))-(b0^2*cos(beta))/(2*(b_up*b_down)^2*sin(beta)^3))/3);
            %have to approximate the integral somehow, use trapezoidal rule
            %to approximate the integral
            energy = 0;
            %start at undeformed length
            link_copy = JointedLink(obj.B,obj.beta0,obj.joint_stiffness); %make a copy to modify
            b0 = obj.init_length;
            b = obj.length;
            n = 100; %number of discretizations
            del = b-b0;
            dl = del/n;
            
            %the first step (not multiplied by 2
            i = 0;
            b = i*dl;
            link_copy.top = [0;0;b+b0];
            energy = energy + link_copy.stiffness*(b-b0);
            for i=1:n-1
                b = i*dl;
                link_copy.top = [0;0;b+b0];
                energy = energy + 2*link_copy.stiffness*(b-b0);
            end
            %last step
            i = n;
            b = i*dl;
            link_copy.top = [0;0;b+b0];
            energy = energy + link_copy.stiffness*(b-b0);
            
            energy = dl/2*energy; %multiply by a factor
        end
            
        
%         function force = get.force(obj)
%             force = -obj.stiffness*obj.strain*obj.normal_vector;
%         end
%         
%         function energy = get.energy(obj)
%             energy = (1/2)*obj.stiffness*obj.strain^2;
%         end
%         
%         function vec = get.normal_vector(obj)
%             vec = (obj.top-obj.bottom)/obj.length;
%         end
%         
        function plot(obj,offset)
            if nargin == 1
                plot3([obj.bottom(1),obj.top(1)],[obj.bottom(2),obj.top(2)],[obj.bottom(3),obj.top(3)],'y');
            else
                top = offset*obj.top;
                bot = offset*obj.bottom;
                plot3([bot(1),top(1)],[bot(2),top(2)],[bot(3),top(3)],'g')
            end
        end
    end
end