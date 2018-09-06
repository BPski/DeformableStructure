classdef Link
    %methods and properties for a link
    %compute 
    properties
        init_length %the initial length
        bottom %bottom position
        top %the top position
        color %the color to draw with
    end
    properties (Access = private)
        stiffness
    end
%     properties (Dependent)
%         length %the current length
%         
%         %things that can be computed, but not set
%         strain
%         force
%         energy
%         normal_vector
%     end
    
    methods
        function obj = Link(init_length,stiffness,color)          
            
            if nargin == 2
                obj.color = [0;0;256]/256;
            else
                obj.color = color;
            end
%             if length(top) ~= 3 || length(bottom) ~= 3
%                 %make sure defined in 3D
%                 error("link expects positions to be 3D");
%             end
            
            obj.init_length = init_length;
            obj.stiffness = stiffness;
            obj.top = [0;0;init_length];
            obj.bottom = [0;0;0];
            
        end
        
        
        function length = length(obj)
            length = norm(obj.top-obj.bottom);
        end
        
        function obj = set.top(obj,val)
            if length(val) ~= 3
                error('top position of link needs to be in 3D')
            end
            obj.top = reshape(val,3,1);
        end
        
        function obj = set.bottom(obj,val)
            if length(val) ~= 3
                error('bottom position of link needs to be in 3D')
            end
            obj.bottom = reshape(val,3,1);
        end
        
        function strain = strain(obj)
            strain = (obj.length-obj.init_length)/obj.init_length;
        end
        
        function force = force(obj)
            force = -obj.stiffness*obj.strain*obj.normal_vector;
        end
        
        function energy = energy(obj)
            energy = (1/2)*obj.stiffness*obj.strain^2;
        end
        
        function vec = normal_vector(obj)
            vec = (obj.top-obj.bottom)/obj.length;
        end
        
        function plot(obj,offset)
            if nargin == 1
                plot3([obj.bottom(1),obj.top(1)],[obj.bottom(2),obj.top(2)],[obj.bottom(3),obj.top(3)],'color',obj.color);
            else
                top = offset*obj.top;
                bot = offset*obj.bottom;
                plot3([bot(1),top(1)],[bot(2),top(2)],[bot(3),top(3)],'color',obj.color)
            end
        end
    end
end