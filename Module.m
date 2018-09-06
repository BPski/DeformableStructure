classdef Module
    %define the modules of a structure
    %the two plates and the links between them
    %for this they are allowed to be general asymmetric modules

    %for now the assumption for drawing is that just draw the polygon formed
    %by the points and draw the line between the centers
    properties
        r_up % the radii for the points on the upper plate
        r_down %radii for lower plate

        theta_up %angles for upper plate
        theta_down %angles lower plate

        g %the configuration matrix between the two plates

        n %number of polygons sides (mainly for reference)

        links %a list of links
    end
    
    properties (SetAccess = private)
        init_delta %the initial delta for resetting purposes
    end

    properties (Dependent)
        delta %from g
    end

    methods
        function obj=Module(r_up,r_down,theta_up,theta_down,links)
            %initialize the object, main objective is to determine delta
            %when undeformed

            %since this is the general case delta is not gauranteed to be
            %solved unless it is a trangular shape (or only 3 points given)

            n = length(r_up);
            if any([length(r_down),length(theta_up),length(theta_down)] ~= n)
                error('not all properties specified are consistent');
            end

            %setting the known properties and making sure they are in the
            %desired shape (column vectors)
            obj.n = n;         
            obj.r_up = reshape(r_up,n,1);
            obj.r_down = reshape(r_down,n,1);
            obj.theta_up = reshape(theta_up,n,1);
            obj.theta_down = reshape(theta_down,n,1);

            obj.g = Configuration();

            %setup the links
            obj.links = links; %this is a cell array, very important
%             for i=1:n
%                 %go in the order of b then a through the iterations
%                 %the positions aren't actually known yet due to delta not
%                 %being determined, know the bottom
%                 b = link(bs(i),Kb(i),Rz(theta_down(i))*[r_down(i);0;0],Rz(theta_up(i))*[r_up(i);0;0]+[0;0;1]);
%                 if i ~= n
%                     a = link(as(i),Ka(i),Rz(theta_down(i+1))*[r_down(i+1);0;0],Rz(theta_up(i))*[r_up(i);0;0]+[0;0;1]);
%                 else
%                     a = link(as(i),Ka(i),Rz(theta_down(1))*[r_down(1);0;0],Rz(theta_up(i))*[r_up(i);0;0]+[0;0;1]);
%                 end
%                 obj.links = [obj.links,b,a];
%             end
            %do the links
            for i=1:2:length(obj.links)
                j = round(i/2);
                %b links
                obj.links{i}.bottom = Rz(theta_down(j))*[r_down(j);0;0];
                obj.links{i}.top = Rz(theta_up(j))*[r_up(j);0;0]+[0;0;obj.links{i}.init_length];
                %a links
                if i+1 ~= 2*n
                    obj.links{i+1}.bottom = Rz(theta_down(j+1))*[r_down(j+1);0;0];
                    obj.links{i+1}.top = Rz(theta_up(j))*[r_up(j);0;0]+[0;0;obj.links{i+1}.init_length];
                else
                    obj.links{i+1}.bottom = Rz(theta_down(1))*[r_down(1);0;0];
                    obj.links{i+1}.top = Rz(theta_up(n))*[r_up(n);0;0]+[0;0;obj.links{i+1}.init_length];
                end
            end

            %setup the initial configuration
            obj.g = [0;0;pi/n;0;0;obj.links{1}.init_length];
            obj = obj.minimize('strain');
            obj.init_delta = obj.delta;
        end
        
        function obj = reset(obj)
            obj.g = obj.init_delta;
        end

        function plot(obj,offset)
            %draw the module

            if nargin == 1 %no offset specified
                g_off = Configuration();
            elseif isa(offset,'Configuration') %gave a configuration
                g_off = offset;
            elseif length(offset) == 6 %gave a delta
                g_off = Configuration(offset);
            else
                error('not sure what to do with the specified offset');
            end
            
            top = [];
            bottom = [];

            %extract info from the b links
            for i=1:2:2*obj.n
                top = [top,g_off*obj.links{i}.top];
                bottom = [bottom,g_off*obj.links{i}.bottom];
            end
            %add in a cycle
            top = [top,top(:,1)];
            bottom = [bottom,bottom(:,1)];
            hold on
            plot3(top(1,:),top(2,:),top(3,:),'r')
            plot3(bottom(1,:),bottom(2,:),bottom(3,:),'r')

            %the center line
            center_bottom = g_off.delta(4:6);
            center_top = g_off*obj.delta(4:6);
            plot3([center_bottom(1);center_top(1)],[center_bottom(2);center_top(2)],[center_bottom(3);center_top(3)],'black');

            %links
            for i=1:length(obj.links)
                plot(obj.links{i},g_off)
            end
            hold off
        end

        function obj = updateLinks(obj)
            %it just so happens that the links need to be updated when a
            %delta is set, but I get some warnings, see if this avoids it

            %update all the top positions
            tops = zeros(3,obj.n);
            for i=1:obj.n
                tops(:,i) = obj.g*(Rz(obj.theta_up(i))*[obj.r_up(i);0;0]);
            end

            for i=1:obj.n*2
                obj.links{i}.top = tops(:,round(i/2));
            end
        end


        function obj = set.g(obj,val)
            %set the delta in g
            if isa(val,'Configuration') %gave a config 
                obj.g = val;
            elseif length(val) ~= 6 %gave a delta
                error('delta requires 6 items to specify g [3x angles, 3x displacements]')
            else
                obj.g.delta = val;
            end

            obj = obj.updateLinks();
        end

        function delta = get.delta(obj)
            delta = obj.g.delta;
        end

        function out = total(obj,name,delta)
            %generalize the total___ computations, though it is dumb to be
            %checking against a string I think

            if nargin == 3
                obj.g = delta;
            end
            out = 0;
            if strcmpi('strain',name) %total strain
                for i=1:length(obj.links)
                    out = out + abs(obj.links{i}.strain);
                end
            elseif strcmpi('force',name) %total force
                for i=1:length(obj.links)
                    out = out + obj.links{i}.force;
                end
            elseif strcmpi('moment',name) %total moment
                for i=1:length(obj.links)
                    out = out + cross(obj.links{i}.bottom,obj.links{i}.force);
                end
            elseif strcmpi('wrench',name) %total wrench
                for i=1:length(obj.links)
                    out = out + [cross(obj.links{i}.bottom,obj.links{i}.force);obj.links{i}.force];
                end
            elseif strcmpi('energy',name) %total energy
                for i=1:length(obj.links)
                    out = out + obj.links{i}.energy;
                end
            else
                error('the provided name is not recognized');
            end
        end


        function obj = balanceWrench(obj,wrench)
            %given an applied wrench (center of top plate) determine delta
            %to balance it
            obj.g = fsolve(@(delta) obj.g.Ad'*wrench-obj.total('wrench',delta),obj.delta);
        end
        
        function obj = balanceForce(obj,force)
            %given an applied wrench (center of top plate) determine delta
            %to balance it
            obj.g = fsolve(@(delta) obj.g.Ad*[0;0;0;force]-[0;0;0;obj.total('force',delta)],obj.delta);
        end

        function obj = minimize(obj,name)
            %generalize computations for the minimization via changing delta
            obj.g = fmincon(@(delta) obj.total(name,delta),obj.delta);
        end

        function obj = minimizePotential(obj,wrench)
            %determine the change to delta that results in the minimum
            %potential = U-W*(delta-delta0)
            obj.g = fmincon(@(delta) abs(obj.total('energy',delta)-dot(wrench,obj.init_delta-delta)),obj.delta,[],[],[],[],[],[],@(delta) obj.potentialConstraints(delta));
        end

        function [c,ceq] = potentialConstraints(obj,delta)
            ceq = [];
            obj_new = obj;
            obj_new.g = delta;
            c = -delta(6);
            

            
            %obj_new.total('energy');
            
            vecs = obj_new.linkNormalVectors();
            for i=1:length(obj.links)
                c(end+1) = -vecs(3,i);
            end
        end
            
        
        function dirs = linkNormalVectors(obj)
            %grab all the normal vectors (tangent, wrong names)
            dirs = [];
            for i=1:length(obj.links)
                dirs(:,i) = obj.links{i}.normal_vector();
            end
        end
        
%         function L = cableLength(obj) %length of cable through centers
%             d = obj.delta;
%             L = norm(d(4:6));
%         end
%         
%         function obj = minimizeEnergyCableConstraint(obj,L)
%             %find the minimum energy configuration given a cable length
%             obj.g = fmincon(@(delta) obj.total('energy',delta),obj.delta,[],[],[],[],[],[],@(delta) obj.cableConstraints(delta,L));
%         end
%         
%         function [c,ceq] = cableConstraints(obj,delta,L)
%             obj_new = obj; %copy object
%             obj_new.g = delta; 
%             
%             c = -delta(6);
%             
%             %bound the z angle
%             
%             
%             ceq = obj_new.cableLength()-L;
%             
%             vecs = obj_new.linkNormalVectors();
%             for i=1:length(obj.links)
%                 c(end+1) = -vecs(3,i);
%             end
%         end
%         
        %the offset cable versions
        function L = cableLength(obj,r, theta)
            if nargin == 1
                r=0;
                theta=0;
            end
            r_bot = Rz(theta)*[r;0;0];
            r_top = Rz(theta-obj.init_delta(3))*[r;0;0];
            d = obj.g*r_top-r_bot;
            L = norm(d);
        end
        
        function obj = minimizeEnergyCableConstraint(obj,L,r,theta)
            %find the minimum energy configuration given a cable length
            if nargin == 2 %assume center points
                r = 0;
                theta = 0;
            end
            obj.g = fmincon(@(delta) obj.total('energy',delta),obj.delta,[],[],[],[],[],[],@(delta) obj.cableConstraints(delta,L,r,theta));
        end
        
        function [c,ceq] = cableConstraints(obj,delta,L,r,theta)
            obj_new = obj; %copy object
            obj_new.g = delta; 
            
            c = -delta(6);
            %the bending angles need to stay between -pi/2 and pi/2
%             %-pi/2<delta(1)
%             c(end+1) = -pi/2-delta(1);
%             c(end+1) = -pi/2-delta(2);
%             %delta(1)<pi/2
%             c(end+1) = delta(1)-pi/2;
%             c(end+1) = delta(2)-pi/2;
%             
%             %the torsion angle can't deviate from the initial angle by more
%             %than pi
%             %delta(3)-init_delta(3)<pi
%             %-pi<delta(3)-init_delta(3)
%             c(end+1) = delta(3)-obj.init_delta(3)-pi;
%             c(end+1) = obj.init_delta(3)-delta(3)-pi;
            %bound the z angle
            
            
            ceq = obj_new.cableLength(r,theta)-L;
            
            vecs = obj_new.linkNormalVectors();
            for i=1:length(obj.links)
                c(end+1) = -vecs(3,i);
            end
        end
        
    end

end