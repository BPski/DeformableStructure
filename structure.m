classdef structure
    %a structure is a collection of modules in series
    %this is for keeping everything well presented an bound together
    %between the modules
    
    properties
        modules %list of the modules
    end
    
    methods
        function obj = structure(modules)
            %constructor now just takes a list of predefined modules 
            %the structure class is then just keeping them organized
            
                        
            obj.modules = modules;

        end
        
        function obj = reset(obj)
            %reset each module
            for i=1:length(obj.modules)
                obj.modules(i) = obj.modules(i).reset();
            end
        end
        
        function deltas = deltas(obj)
            %collect all the deltas
            deltas = [];
            for i=1:length(obj.modules)
                deltas(:,i) = obj.modules(i).delta;
            end
        end
        
        function plot(obj)
            %plot each module seperately
            %have to increment the offset for each one
            hold on
            offset = configuration();
            
            for i=1:length(obj.modules)
                plot(obj.modules(i),offset);
                offset = offset*(obj.modules(i).g);
            end
        end
        
        function out = total(obj,name,deltas)
            %accumulate totals across modules
            if nargin == 2 %only given a name
                out = 0;
                for i=1:length(obj.modules)
                    out = out + obj.modules(i).total(name);
                end
            else
                out = 0;
                for i=1:length(obj.modules)
                    out = out + obj.modules(i).total(name,deltas(:,i));
                end
            end
        end
        
        function obj = minimize(obj,name)
            %minimize values that can be considered independent between
            %modules
            for i=1:length(obj.modules)
                obj.modules(i) = obj.modules(i).minimize(name);
            end
        end
        
        function obj = balanceWrench(obj,wrench)
            %apply a wrench to the top module and have it balance all the
            %way through
            for i=1:length(obj.modules)
                obj.modules(i) = obj.modules(i).balanceWrench(wrench);
                wrench = obj.modules(i).g.Ad'*wrench; %modify wrench
            end
        end
        
        function animateDeformationWrench(obj,wrench)
            %animate the deformation due to an applied wrench
            
            steps = 100;
            for i=1:steps
                w = wrench*i/steps;
                obj = obj.balanceWrench(w);
                obj.total('energy');
                clf;
                plot(obj);
                view(30,40)
                daspect([1;1;1])
                drawnow
            end
        end
                
        function L = cableLength(obj,r,theta) %if there is a cable running through a point what would the length be
            if nargin == 1
                r = 0;
                theta = 0;
            end
            L = 0;
            for i=1:length(obj.modules)
                L = L + obj.modules(i).cableLength(r,theta);
            end
        end
        
        function obj = minimizeEnergyCableDelta(obj,del_L,r,theta)
            %see which increment in del_L of the sub modules gives the
            %smallest energy change
            energy0 = [];
            energy = [];
            updated = [];
            for i=1:length(obj.modules)
                energy0(end+1) = obj.modules(i).total('energy');
                updated = [updated,obj.modules(i).minimizeEnergyCableConstraint(obj.modules(i).cableLength(r,theta)-del_L,r,theta)];
                energy(end+1) = updated(end).total('energy');
            end
            delta_energy = energy-energy0;
            [~,index] = min(delta_energy);
            obj_keep = updated(index);
            obj.modules(index).g = obj_keep.delta;
                
        end
        
        function obj = minimizeEnergyCableConstraint(obj,L,r,theta,file)
            %finding the minimum energy configuration given a certain cable
            %length
            %trying to do this incrementally like the previous version
            %at each step try a small increment to the different modules
            %and then see which change led to the smallest change in energy
            
            %kinda weird at this point the checking of the arguments
            if nargin ~=5
                file = 'cableDeformationAnimation.gif';
            end
            
            
            if nargin == 2
                r = 0;
                theta = 0;
            end
            
            steps = 500;
            grab = 5; %frames to grab
            L0 = obj.cableLength(r,theta);
            for i=1:steps
                [i,steps]
                del_L = (L0-L)/steps;
                obj = obj.minimizeEnergyCableDelta(del_L,r,theta);
                clf;
                plot(obj);
%                 view(30,40)
                xlim([-0.2,0.2])
                ylim([-0.2,0.2])
                zlim([0,0.4])
                daspect([1;1;1])
                drawnow
                obj.deltas;
                
                frame = getframe(gca);
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256);
                if i==1
                    imwrite(imind,cm,file,'gif', 'Loopcount',inf,'DelayTime',10/(steps/grab)); 
                elseif mod(i,grab) == 0
                   imwrite(imind,cm,file,'gif','WriteMode','append','DelayTime',10/(steps/grab)); 
                end 
            end
        end
        
        function obj = minimizeEnergyCableConstraintFunction(obj,f)
            %finding the minimum energy configuration given a certain cable
            %length
            %trying to do this incrementally like the previous version
            %at each step try a small increment to the different modules
            %and then see which change led to the smallest change in energy
            
            file = 'cableDeformationAnimationFunction.gif';
            
                        
            steps = 500;
            grab = 5;
            [L,r,theta] = f(steps,0);
            for i=1:steps
                [i,steps]
                [L_step,r_step,theta_step] = f(steps,i);
                del_L = L-L_step;
                obj = obj.minimizeEnergyCableDelta(del_L,r_step,theta_step);
                L = L_step;
                
                
                clf;
                plot(obj);
                %view(0,0)
                xlim([-0.1,0.1])
                ylim([-0.1,0.1])
                zlim([0,0.2])
                daspect([1;1;1])
                drawnow
                
                frame = getframe(gca);
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256);
                if i==1
                    imwrite(imind,cm,file,'gif', 'Loopcount',inf,'DelayTime',10/(steps/grab)); 
                elseif mod(i,grab) == 0
                   imwrite(imind,cm,file,'gif','WriteMode','append','DelayTime',10/(steps/grab)); 
                end 
            end
            
        end
        
        
    end
    
end
