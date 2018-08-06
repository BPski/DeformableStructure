function out = searchForWalking()
    %try searching for a configuration that gives a walking trajectory
    %hard to quantify motion, so just record all the configurations and
    %watch them to see if they are correct
    
    %can change between triangle, square, etc
    %for base only look one, two, etc soft joints
    %for top have to try each joint as soft and then proper pairs
    
    r = 5e-2;
    l = Link(10e-2,1e9*pi*3e-3^2/4,[0;0;256]/256);
    l_soft = Link(10e-2,1e6*pi*3e-3^2/4,[0;256;0]/256);
    SL = {l_soft,l};
    SR = {l,l_soft};
    H = {l,l};
    for n=3:4 %shape of polygon
        theta = 2*pi/n;
        rs = ones(n,1)*r;
        thetas = (0:n-1)*theta;
        
        base_combinations = allCombinations({SL;H},n-2);
        
        top_left_combinations = allCombinations({SL;H},n);
        top_right_combinations = allCombinations({SR;H},n);
        top_combinations = {top_left_combinations{2:end-1},top_right_combinations{2:end-1}};
        
        for i=1:length(base_combinations)
            comb = base_combinations(i);
            base = {SL,comb{:},H};
            base = horzcat(base{:});
            m1 = module(rs,rs,thetas,thetas,base);

            for j=1:length(top_combinations)
                top = top_combinations{j};
                m2 = module(rs,rs,thetas,thetas,top);
                s = structure([m1,m2]);
                
                %now that the structure is setup run the trajectories and
                %save them
                %store the combination index for the base and the top to
                %recall it later
                file = ['animations/','walkAnimation',num2str(n),'_',num2str(i),'_',num2str(j),'.gif']
                s.minimizeEnergyCableConstraint(s.cableLength(0,0)-0.05,0.0,0,file);
            end
        end
    end
end