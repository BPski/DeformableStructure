function out = allCombinations(elems,n)
    %given the elements to go through and the length of the combinations
    %make all the combinations

    %allCombinations(elems,1) = elems
    %allCombinations(elems,n) = for c in allCombinations(elem,n-1) append
    %elems
    
    m = length(elems);
    
    if n == 1
        out = elems;
    else
        combs = allCombinations(elems,n-1);
        combs_new = {};
        r = length(combs);
        for i=1:r
            for j=1:m
                t = {elems{j},combs{i}};
                combs_new{(i-1)*m+j} = horzcat(t{:});
            end
        end
        out = combs_new;
    end
end

    
    
    