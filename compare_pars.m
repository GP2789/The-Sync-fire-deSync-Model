function [ diff_par ] = compare_pars( filename_1, filename_2 )


load(filename_1)
par_1 = par;
load(filename_2)
par_2 = par;

F1 = fieldnames(par_1);
F2 = fieldnames(par_2);

M1 = F1(~contains(F1,F2));
M2 = F2(~contains(F2,F1));

diff_par = cell(1,3); c = 1;

for i = 1:length(F1)
    is_diff = false;
    F2_i = find(strcmp(F2,F1{i}),1,'first');
    if(~isempty(F2_i))
        if( ischar( par_1.(F1{i}) ) )
            if(~strcmp(par_1.(F1{i}), par_2.(F2{F2_i})))
                is_diff = true;
            end
        elseif( ~iscell( par_1.(F1{i}) ) )
            if(par_1.(F1{i}) ~= par_2.(F2{F2_i}))
                is_diff = true;
            end
        end
        if(is_diff)
            diff_par{c,1} = F1{i};
            diff_par{c,2} = par_1.(F1{i});
            diff_par{c,3} = par_2.(F1{i}); c = c + 1;
        end
    end
end

end

