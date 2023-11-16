function A = mergestructs(x,y)
A=cell2struct([struct2cell(x);struct2cell(y)], ...
              [fieldnames(x);fieldnames(y)]);
end
