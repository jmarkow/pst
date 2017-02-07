function NEWTREE=pst_translate(TREE,SYLLABLES,SYLLABLES_MAP)
%
%
%
%
%
%
%

NEWTREE=TREE;

for i=2:length(TREE)
	for j=1:length(TREE(i).label)

		NEWTREE(i).label{j}='';

		for k=1:length(TREE(i).label{j})
			NEWTREE(i).label{j}=[ NEWTREE(i).label{j} SYLLABLES{strcmp(SYLLABLES_MAP,TREE(i).label{j}(k))} ];
		end
	end
end
