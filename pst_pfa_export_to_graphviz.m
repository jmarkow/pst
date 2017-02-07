function pst_pfa_export_to_graphviz(PFA,varargin)
%pst_pfa_export_to_graphviz exports a PFA to a file suitable for visualization using
%graphviz
%

nparams=length(varargin);

output_dir=pwd;
filename='pfa_graphviz_export';

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'output_dir'
			output_dir=varargin{i+1};
		case 'filename'
			filename=varargin{i+1};
		otherwise
	end
end

dot_id=fopen(fullfile(output_dir,[filename '.dot']),'w');

%fprintf(noa_gsigma_id,'GsigmaS (class=String)\n');
% climb down the tree from the root node

% list all possible walks, specifying connections from one level to the next
% on each line

% first we need to

fprintf(dot_id,'digraph pfa {\n');
fprintf(dot_id,'\tnode [color=lightblue2, style=filled];\n')
fprintf(dot_id,'\trankdir=UD;\nsize="6.5,8.5"'); % changed rankdir from LR to UD

for i=2:length(PFA)

	source=PFA(i).label;

	for j=1:length(PFA(i).arcs)
		target=PFA(i).arcs(j);

		if PFA(i).arcs_p(j)<0.05
			continue;

		end

		% gray 52 was used for paper, and 5 was added to make lines visible for talks

		if PFA(i).arcs_p(j)<.2
			fprintf(dot_id,'\t%s -> %s [penwidth=%g,weight=2,color=gray60]',source,PFA(target).label,5+PFA(i).arcs_p(j));
		else
			fprintf(dot_id,'\t%s -> %s [penwidth=%g,weight=2,color=black]',source,PFA(target).label,PFA(i).arcs_p(j)*40);
		end

		fprintf(dot_id,'\n');

	end
end

for i=2:length(PFA)
	source=PFA(i).label;
	fprintf(dot_id,'\t%s [fontsize=55,ranksep=.5,nodesep=.5];\n',source);
end

%for i=2:length(PFA)
%
%	source=PFA(i).label;
%
%	if length(PFA(i).label)==1
%		fprintf(dot_id,'\t%s [rank=source]\n',source);
%	elseif length(PFA(i).label)>2
%		fprintf(dot_id,'\t%s [rank=sink]\n',source);
%	end
%end

fprintf(dot_id,'}')
fclose(dot_id);
% table for cytoscape 3.0
