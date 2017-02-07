function pst_export_to_cytoscape(TREE,ALPHABET,varargin)
%pst_export_to_cytoscape takes a PST computed by pst_learn
%and generates files suitable for use with Cytoscape (v.2.8.2 was used for
%all figures in the paper).
%
%	pst_export_to_cytoscape(TREE,ALPHABET,varargin)
%
%	TREE
%	structure array returned by fa_pst_learn
%
%	ALPHABET
%	mapping of phrase identities to rows/columns in frequency table, returned
%	by fa_pst_build_matrix
%
%	the following can be passed as parameter/value pairs:
%
%	output_dir
%	directory to store generated files (default: pwd)
%
%	filename
%	rootname for generated files (default: 'cytoscape_output_tree')
%
%The script will generate two files, one ending in .sif and the other .noa, the .sif
%file specifies the network connectivity and .noa information for each node to be loaded
%separately as a table
%
%

nparams=length(varargin);

output_dir=pwd;
filename='cytoscape_output_tree';
thresh=1e-5;

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'output_dir'
			output_dir=varargin{i+1};
		case 'filename'
			filename=varargin{i+1};
		case 'thresh'
			thresh=varargin{i+1};
		otherwise
	end
end

sif_id=fopen(fullfile(output_dir,[filename '.sif']),'w');
noa_gsigma_id=fopen(fullfile(output_dir,[filename '.noa']),'w');
cscript_id=fopen(fullfile(output_dir,[filename '_script.txt']),'w');
TREE(cellfun(@isempty,{TREE(:).string}))=[];

%eda_id=fopen(fullfile(output_dir,[filename '.eda']),'w');

%fprintf(eda_id,'EdgeTransitionProbability (class=Double)\n');

%fprintf(noa_gsigma_id,'GsigmaS (class=String)\n');
% climb down the tree from the root node

% list all possible walks, specifying connections from one level to the next
% on each line

% first we need to


for i=1:length(TREE)-1

	for j=1:length(TREE(i).label)

		source=TREE(i).label{j};

		if isempty(TREE(i+1).parent)
			continue;
		end

		target_idxs=find(TREE(i+1).parent(1,:)==j);
		target_label=TREE(i+1).label(target_idxs);

		fprintf(sif_id,'%s trans',source);

		for k=1:length(target_idxs)
			fprintf(sif_id,' %s',target_label{k});
		end

		fprintf(sif_id,'\n');

	end

end

% if user specifies remove syllables that are rarely transitioned to

if ~isempty(thresh)

	maxvec=zeros(length(ALPHABET),1);

	for i=2:length(TREE)
		tmp=max(TREE(i).g_sigma_s,[],2);
		maxvec=max(tmp,maxvec);
	end

	to_del=maxvec<thresh;

else
	to_del=[];
end

%ALPHABET(to_del)=[];
len=num2str(length(ALPHABET));
%cmap=im2uint8(colormap(['hsv(' len ')']));
cmap=im2uint8(distinguishable_colors(length(ALPHABET),'w'));

% convert the cmap to hex

hexcmap={};
for i=1:size(cmap,1)
	hexcmap{i}=[ '#' dec2hex(cmap(i,1),2) dec2hex(cmap(i,2),2) dec2hex(cmap(i,3),2) ];
end

%for i=1:length(TREE)
%	TREE(i).g_sigma_s(to_del,:)=[];
%end

%


fprintf(noa_gsigma_id,'ID\t');

if ischar(ALPHABET)
	newalph=cell(size(ALPHABET));
	for i=1:length(ALPHABET)
		newalph{i}=ALPHABET(i);
	end
	ALPHABET=newalph;
	clear newalph;
end

for i=1:length(ALPHABET)
	fprintf(noa_gsigma_id,'%s\t',ALPHABET{i});
end

fprintf(noa_gsigma_id,'Frequency\tLogFrequency\tDepth');

internalflag=0;
if isfield(TREE(1),'internal')
	internalflag=1;
	fprintf(noa_gsigma_id,'\tInternal\n');
else
	fprintf(noa_gsigma_id,'\n');
end

for i=1:length(TREE)
	for j=1:length(TREE(i).label)

		fprintf(cscript_id,'nodecharts pie nodelist="%s"',TREE(i).label{j});

		fprintf(noa_gsigma_id,'%s\t',TREE(i).label{j});


		labellist='';
		colorlist='';
		valuelist='';

		for k=1:length(TREE(i).g_sigma_s(:,j))
			fprintf(noa_gsigma_id,'%0.2f\t',TREE(i).g_sigma_s(k,j));

			if TREE(i).g_sigma_s(k,j)>thresh
				labellist=[ labellist sprintf('%s,',ALPHABET{k}) ];
				colorlist=[ colorlist sprintf('%s,',hexcmap{k}) ];
				valuelist=[ valuelist sprintf('%g,',TREE(i).g_sigma_s(k,j)) ];
			end

		end

		%fprintf(cscript_id,' attributelist="%s"',labellist(1:end-1));

		if length(labellist(1:end-1))==2
			labellist(1:end-1)
			labellist=[ labellist 'null,']
			valuelist=[ valuelist '0,']
			colorlist=[ colorlist sprintf('%s,',hexcmap{1}) ];
		end


		fprintf(cscript_id,' labellist="%s"',labellist(1:end-1));
		fprintf(cscript_id,' valuelist="%s"',valuelist(1:end-1));

		fprintf(cscript_id,' colorlist="%s"\n',colorlist(1:end-1));

		fprintf(noa_gsigma_id,'%g\t',TREE(i).f(j));
		fprintf(noa_gsigma_id,'%g\t',log(TREE(i).f(j)));


		fprintf(noa_gsigma_id,'%g',i-1);

		if internalflag
			fprintf(noa_gsigma_id,'\t%i',TREE(i).internal(j));
		end

		fprintf(noa_gsigma_id,'\n');

		%fprintf(noa_gsigma_id,'\n');

	end
end




fclose(noa_gsigma_id);
fclose(sif_id);
fclose(cscript_id);


% create a simple Command Tool script to format everything, including the pie charts
