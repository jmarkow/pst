function TREE=pst_learn(F_MAT,ALPHABET,N,varargin)
%pst_learn takes the output of pst_build_trans_mat
%and build a probabilistic suffix tree per Ron, Singer and Tishby 1996
%"The Power of Amnesia: Learning Probabilistic
% Automata with Variable Memory Length" (algorithm LEARN-PSA)
%
%	TREE=pst_learn(F_MAT,ALPHABET,N,varargin)
%
%	F_MAT
%	cell array of frequency tables returned by pst_build_trans_mat
%
%	ALPHABET
%	string of symbols returned by pst_build_trans_mat
%
%	N
%	total entries per order returned by pst_build_trans_mat
%
% TREE
% array of structures, where the array idx indicates depth, with the following fields:
%
%	string
% indices the indicate the sequence represented by the node
% (e.g. if ABC are the first three characters in ALPHABET, the node w/ ABC
% would have a value of [1 2 3])
%
% parent
% contains the node idx (row 1) and depth (row 2) of the parent
% node
%
%	label
%	the nodel label
%
% internal
% 1 indicates internal node (node added to form a path to a parent node)
%
%	g_sigma_s
% probability of each symbol in ALPHABET conditioned on the node, each row corresponds to
% a given symbol in ALPHABET, and each column a given node index
%
%	p
% probability of the node occurring
%
%	f
% frequency of the node occurring
%
%	the following may be entered as parameter/value pairs:
%
%	L
%	maximum order (default: 7, ensure that you compute a sufficiently high order with
%	fa_build_matrix)
%
%	p_min
%	minimum occurrence probability to check as potential node in tree (default: .0073)
%
%	g_min
%	minimum transition probability (default: .185)
%
%	r
%	minimum divergence between parent and child node (default: 1.8)
%
%	alpha
%	smoothing parameter (default: 0)
%
% See also pst_build_trans_mat

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

% default per PLoS Comp Bio paper

L=7; % maximum order
p_min=.0073; % minimum occurrence probability
g_min=.01; % minimum probability at each node=1.75;
r=1.6; % minimum difference between parent and child node
alpha=17.5;
p_smoothing=0;

delimiter.start='['; % these are deprecated now
delimiter.finish=']';

for i=1:2:nparams
	switch lower(varargin{i})
		case 'p_min'
			p_min=varargin{i+1};
		case 'g_min'
			g_min=varargin{i+1};
		case 'r'
			r=varargin{i+1};
		case 'l'
			L=varargin{i+1};
		case 'alpha'
			alpha=varargin{i+1};
		case 'delimiter'
			delimiter=varargin{i+1};
		case 'p_smoothing'
			p_smoothing=varargin{i+1};
		otherwise
	end
end

% initialize sbar
% exclude symbols whose P<p_min

sbar={};
counter=1;
for i=1:length(ALPHABET)

	if (single(F_MAT{1}(i))./N(1))<p_min,
		continue;
	end

	% no point in looking at the zero order probabilities of
	% start and finish

	if ALPHABET(i)==delimiter.finish | ALPHABET(i)==delimiter.start
		continue;
	end

	sbar{counter}=ALPHABET(i);
	counter=counter+1;
end

% initialize the tree with the empty node

tbar(1).string{1}=[];
tbar(1).parent(:,1)=[0;0];
tbar(1).label{1}='epsilon';
tbar(1).internal(1)=0;

for i=2:L+1
	tbar(i).string={};
	tbar(i).parent=[];
	tbar(i).internal=[];
	tbar(i).label={};
end

while ~isempty(sbar)

	s_char=sbar{1};
	sbar(1)=[];

	%s_char
	s=[];

	for i=1:length(s_char)
		s(i)=findstr(s_char(i),ALPHABET);
	end

	%length(s)
	curr_depth=length(s)+1;

	% find p(sigma|s)>p_min

	f_vec=retrieve_f_sigma(F_MAT,s);

	if length(s)==1
		f_suf=retrieve_f_sigma(F_MAT,[]);
	else
		f_suf=retrieve_f_sigma(F_MAT,s(2:end));
	end

	p_sigma_s=f_vec(:)./(sum(f_vec(:))+eps);
	p_sigma_suf=f_suf(:)./(sum(f_suf(:))+eps);

	ratio=(p_sigma_s+eps)./(p_sigma_suf+eps);
	psize=(p_sigma_s>=(1+alpha)*g_min);

	ratio_test=(ratio>=r)|(ratio<=1/r);
	total=sum(ratio_test&psize);

	%total=sum((ratio>=r)&psize)+sum((ratio<=1/r)&psize);

	% if we have a non-empty total then s gets added to the tree

	if total>0
		tbar(curr_depth).string{end+1}=s;
		[node,depth]=findparent(s,tbar);
		tbar(curr_depth).parent(:,end+1)=[node;depth];
		tbar(curr_depth).label{end+1}=s_char;
		tbar(curr_depth).internal(end+1)=0;
	end

	if length(s)<L

		f_vec=retrieve_f_prime(F_MAT,s);
		p_sigmaprime_s=f_vec./N(curr_depth);
		addnodes=find(p_sigmaprime_s>=p_min);

		for j=1:length(addnodes)
			sbar{end+1}=[ ALPHABET(addnodes(j)) s_char ];
			%ALPHABET(addnodes(j))
		end

	end

end

% lastly we need to walk the tree and smooth probabilities
% may want to include frequency and p(ABCD)

% first, fill out the skeleton

% need to embed this in a while loop until the tree is filled

% or maybe not...

% check for nodes without clear paths

% fix paths in the tree so that there is a clear path
% to every node

tbar=fixpath(tbar);

% compute the smoothed transition probabilities

tbar=findgsig(tbar,F_MAT,g_min,N,p_smoothing);
TREE=tbar;

end


%%%%%%%%%%%%

function f=retrieve_f(f_mat,s)


% the following is an alternative way to implement the retrieve functions,
% eval is typically avoided for speed reasons

%left_edge= [sprintf('f=f_mat{%i}(',length(s))];
%middle='';
%for i=1:length(s)-1
%	middle=[ middle sprintf('s(%i),',i) ];
%end
%
%right_edge=[ sprintf('s(%i));',length(s)) ];
%
%
%[ left_edge middle right_edge ]

switch length(s)

	case 1
		f=f_mat{1}(s);
	case 2
		f=f_mat{2}(s(1),s(2));
	case 3
		f=f_mat{3}(s(1),s(2),s(3));
	case 4
		f=f_mat{4}(s(1),s(2),s(3),s(4));
	case 5
		f=f_mat{5}(s(1),s(2),s(3),s(4),s(5));
	case 6
		f=f_mat{6}(s(1),s(2),s(3),s(4),s(5),s(6));
	case 7
		f=f_mat{7}(s(1),s(2),s(3),s(4),s(5),s(6),s(7));
	case 8
		f=f_mat{8}(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8));
	case 9
		f=f_mat{9}(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9));
	case 10
		f=f_mat{10}(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9));
	otherwise

end

f=single(f);

end

%%%%%%%%%%%%%

function f_vec=retrieve_f_sigma(f_mat,s)


switch length(s)+1
	case 1
		f_vec=f_mat{1};
	case 2
		f_vec=f_mat{2}(s(1),:);
	case 3
		f_vec=squeeze(f_mat{3}(s(1),s(2),:));
	case 4
		f_vec=squeeze(f_mat{4}(s(1),s(2),s(3),:));
	case 5
		f_vec=squeeze(f_mat{5}(s(1),s(2),s(3),s(4),:));
	case 6
		f_vec=squeeze(f_mat{6}(s(1),s(2),s(3),s(4),s(5),:));
	case 7
		f_vec=squeeze(f_mat{7}(s(1),s(2),s(3),s(4),s(5),s(6),:));
	case 8
		f_vec=squeeze(f_mat{8}(s(1),s(2),s(3),s(4),s(5),s(6),s(7),:));
	case 9
		f_vec=squeeze(f_mat{9}(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),:));
	case 10
		f_vec=squeeze(f_mat{10}(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),:));
	otherwise

end

f_vec=single(f_vec);

end

function f_vec=retrieve_f_prime(f_mat,s)

switch	length(s)+1
	case 2
		f_vec=squeeze(f_mat{2}(:,s(1)));
	case 3
		f_vec=squeeze(f_mat{3}(:,s(1),s(2)));
	case 4
		f_vec=squeeze(f_mat{4}(:,s(1),s(2),s(3)));
	case 5
		f_vec=squeeze(f_mat{5}(:,s(1),s(2),s(3),s(4)));
	case 6
		f_vec=squeeze(f_mat{6}(:,s(1),s(2),s(3),s(4),s(5)));
	case 7
		f_vec=squeeze(f_mat{7}(:,s(1),s(2),s(3),s(4),s(5),s(6)));
	case 8
		f_vec=squeeze(f_mat{8}(:,s(1),s(2),s(3),s(4),s(5),s(6),s(7)));
	case 9
		f_vec=squeeze(f_mat{9}(:,s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8)));
	case 10
		f_vec=squeeze(f_mat{10}(:,s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9)));
	otherwise

end

f_vec=single(f_vec);

end

%%%%%%%%%%%%%%%

function tbar=add_node(tbar,sequence,label,parent,depth)

% find the deepest node that is a suffix of subsequence

% if the depth is 1, then simply add the node and all
% suffixes in the path

tbar(depth).string{end+1}=sequence;
tbar(depth).label{end+1}=label;
tbar(depth).parent(:,end+1)=parent;

end


%%%%%%%%%%%%%%%

function tbar=fixpath(tbar)
%
%
%
%

exitflag=0;
while ~exitflag

	changes=0;
	for i=3:length(tbar)
		for j=1:length(tbar(i).string)

			curr_string=tbar(i).string{j};

			% may have new parents

			[node,depth]=findparent(curr_string,tbar);
			parent_depth=tbar(i).parent(2,j);

			if depth>parent_depth
				tbar(i).parent(:,j)=[node;depth];
				parent_depth=tbar(i).parent(2,j);
			end

			% if we don't have a direct path, add
			% the nodes we need to get to the current
			% node

			if parent_depth<i-1
				tbar(i-1).string{end+1}=curr_string(2:end);
				[node,depth]=findparent(curr_string(2:end),tbar);
				tbar(i-1).parent(:,end+1)=[node;depth];
				tbar(i-1).label{end+1}=tbar(i).label{j}(2:end);
				tbar(i-1).internal(end+1)=1;
				changes=changes+1;
				tbar(i).parent(:,j)=[length(tbar(i-1).string);i-1];
			end

		end
	end

	if changes==0
		exitflag=1;
	end

end


end


%%%%%%%%%%%%%%

function tbar=findgsig(tbar,f_mat,g_min,N,p_smoothing)

for i=1:length(tbar)
	for j=1:length(tbar(i).string)

		f_vec=retrieve_f_sigma(f_mat,tbar(i).string{j});
		p_sigma_s=f_vec./(sum(f_vec)+eps);

		if ~isempty(tbar(i).string{j})
			f=retrieve_f(f_mat,tbar(i).string{j});
			p_s=f/N(length(tbar(i).string{j}));
		else
			f=0;
			p_s=1;
		end

		sigma_norm=length(p_sigma_s);

		g_sigma_s=p_sigma_s*(1-sigma_norm*g_min)+g_min;

		% replace with g_sigma_s to get smoothing


		if p_smoothing
			tbar(i).g_sigma_s(:,j)=g_sigma_s;
		else
			tbar(i).g_sigma_s(:,j)=p_sigma_s;
		end
		tbar(i).p(j)=p_s;
		tbar(i).f(j)=f;


	end

end

end

%%%%%%%%%%%%%%%%%%

function [node,depth]=findparent(sequence,tbar)

% given a sequence and tree, find the parent node

node=1;
depth=1;
seq_length=length(sequence);

if seq_length>1

	node=0;
	counter=0;
	for i=2:seq_length

		hits=[];

		for j=1:length(tbar(i).string)
			hits(j)=length(strfind(tbar(i).string{j},sequence(end-counter:end)));
		end

		if sum(hits>0)==1
			node=find(hits>0);
			depth=max(depth,i);
		end

		counter=counter+1;


	end
end

end
