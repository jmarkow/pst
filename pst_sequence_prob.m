function [SEQP LOGL]=pst_sequence_prob(TREE,ALPHABET,BOUT,PI,varargin)
%pst_sequence_prob computes the probability of state sequences
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

if nargin<3
	PI=[];
end

internal=0;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'delimiter'
			delimiter=varargin{i+1};
		case 'internal'
			internal=varargin{i+1};
		otherwise
	end
end

% initialize sbar
% exclude symbols whose P<p_min

% take in the sequence and use the starting state to initialize
% the song, use findparent to locate the parent node and use
% g_sigma_s as the transmission probability

% step through each sequence and calculate the probability of the
% next symbol

SEQP=zeros(length(BOUT),1);
tmp=[];

for i=1:length(BOUT)

	seqtmp=BOUT{i};
	seq=[];
	breakflag=0;

	for j=1:length(seqtmp)

		if isempty(find(seqtmp(j)==ALPHABET))
			warning('Character not in training alphabet, skipping bout');
			SEQP(i)=0;
			breakflag=1;
			break;
		end

		seq(j)=find(seqtmp(j)==ALPHABET);
	end

	if breakflag
		continue;
	end

	seqp=zeros(length(seqtmp)-1,1);

	for j=2:length(seqtmp)

		[node,depth]=findparent(seq(1:j-1),TREE,ALPHABET,internal);
		seqp(j-1)=TREE(depth).g_sigma_s(find(ALPHABET==seqtmp(j)),node);

	end

	SEQP(i)=prod(seqp);
	tmp=[tmp;seqp];

	if ~isempty(PI)
		SEQP(i)=SEQP(i)*PI(seq(1));
		tmp=[tmp;PI(seq(1))];
	end

end

logltmp=log2(tmp+eps);
LOGL=sum(logltmp)./(length(logltmp));

end

%%%%%%%%%%%%%%%%%%

function [node,depth]=findparent(sequence,tbar,ALPHABET,internal)

% given a sequence and tree, find the parent node

node=1;
depth=1;
seq_length=length(sequence);

checklength=min(seq_length,length(tbar)-1);

node=1;
counter=0;

for i=1:checklength

	hits=[];

	for j=1:length(tbar(i+1).string)

		if internal
			if tbar(i+1).internal(j)
				hits(j)=0;
				continue;
			end
		end

		hits(j)=length(strfind(tbar(i+1).string{j},sequence(end-counter:end)));
	end

	if sum(hits>0)==1
		node=find(hits>0);
		depth=max(depth,i+1);
	end

	counter=counter+1;

end

end
