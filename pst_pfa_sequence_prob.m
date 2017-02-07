function [SEQP LOGL]=pst_pfa_sequence_prob(PFA,ALPHABET,BOUT,varargin)
%takes a pfa and computes the log likelihood of the sequence

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'delimiter'
			delimiter=varargin{i+1};
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

	% get the starting state from the first symbol and
	% walk the graph to get our likelihoods

	for j=1:length(PFA)
		if strcmp(PFA(j).label,seqtmp(1))
			startnode=j;
			break;
		end
	end

	% from the startnode walk along

	currnode=PFA(startnode);

	% add log(p) along each edge

	seqp=[];

	for j=2:length(seqtmp)

		idx=strfind(currnode.arcs_states,seqtmp(j));
		seqp(j-1,1)=currnode.arcs_p(idx);
		nextnode=currnode.arcs(idx);
		currnode=PFA(nextnode);

	end

	tmp=[tmp;seqp];
	SEQP(i)=sum(log2(seqp));

end

LOGL=sum(log2(tmp+eps))./(length(tmp));

end
