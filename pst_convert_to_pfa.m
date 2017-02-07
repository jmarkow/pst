function PFA=pst_convert_to_pfa(TREE,ALPHABET)
% Converts a PST (probabilistic suffix tree) to a probabilistic finite automaton
%

% convert a PST to a PFA using the following algorithm:
%
%	1) include all leaves and internal nodes
%	2) starting with the root, flow downward (in terms of order)
%	3) arcs represent generating a given symbol, connecting each state
%          to the state representing the longest suffix of the next generated symbol
%	4) after visiting all nodes, we're finished


% from the tree, start from the bottom and work our way up

counter=1;

for i=1:length(TREE)
	for j=1:length(TREE(i).label)
		%PFA(counter).states=TREE(i).string{j};

		if TREE(i).internal(j)
			continue;
		end

		PFA(counter).label=TREE(i).label{j};
		PFA(counter).trans=TREE(i).g_sigma_s(:,j)';
		PFA(counter).order=i;
		PFA(counter).recurrent=1;
		counter=counter+1;
	end
end

% first we need to ensure that every leaf and its prefix are included...

% add the necessary arcs

% node order
% add the unused one symbol NODES FIRST!!!

for i=1:length(PFA)
	ORDER(i)=length(PFA(i).label);
end

%

exitflag=0;
for i=1:length(PFA)

	conns=find(PFA(i).trans>0);

	% for each connection, get the sequence and PFA

	currlabel=PFA(i).label;

	if strcmp(currlabel,'epsilon')
		currlabel='';
	end

	for j=1:length(conns)

		node=get_suffix([ currlabel ALPHABET(conns(j)) ],PFA);
		node

		% add any new nodes if necessary, should be one symbol since no other
		% suffix exists

		% get the edge of the first order

		if isnan(node)

			order_edge=max(find(ORDER==1));
			insertion_point=order_edge+1;

			NEWPFA=PFA;
			%NEWPFA(insertion_point).states=conns(j);
			NEWPFA(insertion_point).label=ALPHABET(conns(j));
			NEWPFA(insertion_point).trans=PFA(1).trans';
			NEWPFA(insertion_point+1:end+1)=PFA(insertion_point:end);

			PFA=NEWPFA;

			NEWORDER=ORDER;
			NEWORDER(insertion_point)=1;
			NEWORDER(insertion_point+1:end+1)=ORDER(insertion_point:end);

			ORDER=NEWORDER;
			clear NEWORDER;
			clear NEWPFA;

			node=insertion_point;

		end

		PFA(i).arcs(j)=node;
		PFA(i).arcs_p(j)=PFA(i).trans(conns(j));
		PFA(i).arcs_states(j)=ALPHABET(conns(j));
	end


end

% need to add all prefixes here or we'll be caught with orphaned nodes

exitflag=0;

% add all prefixes that don't already exist as states, use the longest suffix for transition probabilities

while ~exitflag

	exitflag=1;
	for i=1:length(PFA)
		if PFA(i).order>2

			% get list of prefixes

			brkflag=0;
			match=[];

			for j=length(PFA(i).label)-1:-1:1

				currprefix=PFA(i).label(1:j)
				match=[];

				for k=1:length(PFA)
					if strcmp(currprefix,PFA(k).label)
						match=k
						break;
					end
				end

				% if match is empty append a new prefix

				if isempty(match)

					newnode=length(PFA)+1;
					node=get_suffix(currprefix,PFA);

					PFA(newnode).label=currprefix

					if isnan(node)
						node=1;
					end

					node

					PFA(newnode).trans=PFA(node).trans;
					PFA(newnode).order=length(currprefix);


					conns=find(PFA(newnode).trans>0);

					for k=1:length(conns)
						node=get_suffix([ currlabel ALPHABET(conns(k)) ],PFA);
						PFA(newnode).arcs(k)=node;
						PFA(newnode).arcs_p(k)=PFA(newnode).trans(conns(k));
						PFA(newnode).arcs_states(k)=ALPHABET(conns(k));
					end

					exitflag=0;
				end
			end


		end

	end
end

for i=1:length(PFA)

	conns=find(PFA(i).trans>0);

	% for each connection, get the sequence and PFA

	currlabel=PFA(i).label;

	if strcmp(currlabel,'epsilon')
		currlabel='';
	end

	for j=1:length(conns)

		node=get_suffix([ currlabel ALPHABET(conns(j)) ],PFA);
		node

		% add any new nodes if necessary, should be one symbol since no other
		% suffix exists

		% get the edge of the first order

		if isnan(node)
			error('Broken connection');
		end

		PFA(i).arcs(j)=node;
		PFA(i).arcs_p(j)=PFA(i).trans(conns(j));
		PFA(i).arcs_states(j)=ALPHABET(conns(j));
	end


end

end

function node=get_suffix(SEQUENCE,PFA)
%
% for a given tree, find the longest suffix of a given sequence
%

% remove the prefix until we find our match (i.e. longest suffix match)

exitflag=0;
curr=SEQUENCE;
while ~exitflag
	if isempty(curr)
		node=NaN;
		exitflag=1;
		break;
	end

	for j=length(PFA):-1:1

		match=strcmp(PFA(j).label,curr);

		if match
			node=j;
			exitflag=1;
		end
	end

	curr=curr(2:end);

end


end
