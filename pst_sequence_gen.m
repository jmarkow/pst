function [SEQUENCE ALPHABET]=pst_sequence_gen(DATA,DELIM,RANDOMIZE);
%pst_sequence_gen
%
%
%
%
%
%


SEQUENCE=[];

if nargin<3, RANDOMIZE=0; end
if nargin<2, DELIM=''; end

for i=1:length(DATA)
	curr_string=DATA{i};

	if RANDOMIZE
		curr_string=curr_string(randperm(length(curr_string)));
	end

	SEQUENCE=[ SEQUENCE DELIM curr_string DELIM ];
end

ALPHABET=sort(unique(SEQUENCE));
