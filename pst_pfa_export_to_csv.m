function pst_pfa_export_to_csv(PFA,varargin)
%pst_pfa_export_to_csv will convert a PFA to Excel format
%
%	pst_pfa_export_to_csv(PFA,ALPHABET,varargin)
%
%	PFA
%
%	filename
%	name of file to write to (default: transitions.csv)
%

%
% See also fa_pst_learn

if nargin<2, error('Need transition matrix and alphabet to continue'); end

nparams=length(varargin);

output_dir=pwd;
filename='transitions';

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

csv_id=fopen(fullfile(output_dir,[filename '.csv']),'w');

n_states=length(PFA)-1;

for i=2:length(PFA)
	fprintf(csv_id,',%s',PFA(i).label);
end

for i=2:length(PFA)

	source=PFA(i).label;
	fprintf(csv_id,'\n%s',source);

	trans=zeros(1,n_states);
	trans(PFA(i).arcs-1)=PFA(i).arcs_p;

	for j=1:n_states
		fprintf(csv_id,',%.5f',trans(j));
	end
end

fclose(csv_id);
