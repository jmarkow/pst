function cv_store=pst_cross_validate(BOUT,varargin)
%runs a grid-search over plausible parameters using cross-validated
%held-out likelihood
%

%g_min=[ .0001 .001 .01 .025 .05 ];
%p_min=[ 1e-4 1e-3 5e-3 7e-3 1e-2 .025 .05 ];

r=[1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];

p_min=[ repmat(.007,[1 length(r)]) ];
%r=[ repmat(1.6,[1 length(p_min)]) ];
g_min=[ repmat(.01,[1 length(p_min)]) ];
alpha=[ repmat(17.5, [1 length(p_min)]) ];
L=[ repmat(7,[1 length(p_min)]) ];

params={ 'p_min','r','g_min','alpha','L' };

repetitions=5;
ncv=10; % number of folds, set to nsongs for leave-one-out


nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end


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
		case 'repetitions'
			repetitions=varargin{i+1};
		case 'ncv'
			ncv=varargin{i+1};
		otherwise
	end
end

nsongs=length(BOUT);

% break dataset into folds
% form the partitions
% what if we simply vary p_min?

counter=1;

for i=1:length(params)
	cv_store.(params{i})=zeros(repetitions,length(p_min),ncv);
end

cv_store.train_logl=zeros(repetitions,length(p_min),ncv);
cv_store.test_logl=zeros(size(cv_store.train_logl));
cv_store.tree={};

for i=1:repetitions

	disp(['Rep ' num2str(i)]);

	[testsamples,trainsamples]=cv_parts(nsongs,ncv,1); % get random CV split

	for j=1:length(p_min)

		for k=1:ncv

			disp(['CV ' num2str(k)]);

			% training samples

			trainbouts=BOUT(trainsamples{k});
			testbouts=BOUT(testsamples{k});

			[f_mat alphabet n pi_dist]=pst_build_trans_mat(trainbouts,7);
			tree=pst_learn(f_mat,alphabet,n,'g_min',g_min(j),'p_min',p_min(j),'r',r(j),'alpha',alpha(j),'L',L(j));

			pi_dist=double(pi_dist+1)./sum(pi_dist+1);

			% uncomment the following line to not use the starting distribution
			pi_dist=[];

			% compute the test and training likelihood

			[qx testlogl]=pst_sequence_prob(tree,alphabet,testbouts,pi_dist);
			[trainqx trainlogl]=pst_sequence_prob(tree,alphabet,trainbouts,pi_dist);

			cv_store.tree{i,j,k}=tree;
			cv_store.train_logl(i,j,k)=trainlogl;
			cv_store.test_logl(i,j,k)=testlogl;

			cv_store.p_min(i,j,k)=p_min(j);
			cv_store.g_min(i,j,k)=g_min(j);
			cv_store.r(i,j,k)=r(j);
			cv_store.alpha(i,j,k)=alpha(j);
			cv_store.L(i,j,k)=L(j);

		end

	end
end

end

function [TEST,TRAIN]=cv_parts(LEN,NCV,RND)
%
% pulls out cv partitions
%
%

partsize=floor(LEN/NCV)-1;

if RND
	pool=randperm(LEN);
else
	pool=1:LEN;
end

TEST={};
TRAIN={};

counter=1;

% pick new splits for each repetition

for ii=1:NCV
	TEST{ii}=counter:counter+partsize;
	TEST{ii}(TEST{ii}>LEN)=[];
	TEST{ii}=pool(TEST{ii});
	TRAIN{ii}=setdiff(pool,TEST{ii});
	counter=counter+1+partsize;
end


end
