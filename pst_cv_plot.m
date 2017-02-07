function pst_cv_plot(CVDATA,varargin)
% Plots the output of pst_cross_validate
%



sort_param='p_min';

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end


for i=1:2:nparams
	switch lower(varargin{i})
		case 'sort_param'
			sort_param=varargin{i+1};
		otherwise
	end
end

% average cv split with the same parameter value

% sort the data by the parameter of choice

sort_param
sort_data=CVDATA.(sort_param);
sort_data=sort_data(:);
param_values=unique(sort_data);

% get the mean value for each parameter value

mean_val=zeros(length(param_values),1);
std_val=zeros(size(mean_val));
n_val=zeros(size(mean_val));

for i=1:length(param_values)


	tmp=-CVDATA.test_logl(sort_data==param_values(i)); % convert to negative logl
	tmp2=-CVDATA.train_logl(sort_data==param_values(i));
	test.mean_val(i)=mean(tmp);
	test.std_val(i)=std(tmp);
	test.n_val(i)=length(tmp);
	train.mean_val(i)=mean(tmp2);
	train.std_val(i)=std(tmp2);
	test.n_val(i)=length(tmp2);

end

test.ci(1,:)=test.mean_val-(test.std_val);
test.ci(2,:)=test.mean_val+(test.std_val);
train.ci(1,:)=train.mean_val-(train.std_val);
train.ci(2,:)=train.mean_val+(train.std_val);

%[val,idx]=sort(sort_data(:));

figure();

plot(param_values,test.mean_val,'--ko','color','c','linewidth',3,'markersize',5)
box off
hold on;

for i=1:length(param_values)
	h2(1)=line([ param_values(i) param_values(i) ],[ test.mean_val(i) test.ci(2,i) ],'color','c','linewidth',3);
	h2(2)=line([ param_values(i) param_values(i) ],[ test.ci(1,i) test.mean_val(i) ],'color','c','linewidth',3);
end

plot(param_values,train.mean_val,'--ko','color','b','linewidth',3,'markersize',5)
box off
hold on;

for i=1:length(param_values)
	h2(1)=line([ param_values(i) param_values(i) ],[ train.mean_val(i) train.ci(2,i) ],'color','b','linewidth',3);
	h2(2)=line([ param_values(i) param_values(i) ],[ train.ci(1,i) train.mean_val(i) ],'color','b','linewidth',3);
end

set(gca,'xdir','rev','xscale','log');
set(gca,'FontSize',15,'FontName','Helvetica');

ylabel('CV Neg. Log Likelihood');
xlabel('Pmin');
