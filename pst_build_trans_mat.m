function [F_MAT ALPHABET N PI]=pst_build_trans_mat(DATA,ORDER)
%pst_build_trans_mat takes the cell array of strings DATA
%and generates a series of transition matrices up to 9th order (10 dimensional matrix). The
%matrices are all uint16, so change the data type if any cell in your frequency table
%might exceed 65535
%
%	[F_MAT ALPHABET N PI]=pst_build_trans_mat(DATA,ORDER)
%
%	DATA
%	cell array, where each element contains a string of alphanumeric characters indicating
% different symbols or syllables.  Elements are interpreted as separate trials.
%
%	ORDER
%	maximum order frequency table to compute (default: 5)
%
%	ALPHABET
%	string of all unique symbols found in the data (i.e. the corpus)
%
%	N
%	vector of total entries for each order
%
%	PI
%	starting distribution
%

if nargin<2, ORDER=5; end
if ORDER<1, error('Order must be at least 1!'); end
if ORDER>9, error('Maximum order is 8, user selected %g',ORDER); end

[sequence ALPHABET]=pst_sequence_gen(DATA);
%ALPHABET=[ALPHABET ']']; % add our end delimiter
ncat=length(ALPHABET); % add the end as an additional category

F_MAT{1}=zeros(1,ncat,'uint16');
F_MAT{2}=zeros(ncat,ncat,'uint16');
PI=zeros(size(F_MAT{1}));

if ORDER>1, F_MAT{3}=zeros(ncat,ncat,ncat,'uint16'); end
if ORDER>2, F_MAT{4}=zeros(ncat,ncat,ncat,ncat,'uint16'); end
if ORDER>3, F_MAT{5}=zeros(ncat,ncat,ncat,ncat,ncat,'uint16'); end
if ORDER>4, F_MAT{6}=zeros(ncat,ncat,ncat,ncat,ncat,ncat,'uint16'); end
if ORDER>5, F_MAT{7}=zeros(ncat,ncat,ncat,ncat,ncat,ncat,ncat,'uint16'); end
if ORDER>6, F_MAT{8}=zeros(ncat,ncat,ncat,ncat,ncat,ncat,ncat,ncat,'uint16'); end
if ORDER>7, F_MAT{9}=zeros(ncat,ncat,ncat,ncat,ncat,ncat,ncat,ncat,ncat,'uint16'); end
if ORDER>8, F_MAT{10}=zeros(ncat,ncat,ncat,ncat,ncat,ncat,ncat,ncat,ncat,ncat,'uint16'); end


% reached largest variable size, could use memmapfile

for i=1:length(F_MAT)
	N(i)=0;
end

for i=1:length(DATA)

	%song=[ DATA{i} ']' ]; % add the end delimiter
	song=[ DATA{i} ];

	if length(song)>1

		A=find(song(1)==ALPHABET);
		PI(A)=PI(A)+1;

	end

	for j=1:length(song)

		A=find(song(j)==ALPHABET);
		F_MAT{1}(A)=F_MAT{1}(A)+1;
		N(1)=N(1)+1;

	end

	for j=1:length(song)-1

		A=find(song(j)==ALPHABET);
		B=find(song(j+1)==ALPHABET);
		F_MAT{2}(A,B)=F_MAT{2}(A,B)+1;
		N(2)=N(2)+1;

	end

	if ORDER>1

		for j=1:length(song)-2

			A=find(song(j)==ALPHABET);
			B=find(song(j+1)==ALPHABET);
			C=find(song(j+2)==ALPHABET);

			F_MAT{3}(A,B,C)=F_MAT{3}(A,B,C)+1;
			N(3)=N(3)+1;

		end

	end

	if ORDER>2

		for j=1:length(song)-3

			A=find(song(j)==ALPHABET);
			B=find(song(j+1)==ALPHABET);
			C=find(song(j+2)==ALPHABET);
			D=find(song(j+3)==ALPHABET);

			F_MAT{4}(A,B,C,D)=F_MAT{4}(A,B,C,D)+1;
			N(4)=N(4)+1;

		end

	end

	if ORDER>3
		for j=1:length(song)-4

			A=find(song(j)==ALPHABET);
			B=find(song(j+1)==ALPHABET);
			C=find(song(j+2)==ALPHABET);
			D=find(song(j+3)==ALPHABET);
			E=find(song(j+4)==ALPHABET);

			F_MAT{5}(A,B,C,D,E)=F_MAT{5}(A,B,C,D,E)+1;
			N(5)=N(5)+1;

		end
	end

	if ORDER>4
		for j=1:length(song)-5

			A=find(song(j)==ALPHABET);
			B=find(song(j+1)==ALPHABET);
			C=find(song(j+2)==ALPHABET);
			D=find(song(j+3)==ALPHABET);
			E=find(song(j+4)==ALPHABET);
			F=find(song(j+5)==ALPHABET);

			F_MAT{6}(A,B,C,D,E,F)=F_MAT{6}(A,B,C,D,E,F)+1;
			N(6)=N(6)+1;

		end
	end

	if ORDER>5
		for j=1:length(song)-6

			A=find(song(j)==ALPHABET);
			B=find(song(j+1)==ALPHABET);
			C=find(song(j+2)==ALPHABET);
			D=find(song(j+3)==ALPHABET);
			E=find(song(j+4)==ALPHABET);
			F=find(song(j+5)==ALPHABET);
			G=find(song(j+6)==ALPHABET);

			F_MAT{7}(A,B,C,D,E,F,G)=F_MAT{7}(A,B,C,D,E,F,G)+1;
			N(7)=N(7)+1;

		end
	end


	if ORDER>6
		for j=1:length(song)-7

			A=find(song(j)==ALPHABET);
			B=find(song(j+1)==ALPHABET);
			C=find(song(j+2)==ALPHABET);
			D=find(song(j+3)==ALPHABET);
			E=find(song(j+4)==ALPHABET);
			F=find(song(j+5)==ALPHABET);
			G=find(song(j+6)==ALPHABET);
			H=find(song(j+7)==ALPHABET);

			F_MAT{8}(A,B,C,D,E,F,G,H)=F_MAT{8}(A,B,C,D,E,F,G,H)+1;
			N(8)=N(8)+1;

		end
	end

	if ORDER>7
		for j=1:length(song)-8

			A=find(song(j)==ALPHABET);
			B=find(song(j+1)==ALPHABET);
			C=find(song(j+2)==ALPHABET);
			D=find(song(j+3)==ALPHABET);
			E=find(song(j+4)==ALPHABET);
			F=find(song(j+5)==ALPHABET);
			G=find(song(j+6)==ALPHABET);
			H=find(song(j+7)==ALPHABET);
			I=find(song(j+8)==ALPHABET);

			F_MAT{9}(A,B,C,D,E,F,G,H,I)=F_MAT{9}(A,B,C,D,E,F,G,H,I)+1;
			N(9)=N(9)+1;

		end

	end


	if ORDER>8
		for j=1:length(song)-9

			A=find(song(j)==ALPHABET);
			B=find(song(j+1)==ALPHABET);
			C=find(song(j+2)==ALPHABET);
			D=find(song(j+3)==ALPHABET);
			E=find(song(j+4)==ALPHABET);
			F=find(song(j+5)==ALPHABET);
			G=find(song(j+6)==ALPHABET);
			H=find(song(j+7)==ALPHABET);
			I=find(song(j+8)==ALPHABET);
			J=find(song(j+9)==ALPHABET);

			F_MAT{10}(A,B,C,D,E,F,G,H,I,J)=F_MAT{10}(A,B,C,D,E,F,G,H,I,J)+1;
			N(10)=N(10)+1;

		end
	end

end
